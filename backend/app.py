import os
import hmac
import json
import base64
import hashlib
import sqlite3
import uuid
import shutil
from datetime import datetime, timedelta
from typing import Optional, List
from fastapi import FastAPI, Depends, HTTPException, status, Request, Response, UploadFile, File, BackgroundTasks
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import httpx
from fastapi.staticfiles import StaticFiles

DATA_DIR = os.getenv("ENZYME_DATA_DIR", os.path.join(os.getcwd(), "data"))
DB_PATH = os.path.join(DATA_DIR, "meta.sqlite")
JWT_SECRET = os.getenv("ENZYME_JWT_SECRET", "dev_secret")
STATIC_DIR = os.getenv("ENZYME_STATIC_DIR")

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

if STATIC_DIR and os.path.isdir(STATIC_DIR):
    app.mount("/ui", StaticFiles(directory=STATIC_DIR, html=True), name="ui")

def ensure_dirs():
    os.makedirs(DATA_DIR, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS users(id TEXT PRIMARY KEY, email TEXT UNIQUE, password_salt BLOB, password_hash BLOB, created_at TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS jobs(id TEXT PRIMARY KEY, uid TEXT, type TEXT, status TEXT, created_at TEXT, inputs_json TEXT, outputs_json TEXT, metrics_json TEXT)")
    conn.commit()
    try:
        cur.execute("ALTER TABLE users ADD COLUMN password_plain TEXT")
        conn.commit()
    except Exception:
        pass
    try:
        # add uid to jobs if missing
        cur.execute("SELECT uid FROM jobs LIMIT 1")
    except Exception:
        try:
            cur.execute("ALTER TABLE jobs ADD COLUMN uid TEXT")
            conn.commit()
        except Exception:
            pass
    try:
        cur.execute("SELECT id FROM users WHERE email=?", ("admin",))
        row = cur.fetchone()
        if not row:
            uid = str(uuid.uuid4())
            salt = os.urandom(16)
            pwd = pbkdf2_hash("admin", salt)
            cur.execute("INSERT INTO users(id,email,password_salt,password_hash,password_plain,created_at) VALUES(?,?,?,?,?,?)", (
                uid, "admin", salt, pwd, "admin", datetime.utcnow().isoformat()
            ))
            conn.commit()
        else:
            cur.execute("UPDATE users SET password_plain=? WHERE email=?", ("admin", "admin"))
            conn.commit()
    except Exception:
        pass
    conn.close()

ensure_dirs()

class RegisterInput(BaseModel):
    email: str
    password: str

class LoginInput(BaseModel):
    email: str
    password: str

class JobCreate(BaseModel):
    type: str
    inputs: dict

class FileArtifactInput(BaseModel):
    job_id: str
    artifact: str

def pbkdf2_hash(password: str, salt: bytes) -> bytes:
    return hashlib.pbkdf2_hmac("sha256", password.encode(), salt, 100_000)

def sign_token(payload: dict, secret: str) -> str:
    body = json.dumps(payload, separators=(",", ":")).encode()
    sig = hmac.new(secret.encode(), body, hashlib.sha256).digest()
    return base64.urlsafe_b64encode(body).decode() + "." + base64.urlsafe_b64encode(sig).decode()

def verify_token(token: str, secret: str) -> Optional[dict]:
    try:
        body_b64, sig_b64 = token.split(".")
        body = base64.urlsafe_b64decode(body_b64.encode())
        sig = base64.urlsafe_b64decode(sig_b64.encode())
        expected = hmac.new(secret.encode(), body, hashlib.sha256).digest()
        if not hmac.compare_digest(sig, expected):
            return None
        return json.loads(body.decode())
    except Exception:
        return None

def get_conn():
    return sqlite3.connect(DB_PATH)

def auth_user(request: Request) -> Optional[str]:
    token = request.cookies.get("token")
    if not token:
        return None
    payload = verify_token(token, JWT_SECRET)
    if not payload:
        return None
    if payload.get("exp") and datetime.utcnow().timestamp() > payload["exp"]:
        return None
    return payload.get("uid")

@app.post("/auth/register")
def register(data: RegisterInput):
    conn = get_conn()
    cur = conn.cursor()
    uid = str(uuid.uuid4())
    salt = os.urandom(16)
    pwd = pbkdf2_hash(data.password, salt)
    try:
        cur.execute("INSERT INTO users(id,email,password_salt,password_hash,password_plain,created_at) VALUES(?,?,?,?,?,?)", (
            uid, data.email, salt, pwd, data.password, datetime.utcnow().isoformat()
        ))
        conn.commit()
    except sqlite3.IntegrityError:
        conn.close()
        raise HTTPException(status_code=400, detail="email_exists")
    conn.close()
    return {"ok": True}

@app.post("/auth/login")
def login(data: LoginInput, response: Response):
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("SELECT id,password_salt,password_hash,password_plain FROM users WHERE email=?", (data.email,))
    row = cur.fetchone()
    conn.close()
    if not row:
        raise HTTPException(status_code=401, detail="invalid_credentials")
    uid, salt, pwd_hash, pwd_plain = row
    ok = False
    if pwd_plain is not None and pwd_plain == data.password:
        ok = True
    else:
        try:
            if pbkdf2_hash(data.password, salt) == pwd_hash:
                ok = True
        except Exception:
            ok = False
    if not ok:
        raise HTTPException(status_code=401, detail="invalid_credentials")
    exp = int((datetime.utcnow() + timedelta(hours=8)).timestamp())
    token = sign_token({"uid": uid, "exp": exp}, JWT_SECRET)
    response.set_cookie("token", token, httponly=True, samesite="lax", secure=False)
    return {"ok": True}

@app.post("/auth/logout")
def logout(response: Response):
    response.delete_cookie("token")
    return {"ok": True}

@app.get("/auth/me")
def me(request: Request):
    uid = auth_user(request)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("SELECT email FROM users WHERE id=?", (uid,))
    row = cur.fetchone()
    conn.close()
    if not row:
        raise HTTPException(status_code=404, detail="not_found")
    return {"email": row[0]}

def type_to_prop(t: str) -> str:
    if t in ("smol", "af2", "dock", "md", "analysis", "mmpbsa"):
        return t
    if t == "md_prepare":
        return "md"
    return t or "misc"

def process_job_background(uid: str, job_id: str, job_type: str, inputs: dict, outputs_dir: str):
    conn = get_conn()
    cur = conn.cursor()
    try:
        module_map = {
            "smol": "http://mod-smol:8001/predict",
            "dock": "http://mod-dock:8003/run",
            "md": "http://mod-md:8004/run",
            "md_prepare": "http://mod-md:8004/prepare",
            "analysis": "http://mod-analysis:8005/run",
            "af2": "http://mod-af2:8002/predict",
            "mmpbsa": "http://mod-mmpbsa:8008/run",
        }
        url = module_map.get(job_type)
        if not url:
            raise RuntimeError("unknown_type")
        # Use no timeout to allow long running jobs in background
        with httpx.Client(timeout=None) as client:
            payload = {"job_id": job_id}
            payload.update(inputs)
            payload["uid"] = uid
            if job_type == "dock":
                # resolve structure file paths to absolute user paths
                def resolve_path(p: Optional[str]) -> Optional[str]:
                    if not p:
                        return None
                    name = os.path.basename(p)
                    d = user_structure_dir(uid)
                    return os.path.join(d, name)
                if payload.get("receptor") and isinstance(payload["receptor"], dict):
                    rp = payload["receptor"].get("path")
                    payload["receptor"]["path"] = resolve_path(rp)
                if payload.get("ligand") and isinstance(payload["ligand"], dict):
                    lp = payload["ligand"].get("path")
                    payload["ligand"]["path"] = resolve_path(lp)
            elif job_type == "md":
                def resolve_path_md(p: Optional[str]) -> Optional[str]:
                    if not p:
                        return None
                    name = os.path.basename(p)
                    d = user_structure_dir(uid)
                    return os.path.join(d, name)
                if payload.get("protein") and isinstance(payload["protein"], dict):
                    pp = payload["protein"].get("path")
                    payload["protein"]["path"] = resolve_path_md(pp)
            elif job_type == "md_prepare":
                def resolve_path_md_prep(p: Optional[str]) -> Optional[str]:
                    if not p:
                        return None
                    name = os.path.basename(p)
                    d = user_structure_dir(uid)
                    return os.path.join(d, name)
                if payload.get("protein") and isinstance(payload["protein"], dict):
                    pp = payload["protein"].get("path")
                    payload["protein"]["path"] = resolve_path_md_prep(pp)
            try:
                with open(os.path.join(outputs_dir, "log.txt"), "a", encoding="utf-8") as lf:
                    lf.write("dispatching_module_background\n")
            except Exception:
                pass
            resp = client.post(url, json=payload)
        if resp.status_code != 200:
            raise RuntimeError(f"module_error:{resp.status_code}")
        out = resp.json()
        outputs = out.get("outputs", {})
        
        job_prop_local = type_to_prop(job_type)
        job_dir = os.path.join(DATA_DIR, uid, job_prop_local, job_id, "outputs")
        os.makedirs(job_dir, exist_ok=True)
        persisted_outputs = {}
        if job_type == "smol":
            written = {}
            props_items = outputs.get("props_items") or []
            with open(os.path.join(job_dir, "props.json"), "w", encoding="utf-8") as f:
                json.dump({"items": props_items}, f, ensure_ascii=False)
            written["props"] = "props.json"
            molblocks = outputs.get("molblocks") or []
            smiles_list = outputs.get("smiles_list") or []
            with open(os.path.join(job_dir, "out.sdf"), "w", encoding="utf-8") as f:
                for mb in molblocks:
                    f.write(mb); f.write("\n$$$$\n")
            written["sdf"] = "out.sdf"
            sdf_list = []
            for i, mb in enumerate(molblocks):
                p = os.path.join(job_dir, f"mol_{i}.sdf")
                with open(p, "w", encoding="utf-8") as f:
                    f.write(mb)
                entry = {"index": i, "sdf": f"mol_{i}.sdf", "smiles": smiles_list[i] if i < len(smiles_list) else ""}
                png_list = outputs.get("png_list") or []
                if i < len(png_list) and png_list[i]:
                    import base64
                    data_b = base64.b64decode(png_list[i])
                    pngp = os.path.join(job_dir, f"mol_{i}.png")
                    with open(pngp, "wb") as pf:
                        pf.write(data_b)
                    entry["png"] = f"mol_{i}.png"
                sdf_list.append(entry)
            persisted_outputs = {"sdf": written.get("sdf"), "props": written.get("props"), "sdf_list": sdf_list}
        elif job_type == "dock":
            # write docking outputs
            pose_sdf = outputs.get("pose_sdf")
            pose_png = outputs.get("pose_png")
            scores = outputs.get("scores") or {}
            if pose_sdf:
                with open(os.path.join(job_dir, "pose.sdf"), "w", encoding="utf-8") as f:
                    f.write(pose_sdf)
                persisted_outputs["pose_sdf"] = "pose.sdf"
            if pose_png:
                import base64
                data_b = base64.b64decode(pose_png)
                with open(os.path.join(job_dir, "pose.png"), "wb") as pf:
                    pf.write(data_b)
                persisted_outputs["pose_png"] = "pose.png"
            with open(os.path.join(job_dir, "scores.json"), "w", encoding="utf-8") as f:
                json.dump(scores, f, ensure_ascii=False)
            persisted_outputs["scores"] = "scores.json"
            modes_json = outputs.get("modes_json")
            if modes_json and isinstance(modes_json, str):
                persisted_outputs["modes_json"] = modes_json
            pose_pdbqt = outputs.get("pose_pdbqt")
            if pose_pdbqt and isinstance(pose_pdbqt, str):
                persisted_outputs["pose_pdbqt"] = pose_pdbqt
            complex_rel = outputs.get("complex_pdb")
            if complex_rel and isinstance(complex_rel, str):
                # Module already wrote file to outputs; persist relative name
                persisted_outputs["complex_pdb"] = complex_rel
            for k in ("receptor_pdb", "pose_pdb"):
                rel = outputs.get(k)
                if rel and isinstance(rel, str):
                    persisted_outputs[k] = rel
        elif job_type == "md":
            for key, rel in outputs.items():
                if rel and isinstance(rel, str):
                    persisted_outputs[key] = rel
            # fallback: attach known files if module omitted keys
            known = {
                "traj_pdb": "md.pdb",
                "report": "report.csv",
                "em_log": "em.log",
                "npt_log": "npt.log",
                "prod_log": "prod.log",
                "em_plot": "em_potential.png",
                "npt_temp_plot": "npt_temp.png",
                "npt_press_plot": "npt_pressure.png",
                "prod_temp_plot": "prod_temp.png",
                "prod_press_plot": "prod_pressure.png",
            }
            for k, fname in known.items():
                try:
                    if k not in persisted_outputs and os.path.exists(os.path.join(job_dir, fname)):
                        persisted_outputs[k] = fname
                except Exception:
                    pass
        elif job_type == "analysis":
            for key, rel in outputs.items():
                if rel and isinstance(rel, str):
                    persisted_outputs[key] = rel
        elif job_type == "mmpbsa":
            for key, rel in outputs.items():
                if rel and isinstance(rel, str):
                    persisted_outputs[key] = rel
        elif job_type == "md_prepare":
            log_rel = outputs.get("log")
            top_rel = outputs.get("topol")
            gro_rel = outputs.get("gro")
            posre_rel = outputs.get("posre_itp")
            lig_itp_rel = outputs.get("ligand_itp")
            lig_atomtypes_rel = outputs.get("ligand_atomtypes_itp")
            for key, rel in [("log", log_rel), ("topol", top_rel), ("gro", gro_rel), ("posre_itp", posre_rel), ("ligand_itp", lig_itp_rel), ("ligand_atomtypes_itp", lig_atomtypes_rel)]:
                if rel and isinstance(rel, str):
                    persisted_outputs[key] = rel
        cur.execute("UPDATE jobs SET outputs_json=?, status=? WHERE id=?", (
            json.dumps(persisted_outputs), "succeeded", job_id
        ))
        conn.commit()
    except Exception as e:
        try:
            with open(os.path.join(outputs_dir, "log.txt"), "a", encoding="utf-8") as lf:
                if job_type == "md":
                    lf.write("terminated:md\n")
                lf.write(f"backend_error:{str(e)}\n")
        except Exception:
            pass
        cur.execute("UPDATE jobs SET status=?, metrics_json=? WHERE id=?", (
            "failed", json.dumps({"error": str(e)}), job_id
        ))
        conn.commit()
    finally:
        conn.close()

@app.post("/jobs")
def create_job(req: Request, data: JobCreate, background_tasks: BackgroundTasks):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    job_id = str(uuid.uuid4())
    
    job_prop = type_to_prop(data.type)
    job_dir = os.path.join(DATA_DIR, uid, job_prop, job_id)
    inputs_dir = os.path.join(job_dir, "inputs")
    outputs_dir = os.path.join(job_dir, "outputs")
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)
    
    with open(os.path.join(inputs_dir, "inputs.json"), "w", encoding="utf-8") as f:
        json.dump(data.inputs, f, ensure_ascii=False)
        
    try:
        with open(os.path.join(outputs_dir, "log.txt"), "a", encoding="utf-8") as lf:
            lf.write("job_created\n")
    except Exception:
        pass
        
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("INSERT INTO jobs(id,uid,type,status,created_at,inputs_json,outputs_json,metrics_json) VALUES(?,?,?,?,?,?,?,?)", (
        job_id, uid, data.type, "running", datetime.utcnow().isoformat(), json.dumps(data.inputs), json.dumps({}), json.dumps({})
    ))
    conn.commit()
    conn.close()
    
    background_tasks.add_task(process_job_background, uid, job_id, data.type, data.inputs, outputs_dir)
    
    return {"job_id": job_id}

@app.get("/jobs/{job_id}")
def get_job(req: Request, job_id: str):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("SELECT id,uid,type,status,created_at,inputs_json,outputs_json,metrics_json FROM jobs WHERE id=?", (job_id,))
    row = cur.fetchone()
    conn.close()
    if not row:
        raise HTTPException(status_code=404, detail="not_found")
    if row[1] != uid:
        raise HTTPException(status_code=403, detail="forbidden")
    return {
        "id": row[0],
        "type": row[2],
        "status": row[3],
        "created_at": row[4],
        "inputs": json.loads(row[5] or "{}"),
        "outputs": json.loads(row[6] or "{}"),
        "metrics": json.loads(row[7] or "{}"),
    }

@app.get("/jobs/md/dir")
def list_md_dirs(req: Request):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    base_dir = os.path.join(DATA_DIR, uid, "md")
    jobs = []
    try:
        if os.path.isdir(base_dir):
            for name in os.listdir(base_dir):
                d = os.path.join(base_dir, name)
                if os.path.isdir(d):
                    outputs_d = os.path.join(d, "outputs")
                    work_d = os.path.join(outputs_d, "work")
                    has_work = os.path.isdir(work_d)
                    jobs.append({"id": name, "has_work": has_work})
    except Exception:
        pass
    conn = get_conn()
    cur = conn.cursor()
    enriched = []
    for item in jobs:
        try:
            cur.execute("SELECT type,status,created_at FROM jobs WHERE id=? AND uid=?", (item["id"], uid))
            row = cur.fetchone()
            if row and row[0] == "md":
                enriched.append({"id": item["id"], "status": row[1], "created_at": row[2], "has_work": item["has_work"]})
            else:
                enriched.append(item)
        except Exception:
            enriched.append(item)
    conn.close()
    return {"jobs": enriched}

@app.get("/jobs")
def list_jobs(req: Request):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("SELECT id,type,status,created_at FROM jobs WHERE uid=? ORDER BY created_at DESC LIMIT 100", (uid,))
    rows = cur.fetchall()
    conn.close()
    return [{"id": r[0], "type": r[1], "status": r[2], "created_at": r[3]} for r in rows]

@app.get("/files/{job_id}/{artifact:path}")
def download_file(req: Request, job_id: str, artifact: str, download: bool = False):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("SELECT type FROM jobs WHERE id=?", (job_id,))
    row = cur.fetchone()
    conn.close()
    if not row:
        raise HTTPException(status_code=404, detail="not_found")
    job_type = row[0]
    def type_to_prop(t: str) -> str:
        if t in ("smol", "af2", "dock", "md", "analysis"):
            return t
        if t == "md_prepare":
            return "md"
        return t or "misc"
    job_prop = type_to_prop(job_type)
    path = os.path.join(DATA_DIR, uid, job_prop, job_id, "outputs", artifact)
    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail="not_found")
    ext = os.path.splitext(path)[1].lower()
    mt = "application/octet-stream"
    if ext == ".json":
        mt = "application/json"
    elif ext == ".png":
        mt = "image/png"
    elif ext in (".sdf", ".mol"):
        mt = "chemical/x-mdl-sdfile"
    with open(path, "rb") as f:
        data = f.read()
    headers = {}
    if download:
        headers["Content-Disposition"] = f"attachment; filename=\"{os.path.basename(path)}\""
    return Response(content=data, media_type=mt, headers=headers)

@app.post("/files/by-artifact")
def download_file_by_artifact(req: Request, data: FileArtifactInput):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    job_id = data.job_id
    artifact = data.artifact
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("SELECT type FROM jobs WHERE id=?", (job_id,))
    row = cur.fetchone()
    conn.close()
    if not row:
        raise HTTPException(status_code=404, detail="not_found")
    job_type = row[0]
    def type_to_prop(t: str) -> str:
        if t in ("smol", "af2", "dock", "md", "analysis"):
            return t
        if t == "md_prepare":
            return "md"
        return t or "misc"
    job_prop = type_to_prop(job_type)
    base = os.getenv("ENZYME_DATA_DIR", DATA_DIR)
    path = os.path.join(base, uid, job_prop, job_id, "outputs", artifact)
    if not os.path.exists(path):
        alt_base = os.getenv("ALT_DATA_DIR", "/data")
        alt_path = os.path.join(alt_base, uid, job_prop, job_id, "outputs", artifact)
        if os.path.exists(alt_path):
            path = alt_path
        else:
            raise HTTPException(status_code=404, detail="not_found")
    ext = os.path.splitext(path)[1].lower()
    mt = "application/octet-stream"
    if ext == ".json":
        mt = "application/json"
    elif ext == ".png":
        mt = "image/png"
    elif ext in (".sdf", ".mol"):
        mt = "chemical/x-mdl-sdfile"
    elif ext in (".pdb", ".ent"):
        mt = "chemical/x-pdb"
    with open(path, "rb") as f:
        data_b = f.read()
    return Response(content=data_b, media_type=mt)

@app.get("/jobs/{job_id}/files")
def list_job_files(req: Request, job_id: str):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("SELECT type FROM jobs WHERE id=?", (job_id,))
    row = cur.fetchone()
    conn.close()
    if not row:
        raise HTTPException(status_code=404, detail="not_found")
    job_type = row[0]
    def type_to_prop(t: str) -> str:
        if t in ("smol", "af2", "dock", "md", "analysis"):
            return t
        if t == "md_prepare":
            return "md"
        return t or "misc"
    job_prop = type_to_prop(job_type)
    job_dir = os.path.join(DATA_DIR, uid, job_prop, job_id, "outputs")
    if not os.path.isdir(job_dir):
        raise HTTPException(status_code=404, detail="not_found")
    files = []
    for name in os.listdir(job_dir):
        fp = os.path.join(job_dir, name)
        try:
            sz = os.path.getsize(fp)
        except Exception:
            sz = 0
        files.append({"name": name, "size": sz})
    return {"files": files}

@app.get("/user/structure/file/{name}")
def get_user_structure(req: Request, name: str):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    d = user_structure_dir(uid)
    path = os.path.join(d, os.path.basename(name))
    if not os.path.exists(path):
        raise HTTPException(status_code=404, detail="not_found")
    ext = os.path.splitext(path)[1].lower()
    mt = "application/octet-stream"
    if ext in (".pdb", ".ent"):
        mt = "chemical/x-pdb"
    elif ext in (".sdf", ".mol"):
        mt = "chemical/x-mdl-sdfile"
    elif ext in (".cif", ".mmcif"):
        mt = "chemical/x-cif"
    with open(path, "rb") as f:
        data = f.read()
    return Response(content=data, media_type=mt)

@app.get("/public/structure/{name}")
def get_public_structure(name: str):
    target = None
    base = os.getenv("UPLOAD_ASSETS", DATA_DIR)
    for uid in os.listdir(base):
        d = os.path.join(base, uid, "structure")
        p = os.path.join(d, os.path.basename(name))
        if os.path.exists(p):
            target = p
            break
    if not target:
        raise HTTPException(status_code=404, detail="not_found")
    ext = os.path.splitext(target)[1].lower()
    mt = "application/octet-stream"
    if ext in (".pdb", ".ent"):
        mt = "chemical/x-pdb"
    elif ext in (".sdf", ".mol"):
        mt = "chemical/x-mdl-sdfile"
    elif ext in (".cif", ".mmcif"):
        mt = "chemical/x-cif"
    with open(target, "rb") as f:
        data = f.read()
    return Response(content=data, media_type=mt)
def user_structure_dir(uid: str) -> str:
    # Use UPLOAD_ASSETS environment variable if available, otherwise fallback to default
    base_dir = os.getenv("UPLOAD_ASSETS")
    if base_dir:
        d = os.path.join(base_dir, uid, "structure")
    else:
        d = os.path.join(DATA_DIR, uid, "structure")
    os.makedirs(d, exist_ok=True)
    return d

@app.get("/user/structure/list")
def structure_list(req: Request):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    d = user_structure_dir(uid)
    files = []
    for name in os.listdir(d):
        fp = os.path.join(d, name)
        if os.path.isfile(fp):
            try:
                sz = os.path.getsize(fp)
            except Exception:
                sz = 0
            files.append({"name": name, "size": sz})
    return {"files": files}

@app.post("/user/structure/upload")
def structure_upload(req: Request, file: UploadFile = File(...)):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    d = user_structure_dir(uid)
    dst = os.path.join(d, os.path.basename(file.filename))
    data = file.file.read()
    with open(dst, "wb") as f:
        f.write(data)
    return {"ok": True, "name": os.path.basename(file.filename)}
