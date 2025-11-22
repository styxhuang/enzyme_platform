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
from fastapi import FastAPI, Depends, HTTPException, status, Request, Response
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
import httpx

DATA_DIR = os.getenv("ENZYME_DATA_DIR", os.path.join(os.getcwd(), "data"))
DB_PATH = os.path.join(DATA_DIR, "meta.sqlite")
JWT_SECRET = os.getenv("ENZYME_JWT_SECRET", "dev_secret")

app = FastAPI()
app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://localhost:3000",
        "http://127.0.0.1:3000",
        "http://localhost:5500",
        "http://127.0.0.1:5500",
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

def ensure_dirs():
    os.makedirs(DATA_DIR, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    cur = conn.cursor()
    cur.execute("CREATE TABLE IF NOT EXISTS users(id TEXT PRIMARY KEY, email TEXT UNIQUE, password_salt BLOB, password_hash BLOB, created_at TEXT)")
    cur.execute("CREATE TABLE IF NOT EXISTS jobs(id TEXT PRIMARY KEY, type TEXT, status TEXT, created_at TEXT, inputs_json TEXT, outputs_json TEXT, metrics_json TEXT)")
    conn.commit()
    try:
        cur.execute("ALTER TABLE users ADD COLUMN password_plain TEXT")
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
    response.set_cookie("token", token, httponly=True, samesite="lax")
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

@app.post("/jobs")
def create_job(req: Request, data: JobCreate):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    job_id = str(uuid.uuid4())
    job_dir = os.path.join(DATA_DIR, job_id)
    inputs_dir = os.path.join(job_dir, "inputs")
    outputs_dir = os.path.join(job_dir, "outputs")
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)
    with open(os.path.join(inputs_dir, "inputs.json"), "w", encoding="utf-8") as f:
        json.dump(data.inputs, f, ensure_ascii=False)
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("INSERT INTO jobs(id,type,status,created_at,inputs_json,outputs_json,metrics_json) VALUES(?,?,?,?,?,?,?)", (
        job_id, data.type, "running", datetime.utcnow().isoformat(), json.dumps(data.inputs), json.dumps({}), json.dumps({})
    ))
    conn.commit()
    try:
        module_map = {
            "smol": "http://mod-smol:8001/predict",
            "dock": "http://mod-dock:8003/run",
            "md": "http://mod-md:8004/run",
            "analysis": "http://mod-analysis:8005/run",
            "af2": "http://mod-af2:8002/predict",
        }
        url = module_map.get(data.type)
        if not url:
            raise RuntimeError("unknown_type")
        with httpx.Client(timeout=300) as client:
            payload = {"job_id": job_id}
            payload.update(data.inputs)
            resp = client.post(url, json=payload)
        if resp.status_code != 200:
            raise RuntimeError(f"module_error:{resp.status_code}")
        out = resp.json()
        outputs = out.get("outputs", {})
        job_dir = os.path.join(DATA_DIR, job_id, "outputs")
        # write outputs from module into job_dir (single source of truth)
        os.makedirs(job_dir, exist_ok=True)
        written = {}
        # props
        props_items = outputs.get("props_items") or []
        props_path = os.path.join(job_dir, "props.json")
        with open(props_path, "w", encoding="utf-8") as f:
            json.dump({"items": props_items}, f, ensure_ascii=False)
        written["props"] = "props.json"
        # molblocks -> out.sdf and per-molecule sdf
        molblocks = outputs.get("molblocks") or []
        smiles_list = outputs.get("smiles_list") or []
        out_sdf_path = os.path.join(job_dir, "out.sdf")
        with open(out_sdf_path, "w", encoding="utf-8") as f:
            for mb in molblocks:
                f.write(mb)
                f.write("\n$$$$\n")
        written["sdf"] = "out.sdf"
        sdf_list = []
        for i, mb in enumerate(molblocks):
            p = os.path.join(job_dir, f"mol_{i}.sdf")
            with open(p, "w", encoding="utf-8") as f:
                f.write(mb)
            sdf_list.append({"index": i, "sdf": f"mol_{i}.sdf", "smiles": smiles_list[i] if i < len(smiles_list) else ""})
        # png_list -> per-molecule png
        png_list = outputs.get("png_list") or []
        for i, b64 in enumerate(png_list):
            if not b64:
                continue
            import base64
            data = base64.b64decode(b64)
            with open(os.path.join(job_dir, f"mol_{i}.png"), "wb") as f:
                f.write(data)
            # attach png path
            if i < len(sdf_list):
                sdf_list[i]["png"] = f"mol_{i}.png"
        # persist outputs json returned to client
        persisted_outputs = {"sdf": written.get("sdf"), "props": written.get("props"), "sdf_list": sdf_list}
        cur.execute("UPDATE jobs SET outputs_json=?, status=? WHERE id=?", (
            json.dumps(persisted_outputs), "succeeded", job_id
        ))
        conn.commit()
    except Exception as e:
        cur.execute("UPDATE jobs SET status=?, metrics_json=? WHERE id=?", (
            "failed", json.dumps({"error": str(e)}), job_id
        ))
        conn.commit()
    finally:
        conn.close()
    return {"job_id": job_id}

@app.get("/jobs/{job_id}")
def get_job(req: Request, job_id: str):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("SELECT id,type,status,created_at,inputs_json,outputs_json,metrics_json FROM jobs WHERE id=?", (job_id,))
    row = cur.fetchone()
    conn.close()
    if not row:
        raise HTTPException(status_code=404, detail="not_found")
    return {
        "id": row[0],
        "type": row[1],
        "status": row[2],
        "created_at": row[3],
        "inputs": json.loads(row[4] or "{}"),
        "outputs": json.loads(row[5] or "{}"),
        "metrics": json.loads(row[6] or "{}"),
    }

@app.get("/jobs")
def list_jobs(req: Request):
    uid = auth_user(req)
    if not uid:
        raise HTTPException(status_code=401, detail="unauthorized")
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("SELECT id,type,status,created_at FROM jobs ORDER BY created_at DESC LIMIT 100")
    rows = cur.fetchall()
    conn.close()
    return [{"id": r[0], "type": r[1], "status": r[2], "created_at": r[3]} for r in rows]

@app.get("/files/{job_id}/{artifact}")
def download_file(req: Request, job_id: str, artifact: str):
    path = os.path.join(DATA_DIR, job_id, "outputs", artifact)
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
    return Response(content=data, media_type=mt)

@app.get("/jobs/{job_id}/files")
def list_job_files(job_id: str):
    job_dir = os.path.join(DATA_DIR, job_id, "outputs")
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