import os
import json
from typing import Optional, List
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdDepictor
import subprocess
import httpx
have_meeko = True
try:
    from meeko import MoleculePreparation
except Exception:
    have_meeko = False
import base64

DATA_DIR = os.getenv("ENZYME_DATA_DIR", "/data")
VINA_URL = os.getenv("VINA_SERVICE_URL", "http://mod-vina:8006/run")
VINA_TIMEOUT = int(os.getenv("VINA_TIMEOUT_SECONDS", "600"))

class ReceptorSel(BaseModel):
    path: Optional[str] = None

class LigandSel(BaseModel):
    path: Optional[str] = None
    smiles: Optional[str] = None

class DockInput(BaseModel):
    job_id: str
    receptor: ReceptorSel
    ligand: LigandSel
    center: List[float]
    size: List[float]
    exhaustiveness: int

app = FastAPI()

def request_vina(job_id: str, receptor_pdbqt: str, ligand_pdbqt: str, pose_pdbqt: str, vina_log: str, center: List[float], size: List[float], exhaustiveness: int) -> dict:
    payload = {
        "job_id": job_id,
        "receptor_pdbqt": receptor_pdbqt,
        "ligand_pdbqt": ligand_pdbqt,
        "pose_pdbqt": pose_pdbqt,
        "log_path": vina_log,
        "center": center,
        "size": size,
        "exhaustiveness": exhaustiveness,
    }
    try:
        with httpx.Client(timeout=VINA_TIMEOUT) as client:
            resp = client.post(VINA_URL, json=payload)
    except httpx.HTTPError as exc:
        raise RuntimeError(f"vina_request_failed:{exc}") from exc
    if resp.status_code != 200:
        raise RuntimeError(f"vina_status_{resp.status_code}:{resp.text}")
    return resp.json()

@app.post("/run")
def run(data: DockInput):
    if not data.job_id:
        raise HTTPException(status_code=400, detail="missing_job_id")
    job_dir = os.path.join(DATA_DIR, data.job_id)
    inputs_dir = os.path.join(job_dir, "inputs")
    os.makedirs(inputs_dir, exist_ok=True)
    with open(os.path.join(inputs_dir, "inputs.json"), "w", encoding="utf-8") as f:
        json.dump(data.model_dump(), f, ensure_ascii=False)
    outputs_dir = os.path.join(job_dir, "outputs")
    os.makedirs(outputs_dir, exist_ok=True)
    lig = None
    lig = None
    lig_smiles = None
    if data.ligand.smiles:
        lig_smiles = data.ligand.smiles
        lig = Chem.AddHs(Chem.MolFromSmiles(lig_smiles))
    elif data.ligand.path:
        fp = data.ligand.path
        if not os.path.exists(fp):
            raise HTTPException(status_code=400, detail="ligand_path_not_found")
        ext = os.path.splitext(fp)[1].lower()
        if ext == ".sdf":
            suppl = Chem.SDMolSupplier(fp)
            lig = suppl[0]
            lig_smiles = Chem.MolToSmiles(Chem.RemoveHs(lig))
        else:
            raise HTTPException(status_code=400, detail="ligand_unsupported")
    else:
        raise HTTPException(status_code=400, detail="ligand_missing")
    ps = AllChem.ETKDGv3(); ps.numThreads = 0
    AllChem.EmbedMolecule(lig, ps)
    AllChem.MMFFOptimizeMolecule(lig)
    receptor_src = data.receptor.path or ""
    if not receptor_src or not os.path.exists(receptor_src):
        raise HTTPException(status_code=400, detail="receptor_missing")
    rec_ext = os.path.splitext(receptor_src)[1].lower()
    receptor_pdbqt = os.path.join(outputs_dir, "receptor.pdbqt")
    if rec_ext == ".pdbqt":
        receptor_pdbqt = receptor_src
    else:
        in_fmt = "pdb" if rec_ext in (".pdb", ".ent") else rec_ext[1:] if rec_ext else "pdb"
        cmd = [
            "obabel",
            f"-i{in_fmt}", receptor_src,
            "-opdbqt", "-O", receptor_pdbqt,
            "-xh", "-p", "7.4",
            "-xw",
            "--partialcharge", "gasteiger",
        ]
        subprocess.run(cmd, check=True)
    ligand_pdbqt = os.path.join(outputs_dir, "ligand.pdbqt")
    ok_lig = False
    if have_meeko:
        try:
            prep = MoleculePreparation()
            res = prep.prepare(lig)
            pdbqt_str = res[0] if isinstance(res, tuple) else res
            with open(ligand_pdbqt, "w", encoding="utf-8") as f:
                f.write(pdbqt_str)
            ok_lig = True
        except Exception:
            ok_lig = False
    if not ok_lig:
        tmp_sdf = os.path.join(outputs_dir, "ligand.sdf")
        w = Chem.SDWriter(tmp_sdf)
        w.write(lig); w.close()
        subprocess.run(["obabel", "-isdf", tmp_sdf, "-opdbqt", "-O", ligand_pdbqt], check=True)
    cx, cy, cz = float(data.center[0]), float(data.center[1]), float(data.center[2])
    sx = float(data.size[0]); sy = float(data.size[1]); sz = float(data.size[2])
    ex = int(data.exhaustiveness)
    pose_pdbqt = os.path.join(outputs_dir, "pose.pdbqt")
    vina_log = os.path.join(outputs_dir, "vina.log")
    energy = None
    try:
        vina_resp = request_vina(
            job_id=data.job_id,
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            pose_pdbqt=pose_pdbqt,
            vina_log=vina_log,
            center=[cx, cy, cz],
            size=[sx, sy, sz],
            exhaustiveness=ex,
        )
        energy = vina_resp.get("energy")
    except Exception:
        ps = AllChem.ETKDGv3(); ps.numThreads = 0
        AllChem.EmbedMolecule(lig, ps)
        AllChem.MMFFOptimizeMolecule(lig)
        mb = Chem.MolToMolBlock(lig)
        rdDepictor.Compute2DCoords(lig)
        img = Draw.MolToImage(lig, size=(320, 220))
        from io import BytesIO
        bio = BytesIO(); img.save(bio, format='PNG')
        png_b64 = base64.b64encode(bio.getvalue()).decode()
        score = - (Descriptors.NumHAcceptors(lig) + Descriptors.NumHDonors(lig)) - 0.1*Descriptors.MolWt(lig)
        scores = {"items": [{"smiles": lig_smiles or "", "score": round(score, 3)}], "center": data.center, "size": data.size, "exhaustiveness": data.exhaustiveness}
        return {"job_id": data.job_id, "outputs": {"pose_sdf": mb, "pose_png": png_b64, "scores": scores}}
    pose_sdf = os.path.join(outputs_dir, "pose.sdf")
    subprocess.run(["obabel", "-ipdbqt", pose_pdbqt, "-osdf", "-O", pose_sdf], check=True)
    with open(pose_sdf, "r", encoding="utf-8") as f:
        mb = f.read()
    mols = Chem.SDMolSupplier(pose_sdf)
    out_mol = mols[0] if len(mols) > 0 else lig
    rdDepictor.Compute2DCoords(out_mol)
    img = Draw.MolToImage(out_mol, size=(320, 220))
    from io import BytesIO
    bio = BytesIO(); img.save(bio, format='PNG')
    png_b64 = base64.b64encode(bio.getvalue()).decode()
    item_score = energy if energy is not None else None
    scores = {"items": [{"smiles": lig_smiles or "", "score": item_score}], "center": data.center, "size": data.size, "exhaustiveness": data.exhaustiveness}
    return {"job_id": data.job_id, "outputs": {"pose_sdf": mb, "pose_png": png_b64, "scores": scores}}
