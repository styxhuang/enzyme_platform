import os
import json
from typing import List, Optional
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdDepictor
import base64
import uuid
import traceback

DATA_DIR = os.getenv("ENZYME_DATA_DIR", "/data")

class SmolInput(BaseModel):
    smiles: List[str]
    num_confs: Optional[int] = 3
    minimize: Optional[str] = "MMFF94"
    job_id: Optional[str] = None

app = FastAPI()

def ensure_job_dirs(job_id: str):
    job_dir = os.path.join(DATA_DIR, job_id)
    inputs_dir = os.path.join(job_dir, "inputs")
    outputs_dir = os.path.join(job_dir, "outputs")
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)
    return inputs_dir, outputs_dir

@app.post("/predict")
def predict(data: SmolInput):
    if not data.job_id:
        raise HTTPException(status_code=400, detail="missing_job_id")
    job_id = data.job_id
    inputs_dir, outputs_dir = ensure_job_dirs(job_id)
    with open(os.path.join(inputs_dir, "inputs.json"), "w", encoding="utf-8") as f:
        json.dump(data.model_dump(), f, ensure_ascii=False)
    log_path = os.path.join(outputs_dir, "log.txt")
    def log(msg: str):
        try:
            with open(log_path, "a", encoding="utf-8") as lf:
                lf.write(msg + "\n")
        except Exception:
            pass
    props_rows = []
    molblocks = []
    png_list = []
    smiles_list = []
    for smi in data.smiles:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                log(f"invalid_smiles: {smi}")
                continue
            mol = Chem.AddHs(mol)
            ps = AllChem.ETKDGv3()
            ps.numThreads = 0
            cid = AllChem.EmbedMultipleConfs(mol, numConfs=data.num_confs, params=ps)
            if data.minimize.upper() == "MMFF94":
                for c in cid:
                    AllChem.MMFFOptimizeMolecule(mol, confId=int(c))
            else:
                for c in cid:
                    AllChem.UFFOptimizeMolecule(mol, confId=int(c))
            mb = Chem.MolToMolBlock(mol)
            logp = Descriptors.MolLogP(mol)
            tpsa = Descriptors.TPSA(mol)
            mw = Descriptors.MolWt(mol)
            hbd = Descriptors.NumHDonors(mol)
            hba = Descriptors.NumHAcceptors(mol)
            rot = Descriptors.NumRotatableBonds(mol)
            props_rows.append({"smiles": smi, "logP": logp, "TPSA": tpsa, "MW": mw, "HBD": hbd, "HBA": hba, "RotB": rot})
            # 2D depiction on a copy
            try:
                mol2d = Chem.Mol(mol)
                rdDepictor.Compute2DCoords(mol2d)
                img = Draw.MolToImage(mol2d, size=(300, 200))
                from io import BytesIO
                bio = BytesIO()
                img.save(bio, format='PNG')
                png_b64 = base64.b64encode(bio.getvalue()).decode()
            except Exception:
                log(f"png_render_fail: {smi}")
                png_b64 = ""
            molblocks.append(mb)
            png_list.append(png_b64)
            smiles_list.append(smi)
        except Exception as e:
            log(f"proc_error: {smi} -> {e}\n{traceback.format_exc()}")
    return {"job_id": job_id, "outputs": {"props_items": props_rows, "molblocks": molblocks, "png_list": png_list, "smiles_list": smiles_list}}