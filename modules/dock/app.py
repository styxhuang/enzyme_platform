import os
import json
from typing import Optional, List
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Draw, rdDepictor
import base64

DATA_DIR = os.getenv("ENZYME_DATA_DIR", "/data")

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

@app.post("/run")
def run(data: DockInput):
    if not data.job_id:
        raise HTTPException(status_code=400, detail="missing_job_id")
    job_dir = os.path.join(DATA_DIR, data.job_id)
    inputs_dir = os.path.join(job_dir, "inputs")
    os.makedirs(inputs_dir, exist_ok=True)
    with open(os.path.join(inputs_dir, "inputs.json"), "w", encoding="utf-8") as f:
        json.dump(data.model_dump(), f, ensure_ascii=False)
    # prepare ligand
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
    # embed 3D conformer
    ps = AllChem.ETKDGv3(); ps.numThreads = 0
    AllChem.EmbedMolecule(lig, ps)
    AllChem.MMFFOptimizeMolecule(lig)
    # simple score (placeholder if vina not available)
    score = - (Descriptors.NumHAcceptors(lig) + Descriptors.NumHDonors(lig)) - 0.1*Descriptors.MolWt(lig)
    # outputs
    mb = Chem.MolToMolBlock(lig)
    # 2D depiction
    rdDepictor.Compute2DCoords(lig)
    img = Draw.MolToImage(lig, size=(320, 220))
    from io import BytesIO
    bio = BytesIO(); img.save(bio, format='PNG')
    png_b64 = base64.b64encode(bio.getvalue()).decode()
    scores = {"items": [{"smiles": lig_smiles or "", "score": round(score, 3)}], "center": data.center, "size": data.size, "exhaustiveness": data.exhaustiveness}
    return {"job_id": data.job_id, "outputs": {"pose_sdf": mb, "pose_png": png_b64, "scores": scores}}