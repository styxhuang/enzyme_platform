import os
import json
import uuid
from fastapi import FastAPI
from pydantic import BaseModel

DATA_DIR = os.getenv("ENZYME_DATA_DIR", "/data")

class DockInput(BaseModel):
    receptor_pdb: str
    ligand_smiles: str
    center: list
    size: list
    exhaustiveness: int

app = FastAPI()

def ensure_job_dirs(job_id: str):
    job_dir = os.path.join(DATA_DIR, job_id)
    inputs_dir = os.path.join(job_dir, "inputs")
    outputs_dir = os.path.join(job_dir, "outputs")
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)
    return inputs_dir, outputs_dir

@app.post("/run")
def run(data: DockInput):
    job_id = str(uuid.uuid4())
    inputs_dir, outputs_dir = ensure_job_dirs(job_id)
    with open(os.path.join(inputs_dir, "inputs.json"), "w", encoding="utf-8") as f:
        json.dump(data.model_dump(), f, ensure_ascii=False)
    with open(os.path.join(outputs_dir, "README.txt"), "w", encoding="utf-8") as f:
        f.write("dock placeholder")
    return {"job_id": job_id, "outputs": {"pdb": "best.pdb", "score": "scores.json"}}