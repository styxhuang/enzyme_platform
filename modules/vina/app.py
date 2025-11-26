import os
import subprocess
from typing import List, Optional
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field

DATA_DIR = os.getenv("ENZYME_DATA_DIR", "/data")

class VinaRunInput(BaseModel):
    job_id: str
    receptor_pdbqt: str
    ligand_pdbqt: str
    pose_pdbqt: str
    log_path: str
    center: List[float]
    size: List[float]
    exhaustiveness: int = Field(default=8, ge=1, le=64)
    n_poses: int = Field(default=9, ge=1, le=20)
    energy_range: int = Field(default=3, ge=1, le=12)

app = FastAPI()

def ensure_data_path(path: str) -> str:
    if not path:
        raise HTTPException(status_code=400, detail="path_missing")
    abs_path = os.path.abspath(path)
    data_root = os.path.abspath(DATA_DIR)
    if not abs_path.startswith(data_root):
        raise HTTPException(status_code=400, detail="path_outside_data_dir")
    os.makedirs(os.path.dirname(abs_path), exist_ok=True)
    return abs_path

def parse_vina_energy(log_path: str) -> Optional[float]:
    if not os.path.exists(log_path):
        return None
    capture = False
    with open(log_path, "r", encoding="utf-8") as fh:
        for line in fh:
            txt = line.strip()
            if not txt:
                continue
            if txt.startswith("-----+"):
                capture = True
                continue
            if not capture:
                continue
            parts = txt.split()
            if parts and parts[0].isdigit() and len(parts) >= 2:
                try:
                    return float(parts[1])
                except ValueError:
                    continue
    return None

def ensure_vector(name: str, values: List[float]) -> List[float]:
    if len(values) != 3:
        raise HTTPException(status_code=400, detail=f"{name}_invalid")
    return [float(values[0]), float(values[1]), float(values[2])]

@app.post("/run")
def run(payload: VinaRunInput):
    receptor = ensure_data_path(payload.receptor_pdbqt)
    ligand = ensure_data_path(payload.ligand_pdbqt)
    pose_out = ensure_data_path(payload.pose_pdbqt)
    log_path = ensure_data_path(payload.log_path)
    center = ensure_vector("center", payload.center)
    size = ensure_vector("size", payload.size)
    cmd = [
        "vina",
        "--receptor", receptor,
        "--ligand", ligand,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--exhaustiveness", str(payload.exhaustiveness),
        "--out", pose_out,
        "--log", log_path,
        "--num_modes", str(payload.n_poses),
        "--energy_range", str(payload.energy_range),
    ]
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as exc:
        raise HTTPException(status_code=500, detail="vina_failed") from exc
    energy = parse_vina_energy(log_path)
    return {"job_id": payload.job_id, "pose_pdbqt": pose_out, "log_path": log_path, "energy": energy}

