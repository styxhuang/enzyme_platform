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
VINA_TIMEOUT = int(os.getenv("VINA_TIMEOUT_SECONDS", "600"))

class ReceptorSel(BaseModel):
    path: Optional[str] = None

class LigandSel(BaseModel):
    path: Optional[str] = None
    smiles: Optional[str] = None

class DockInput(BaseModel):
    job_id: str
    uid: str
    receptor: ReceptorSel
    ligand: LigandSel
    center: List[float]
    size: List[float]
    exhaustiveness: int

app = FastAPI()

def which(cmd: str) -> Optional[str]:
    for p in os.environ.get("PATH", "").split(":"):
        fp = os.path.join(p, cmd)
        if os.path.isfile(fp) and os.access(fp, os.X_OK):
            return fp
    return None

def append_log(path: str, text: str):
    try:
        with open(path, "a", encoding="utf-8") as f:
            f.write(text + "\n")
    except Exception:
        pass

def ensure_obabel_ready(log_path: str) -> None:
    ob = obabel_bin()
    if not ob:
        append_log(log_path, "obabel_missing")
        raise HTTPException(status_code=500, detail="openbabel_cli_missing")
    try:
        p = subprocess.run([ob, "-V"], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
        append_log(log_path, "obabel_version " + (p.stdout or ""))
    except Exception as e:
        append_log(log_path, "obabel_unavailable " + str(e))
        raise HTTPException(status_code=500, detail="openbabel_runtime_error")

def obabel_bin() -> Optional[str]:
    prefer = os.environ.get("OBABEL_BIN")
    if prefer and os.path.isfile(prefer) and os.access(prefer, os.X_OK):
        return prefer
    for cand in ["/usr/bin/obabel", "/usr/local/bin/obabel"]:
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return which("obabel")

def sanitize_ligand_pdb_file(path: str) -> None:
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.readlines()
        atoms = []
        conect = []
        for ln in lines:
            t = ln.strip()
            if not t:
                continue
            if t.startswith("ATOM") or t.startswith("HETATM"):
                toks = ln.split()
                if len(toks) >= 8:
                    try:
                        serial = int(toks[1]) if toks[1].isdigit() else len(atoms) + 1
                    except Exception:
                        serial = len(atoms) + 1
                    name = toks[2]
                    try:
                        x = float(toks[5]); y = float(toks[6]); z = float(toks[7])
                    except Exception:
                        continue
                    occ = float(toks[8]) if len(toks) > 8 else 1.00
                    tf = float(toks[9]) if len(toks) > 9 else 0.00
                    elem = toks[10] if len(toks) > 10 else (name.strip()[0] if name.strip() else "C")
                    res_seq = int(toks[4]) if len(toks) > 4 and toks[4].isdigit() else 1
                    line = (
                        f"HETATM{serial:5d} {name:<4} {'':1}{'LIG':>3} {'':1}{'':1}{res_seq:4d}{'':4}   "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}{occ:6.2f}{tf:6.2f}          {elem:>2}"
                    )
                    atoms.append(line + "\n")
            elif t.startswith("CONECT"):
                conect.append(t + "\n")
        out = []
        out.append("HEADER    LIGAND\n")
        out.append("CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n")
        out.extend(atoms)
        out.extend(conect)
        out.append("TER\n")
        out.append("END\n")
        with open(path, "w", encoding="utf-8") as f:
            f.writelines(out)
    except Exception:
        pass

class PDBAtom:
    def __init__(self, record: str, serial: int, name: str, resName: str, chainID: str, resSeq: int, x: float, y: float, z: float, occupancy: float, tempFactor: float, element: str):
        self.record = record
        self.serial = serial
        self.name = name
        self.resName = resName
        self.chainID = chainID
        self.resSeq = resSeq
        self.x = x
        self.y = y
        self.z = z
        self.occupancy = occupancy
        self.tempFactor = tempFactor
        self.element = element

def _parse_atom_line_fixed(ln: str) -> Optional[PDBAtom]:
    try:
        record = ln[0:6].strip()
        serial = int(ln[6:11].strip()) if ln[6:11].strip() else 0
        name = ln[12:16].strip() or "C"
        resName = ln[17:20].strip() or "RES"
        chainID = ln[21:22].strip() or ""
        resSeq = int(ln[22:26].strip()) if ln[22:26].strip() else 1
        x = float(ln[30:38].strip())
        y = float(ln[38:46].strip())
        z = float(ln[46:54].strip())
        occupancy = float(ln[54:60].strip()) if ln[54:60].strip() else 1.00
        tempFactor = float(ln[60:66].strip()) if ln[60:66].strip() else 0.00
        element = (ln[76:78].strip() or name.strip()[0]).upper()
        return PDBAtom(record, serial, name, resName, chainID, resSeq, x, y, z, occupancy, tempFactor, element)
    except Exception:
        return None

def _parse_atom_line_split(ln: str) -> Optional[PDBAtom]:
    try:
        toks = ln.split()
        record = toks[0]
        serial = int(toks[1]) if len(toks) > 1 and toks[1].isdigit() else 0
        name = toks[2] if len(toks) > 2 else "C"
        resName = (toks[3] if len(toks) > 3 else "RES")
        resSeq = int(toks[4]) if len(toks) > 4 and toks[4].isdigit() else 1
        x = float(toks[5]); y = float(toks[6]); z = float(toks[7])
        occupancy = float(toks[8]) if len(toks) > 8 else 1.00
        tempFactor = float(toks[9]) if len(toks) > 9 else 0.00
        element = toks[10] if len(toks) > 10 else name.strip()[0]
        return PDBAtom(record, serial, name, resName, "", resSeq, x, y, z, occupancy, tempFactor, element)
    except Exception:
        return None

def parse_pdb_file(path: str) -> dict:
    model = {"header": None, "cryst1": None, "atoms": [], "conect": {}}
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if ln.startswith("HEADER"):
                    model["header"] = ln.rstrip("\n")
                elif ln.startswith("CRYST1"):
                    model["cryst1"] = ln.rstrip("\n")
                elif ln.startswith("ATOM") or ln.startswith("HETATM"):
                    a = _parse_atom_line_fixed(ln) or _parse_atom_line_split(ln)
                    if a:
                        model["atoms"].append(a)
                elif ln.startswith("CONECT"):
                    parts = ln.split()
                    if len(parts) >= 2 and parts[1].isdigit():
                        src = int(parts[1])
                        nbrs = []
                        for p in parts[2:]:
                            if p.isdigit():
                                nbrs.append(int(p))
                        if src not in model["conect"]:
                            model["conect"][src] = set()
                        for n in nbrs:
                            model["conect"][src].add(n)
    except Exception:
        pass
    return model

def _format_atom_line(atom: PDBAtom) -> str:
    return (
        f"{atom.record:<6}{atom.serial:5d} {atom.name:<4} "
        f"{atom.resName:>3} {atom.chainID:1}{atom.resSeq:4d}    "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}{atom.occupancy:6.2f}{atom.tempFactor:6.2f}          {atom.element:>2}\n"
    )

def write_complex_pdb(out_path: str, receptor_path: str, ligand_path: str) -> None:
    rec = parse_pdb_file(receptor_path)
    lig = parse_pdb_file(ligand_path)
    # Renumber atoms and adjust CONECT
    merged_atoms: List[PDBAtom] = []
    serial_map = {}
    cur = 1
    # receptor atoms
    for a in rec["atoms"]:
        a2 = PDBAtom(a.record, cur, a.name, a.resName, a.chainID, a.resSeq, a.x, a.y, a.z, a.occupancy, a.tempFactor, a.element)
        merged_atoms.append(a2)
        serial_map[("rec", a.serial)] = cur
        cur += 1
    # ligand atoms (force HETATM, resName=LIG)
    for a in lig["atoms"]:
        a2 = PDBAtom("HETATM", cur, a.name, "LIG", a.chainID or "", 1, a.x, a.y, a.z, a.occupancy, a.tempFactor, a.element)
        merged_atoms.append(a2)
        serial_map[("lig", a.serial)] = cur
        cur += 1
    # Build merged CONECT
    merged_conect = {}
    def add_con(src, dst):
        if src not in merged_conect:
            merged_conect[src] = set()
        merged_conect[src].add(dst)
    for src, nbrs in rec["conect"].items():
        if ("rec", src) in serial_map:
            s_new = serial_map[("rec", src)]
            for n in nbrs:
                if ("rec", n) in serial_map:
                    add_con(s_new, serial_map[("rec", n)])
    for src, nbrs in lig["conect"].items():
        if ("lig", src) in serial_map:
            s_new = serial_map[("lig", src)]
            for n in nbrs:
                if ("lig", n) in serial_map:
                    add_con(s_new, serial_map[("lig", n)])
    # Write output in requested order: HEADER, CRYST1, ATOM/HETATM, CONECT, TER, END
    header = rec.get("header") or "HEADER    COMPLEX"
    cryst1 = rec.get("cryst1") or "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1"
    with open(out_path, "w", encoding="utf-8") as f:
        f.write(header + "\n")
        f.write(cryst1 + "\n")
        rec_n = len(rec["atoms"]) if rec and rec.get("atoms") is not None else 0
        for i, a in enumerate(merged_atoms):
            if i == rec_n:
                f.write("TER\n")
            f.write(_format_atom_line(a))
        # write conect lines, sorted
        for src in sorted(merged_conect.keys()):
            nbrs = sorted(merged_conect[src])
            if nbrs:
                f.write("CONECT" + f"{src:5d}" + "".join(f"{n:5d}" for n in nbrs) + "\n")
        f.write("END\n")

def run_vina_local(receptor_pdbqt: str, ligand_pdbqt: str, pose_pdbqt: str, vina_log: str, center: List[float], size: List[float], exhaustiveness: int, n_poses: int = 9, energy_range: int = 3) -> None:
    cmd = [
        "vina",
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--exhaustiveness", str(exhaustiveness),
        "--out", pose_pdbqt,
        "--num_modes", str(n_poses),
        "--energy_range", str(energy_range),
    ]
    with open(vina_log, "a", encoding="utf-8") as lf:
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        if p.stdout is not None:
            for line in p.stdout:
                lf.write(line)
        rc = p.wait()
        if rc != 0:
            raise subprocess.CalledProcessError(rc, cmd)

def parse_vina_log_modes(log_path: str) -> List[dict]:
    modes = []
    if not os.path.exists(log_path):
        return modes
    capture = False
    with open(log_path, "r", encoding="utf-8", errors="ignore") as fh:
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
            if parts and parts[0].isdigit():
                try:
                    mode = int(parts[0])
                    affinity = float(parts[1])
                    rmsd_lb = float(parts[2]) if len(parts) > 2 else None
                    rmsd_ub = float(parts[3]) if len(parts) > 3 else None
                    modes.append({"mode": mode, "affinity": affinity, "rmsd_lb": rmsd_lb, "rmsd_ub": rmsd_ub})
                except Exception:
                    continue
    return modes

@app.post("/run")
def run(data: DockInput):
    if not data.job_id:
        raise HTTPException(status_code=400, detail="missing_job_id")
    job_dir = os.path.join(DATA_DIR, data.uid, "dock", data.job_id)
    inputs_dir = os.path.join(job_dir, "inputs")
    os.makedirs(inputs_dir, exist_ok=True)
    with open(os.path.join(inputs_dir, "inputs.json"), "w", encoding="utf-8") as f:
        json.dump(data.model_dump(), f, ensure_ascii=False)
    outputs_dir = os.path.join(job_dir, "outputs")
    os.makedirs(outputs_dir, exist_ok=True)
    log_path = os.path.join(outputs_dir, "log.txt")
    append_log(log_path, "dock_start")
    ensure_obabel_ready(log_path)
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
        # Prefer MGLTools prepare_receptor4.py
        pysh = which("pythonsh") or "/usr/local/mgltools/bin/pythonsh"
        prep_py = os.environ.get("ADT_PREP_RECEPTOR", "/usr/local/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py")
        if pysh and os.path.exists(prep_py):
            try:
                cmd = [pysh, prep_py, "-r", receptor_src, "-o", receptor_pdbqt, "-A", "hydrogens", "-U", "nphs_lps_waters"]
                p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
                append_log(log_path, "prepare_receptor4 " + " ".join(cmd))
                append_log(log_path, p.stdout or "")
            except Exception as e:
                append_log(log_path, "prepare_receptor4_error " + str(e))
                raise HTTPException(status_code=500, detail=f"prepare_receptor4_failed: {e}")
        else:
            in_fmt = "pdb" if rec_ext in (".pdb", ".ent") else rec_ext[1:] if rec_ext else "pdb"
            cmd = [
                obabel_bin() or "obabel",
                f"-i{in_fmt}", receptor_src,
                "-opdbqt", "-O", receptor_pdbqt,
                "-xh", "-p", "7.4",
                "-xw",
                "--partialcharge", "gasteiger",
            ]
            try:
                p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
                append_log(log_path, "obabel_receptor_pdbqt " + " ".join(cmd))
                append_log(log_path, p.stdout or "")
            except Exception as e:
                append_log(log_path, "obabel_receptor_pdbqt_error " + str(e))
                raise
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
        try:
            cmd = [obabel_bin() or "obabel", "-isdf", tmp_sdf, "-opdbqt", "-O", ligand_pdbqt]
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
            append_log(log_path, "obabel_ligand_pdbqt " + " ".join(cmd))
            append_log(log_path, p.stdout or "")
        except Exception as e:
            append_log(log_path, "obabel_ligand_pdbqt_error " + str(e))
            raise
    cx, cy, cz = float(data.center[0]), float(data.center[1]), float(data.center[2])
    sx = float(data.size[0]); sy = float(data.size[1]); sz = float(data.size[2])
    ex = int(data.exhaustiveness)
    pose_pdbqt = os.path.join(outputs_dir, "pose.pdbqt")
    vina_log = os.path.join(outputs_dir, "vina.log")
    energy = None
    append_log(log_path, "vina_local_start")
    try:
        run_vina_local(
            receptor_pdbqt=receptor_pdbqt,
            ligand_pdbqt=ligand_pdbqt,
            pose_pdbqt=pose_pdbqt,
            vina_log=vina_log,
            center=[cx, cy, cz],
            size=[sx, sy, sz],
            exhaustiveness=ex,
        )
        append_log(log_path, "vina_local_done")
    except Exception as e:
        append_log(log_path, f"vina_request_failed {e}")
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
        # Fallback: build complex using receptor PDB and ligand translated to center
        try:
            receptor_pdb = os.path.join(outputs_dir, "receptor.pdb")
            if rec_ext == ".pdb":
                try:
                    with open(receptor_src, "r", encoding="utf-8", errors="ignore") as rf, open(receptor_pdb, "w", encoding="utf-8") as wf:
                        wf.write(rf.read())
                    append_log(log_path, "receptor_pdb_copied")
                except Exception as ce:
                    in_fmt = "pdb" if rec_ext in (".pdb", ".ent") else rec_ext[1:] if rec_ext else "pdb"
                    cmd = [obabel_bin() or "obabel", f"-i{in_fmt}", receptor_src, "-opdb", "-O", receptor_pdb]
                    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
                    append_log(log_path, "obabel_receptor_pdb " + " ".join(cmd))
                    append_log(log_path, p.stdout or "")
            else:
                in_fmt = "pdb" if rec_ext in (".pdb", ".ent") else rec_ext[1:] if rec_ext else "pdb"
                cmd = [obabel_bin() or "obabel", f"-i{in_fmt}", receptor_src, "-opdb", "-O", receptor_pdb]
                p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
                append_log(log_path, "obabel_receptor_pdb " + " ".join(cmd))
                append_log(log_path, p.stdout or "")

            # Translate ligand coordinates to box center (visual placement)
            conf = lig.GetConformer()
            import numpy as np
            pts = np.array([list(conf.GetAtomPosition(i)) for i in range(lig.GetNumAtoms())])
            ctr = pts.mean(axis=0)
            delta = np.array([cx, cy, cz]) - ctr
            for i in range(lig.GetNumAtoms()):
                p = conf.GetAtomPosition(i)
                conf.SetAtomPosition(i, Chem.rdGeometry.Point3D(p.x + float(delta[0]), p.y + float(delta[1]), p.z + float(delta[2])))
            pose_pdb = os.path.join(outputs_dir, "pose.pdb")
            with open(pose_pdb, "w", encoding="utf-8") as pf:
                pf.write(Chem.MolToPDBBlock(lig))
            complex_pdb = os.path.join(outputs_dir, "complex.pdb")
            with open(complex_pdb, "w", encoding="utf-8") as outp:
                with open(receptor_pdb, "r", encoding="utf-8", errors="ignore") as rp:
                    rtxt = rp.read(); outp.write(rtxt);
                    if not rtxt.endswith("\n"): outp.write("\n")
                with open(pose_pdb, "r", encoding="utf-8", errors="ignore") as lp:
                    outp.write(lp.read())
            append_log(log_path, "fallback_complex_written")
            outputs_obj = {"pose_sdf": mb, "pose_png": png_b64, "scores": scores, "receptor_pdb": "receptor.pdb", "pose_pdb": "pose.pdb", "complex_pdb": "complex.pdb"}
        except Exception as fe:
            append_log(log_path, f"fallback_complex_error {fe}")
            outputs_obj = {"pose_sdf": mb, "pose_png": png_b64, "scores": scores}
        return {"job_id": data.job_id, "outputs": outputs_obj}
    # Split all modes and build per-mode complexes
    modes_info = parse_vina_log_modes(vina_log)
    poses_dir = os.path.join(outputs_dir, "poses")
    os.makedirs(poses_dir, exist_ok=True)
    # Split ligand models to PDB using obabel (-m)
    try:
        cmd = [obabel_bin() or "obabel", "-ipdbqt", pose_pdbqt, "-opdb", "-O", os.path.join(poses_dir, "ligand_mode.pdb"), "-m"]
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
        append_log(log_path, "obabel_split_modes " + " ".join(cmd))
        append_log(log_path, p.stdout or "")
    except Exception as e:
        append_log(log_path, "obabel_split_modes_error " + str(e))
    # Build complexes per mode
    modes_manifest = []
    try:
        receptor_pdb_for_merge = os.path.join(outputs_dir, "receptor.pdb")
        # ensure receptor PDB exists
        if not os.path.exists(receptor_pdb_for_merge):
            in_fmt = "pdb" if rec_ext in (".pdb", ".ent") else rec_ext[1:] if rec_ext else "pdb"
            subprocess.run([obabel_bin() or "obabel", f"-i{in_fmt}", receptor_src, "-opdb", "-O", receptor_pdb_for_merge], check=True)
        # add hydrogens to receptor once
        receptor_pdb_h = os.path.join(poses_dir, "receptor_h.pdb")
        try:
            cmd_hr = [obabel_bin() or "obabel", "-ipdb", receptor_pdb_for_merge, "-opdb", "-O", receptor_pdb_h, "--add-hydrogens", "-p", "7.4"]
            p = subprocess.run(cmd_hr, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
            append_log(log_path, "obabel_add_h_receptor " + " ".join(cmd_hr))
            append_log(log_path, p.stdout or "")
        except Exception as hre:
            append_log(log_path, "obabel_add_h_receptor_error " + str(hre))
            receptor_pdb_h = receptor_pdb_for_merge
        lig_files = sorted([fn for fn in os.listdir(poses_dir) if fn.startswith("ligand_mode") and fn.endswith(".pdb")])
        for idx, fn in enumerate(lig_files, start=1):
            lig_pdb_path = os.path.join(poses_dir, fn)
            complex_path = os.path.join(poses_dir, f"complex_mode{idx}.pdb")
            try:
                # add hydrogens to ligand before sanitizing and merging
                lig_h_path = os.path.join(poses_dir, f"{os.path.splitext(fn)[0]}_h.pdb")
                try:
                    cmd_hl = [obabel_bin() or "obabel", "-ipdb", lig_pdb_path, "-opdb", "-O", lig_h_path, "--add-hydrogens", "-p", "7.4"]
                    p = subprocess.run(cmd_hl, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True)
                    append_log(log_path, "obabel_add_h_ligand " + " ".join(cmd_hl))
                    append_log(log_path, p.stdout or "")
                except Exception as hle:
                    append_log(log_path, "obabel_add_h_ligand_error " + str(hle))
                    lig_h_path = lig_pdb_path
                sanitize_ligand_pdb_file(lig_h_path)
                write_complex_pdb(complex_path, receptor_pdb_h, lig_h_path)
                modes_manifest.append({"mode": idx, "ligand_pdb": f"poses/{fn}", "complex_pdb": f"poses/complex_mode{idx}.pdb"})
            except Exception as me:
                append_log(log_path, f"complex_mode{idx}_error {me}")
    except Exception as ee:
        append_log(log_path, f"modes_manifest_error {ee}")
    # Pose SDF (best mode) for backward compatibility
    # pose_sdf = os.path.join(outputs_dir, "pose.sdf")
    # mb = ""
    # try:
    #     subprocess.run([obabel_bin() or "obabel", "-ipdbqt", pose_pdbqt, "-osdf", "-O", pose_sdf], check=True)
    #     with open(pose_sdf, "r", encoding="utf-8") as f:
    #         mb = f.read()
    # except Exception:
    #     pass
    
    
    if modes_info:
        items = [{"mode": m["mode"], "affinity": m["affinity"], "rmsd_lb": m["rmsd_lb"], "rmsd_ub": m["rmsd_ub"]} for m in modes_info]
        scores = {"items": items, "center": data.center, "size": data.size, "exhaustiveness": data.exhaustiveness}
    else:
        scores = {"items": [{"mode": 1, "affinity": energy if energy is not None else None}], "center": data.center, "size": data.size, "exhaustiveness": data.exhaustiveness}
    # Write modes manifest
    try:
        with open(os.path.join(outputs_dir, "modes.json"), "w", encoding="utf-8") as jf:
            json.dump({"modes": modes_manifest}, jf, ensure_ascii=False)
    except Exception:
        pass
    outputs_obj = {"scores": scores, "pose_pdbqt": "pose.pdbqt", "receptor_pdb": "receptor.pdb", "modes_json": "modes.json"}
    return {"job_id": data.job_id, "outputs": outputs_obj}
def which(cmd: str) -> Optional[str]:
    for p in os.environ.get("PATH", "").split(":"):
        fp = os.path.join(p, cmd)
        if os.path.isfile(fp) and os.access(fp, os.X_OK):
            return fp
    return None
