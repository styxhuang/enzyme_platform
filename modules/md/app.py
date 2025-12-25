import os
import json
import subprocess
from typing import Optional
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field
import re
import math
import shutil

try:
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except Exception:
    np = None

DATA_DIR = os.getenv("ENZYME_DATA_DIR", "/data")

class MoleculeRef(BaseModel):
    path: str

class MdSteps(BaseModel):
    em_steps: int = Field(default=5000, ge=1, le=1000000)
    em_tol: float = Field(default=1000.0, ge=0.0, le=100000.0)
    npt_steps: int = Field(default=1000, ge=1, le=5000000)
    npt_temp: float = Field(default=300.0, ge=0.0, le=1000.0)
    prod_steps: int = Field(default=1000, ge=1, le=10000000)
    prod_temp: float = Field(default=300.0, ge=0.0, le=1000.0)
    prod_press: float = Field(default=1.0, ge=0.0, le=1000.0)
    box_type: str = Field(default="cubic")
    box_dist: float = Field(default=1.0, ge=0.1, le=5.0)
    ion_neutral: bool = Field(default=True)
    ion_pname: str = Field(default="NA")
    ion_nname: str = Field(default="CL")
    max_frames: int = Field(default=100, ge=1, le=200)

class MdRunInput(BaseModel):
    job_id: str
    uid: str
    protein: MoleculeRef
    forcefield: str = Field(default="amber99sb")
    water: str = Field(default="tip3p")
    steps: MdSteps = Field(default_factory=MdSteps)

app = FastAPI()

def ensure_job_dirs(job_id: str, uid: str):
    job_dir = os.path.join(DATA_DIR, uid, "md", job_id)
    inputs_dir = os.path.join(job_dir, "inputs")
    outputs_dir = os.path.join(job_dir, "outputs")
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)
    return job_dir, inputs_dir, outputs_dir

def ensure_data_path(path: str) -> str:
    if not path:
        raise HTTPException(status_code=400, detail="path_missing")
    abs_path = os.path.abspath(path)
    data_root = os.path.abspath(DATA_DIR)
    if not abs_path.startswith(data_root):
        raise HTTPException(status_code=400, detail="path_outside_data_dir")
    return abs_path

def write_file(path: str, text: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)

def append_log(log_path: str, text: str):
    with open(log_path, "a", encoding="utf-8") as f:
        f.write(text + "\n")

def safe_remove(path: str):
    try:
        if path and os.path.exists(path):
            os.remove(path)
    except Exception:
        pass

def copy_file(src: str, dst: str) -> bool:
    try:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        shutil.copy2(src, dst)
        return True
    except Exception:
        return False

def find_prepared_ff(uid: str, prot_abs: str) -> Optional[dict]:
    base = os.path.join(DATA_DIR, uid, "md")
    try:
        if not os.path.isdir(base):
            return None
        candidates = []
        for name in os.listdir(base):
            d = os.path.join(base, name)
            if not os.path.isdir(d):
                continue
            inp = os.path.join(d, "inputs", "inputs.json")
            out_work = os.path.join(d, "outputs", "work")
            top = os.path.join(out_work, "topol.top")
            gro = os.path.join(out_work, "processed.gro")
            logp = os.path.join(d, "outputs", "log.txt")
            if not (os.path.isfile(inp) and os.path.isfile(top) and os.path.isfile(gro) and os.path.isfile(logp)):
                continue
            try:
                with open(inp, "r", encoding="utf-8") as f:
                    j = json.load(f)
                p = j.get("protein", {}).get("path")
                logtxt = ""
                try:
                    with open(logp, "r", encoding="utf-8", errors="ignore") as lf:
                        logtxt = lf.read()
                except Exception:
                    logtxt = ""
                is_prepare = ("prepare_start" in logtxt) and ("md_start" not in logtxt)
                if p and os.path.abspath(p) == os.path.abspath(prot_abs) and is_prepare:
                    mtime = os.path.getmtime(top)
                    candidates.append({"job": name, "work": out_work, "mtime": mtime})
            except Exception:
                pass
        if not candidates:
            return None
        candidates.sort(key=lambda x: x["mtime"], reverse=True)
        return candidates[0]
    except Exception:
        return None

def sanitize_top_remove_sol(top_path: str):
    try:
        if not os.path.isfile(top_path):
            return
        with open(top_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.read().splitlines()
        idx = None
        for i, ln in enumerate(lines):
            if ln.strip().lower().startswith("[ molecules ]"):
                idx = i
                break
        if idx is None:
            return
        head = lines[:idx+1]
        tail = []
        j = idx + 1
        while j < len(lines):
            s = lines[j].strip()
            if s.startswith("[") and "]" in s:
                tail = lines[j:]
                break
            keep = True
            if s and not s.startswith(";") and not s.startswith("#"):
                first = re.split(r"\s+", s)[0]
                if first.upper() == "SOL":
                    keep = False
            if keep:
                head.append(lines[j])
            j += 1
        text = "\n".join(head + tail)
        if not text.endswith("\n"):
            text += "\n"
        write_file(top_path, text)
    except Exception:
        pass

def merge_gro_with_ligand(prot_gro: str, lig_pdb: str, out_gro: str, work_dir: str, log_path: str) -> bool:
    try:
        lig_gro = os.path.join(work_dir, "lig_amber.gro")
        run_cmd(["gmx", "editconf", "-f", lig_pdb, "-o", lig_gro], work_dir, log_path)
        with open(prot_gro, "r", encoding="utf-8", errors="ignore") as f:
            prot_lines = f.read().splitlines()
        with open(lig_gro, "r", encoding="utf-8", errors="ignore") as f:
            lig_lines = f.read().splitlines()
        prot_title = prot_lines[0] if prot_lines else ""
        prot_n = int(prot_lines[1].strip()) if len(prot_lines) > 1 else 0
        prot_atoms = prot_lines[2:-1]
        box_line = prot_lines[-1] if prot_lines else "0 0 0"
        lig_n = int(lig_lines[1].strip()) if len(lig_lines) > 1 else 0
        lig_atoms = lig_lines[2:-1]
        def renumber(line: str, resid_off: int, atom_off: int) -> str:
            try:
                resid = int(line[0:5])
                atomnr = int(line[15:20])
                new_resid = f"{resid + resid_off:5d}"
                new_atomnr = f"{atomnr + atom_off:5d}"
                return new_resid + line[5:15] + new_atomnr + line[20:]
            except Exception:
                return line
        # Compute residue offset as max resid in protein
        resid_max = 0
        for ln in prot_atoms:
            try:
                resid_max = max(resid_max, int(ln[0:5]))
            except Exception:
                pass
        renum_lig = [renumber(ln, resid_max, prot_n) for ln in lig_atoms]
        total_n = prot_n + lig_n
        out = [prot_title, f"{total_n:5d}"] + prot_atoms + renum_lig + [box_line]
        write_file(out_gro, "\n".join(out))
        append_log(log_path, "merge_coords:protein+ligand")
        return True
    except Exception:
        return False

def run_cmd(cmd: list, cwd: Optional[str], log_path: str, inp: Optional[str] = None, extra_log: Optional[str] = None):
    try:
        p = subprocess.run(cmd, cwd=cwd, input=inp, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        append_log(log_path, " ".join(cmd))
        append_log(log_path, p.stdout or "")
        if extra_log:
            try:
                with open(extra_log, "a", encoding="utf-8") as ef:
                    ef.write(" ".join(cmd) + "\n")
                    ef.write((p.stdout or "") + "\n")
            except Exception:
                pass
        if p.returncode != 0:
            raise RuntimeError("cmd_failed")
    except Exception as e:
        append_log(log_path, f"error:{e}")
        raise HTTPException(status_code=500, detail="md_failed")

def has_lig(path: str) -> bool:
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                if ln.startswith(("ATOM", "HETATM")):
                    res = ln[17:20].strip() if len(ln) >= 21 else ""
                    if res == "LIG" or re.search(r"\bLIG\b", ln):
                        return True
        return False
    except Exception:
        return False

def split_complex_pdb(src: str, lig_out: str, prot_out: str):
    lig_ids = set()
    try:
        with open(src, "r", encoding="utf-8", errors="ignore") as fi, open(lig_out, "w", encoding="utf-8") as fl, open(prot_out, "w", encoding="utf-8") as fp:
            for ln in fi:
                if ln.startswith(("ATOM", "HETATM")):
                    res = ln[17:20].strip() if len(ln) >= 21 else ""
                    if res == "LIG" or re.search(r"\bLIG\b", ln):
                        fl.write(ln)
                        try:
                            sid = int(ln[6:11].strip()) if len(ln) >= 12 else None
                            if sid is not None:
                                lig_ids.add(sid)
                        except Exception:
                            pass
                    else:
                        fp.write(ln)
                elif ln.startswith("CONECT"):
                    parts = re.split(r"\s+", ln.strip())
                    nums = []
                    for p in parts[1:]:
                        try:
                            nums.append(int(p))
                        except Exception:
                            pass
                    if any(n in lig_ids for n in nums):
                        fl.write(ln)
                    else:
                        fp.write(ln)
                else:
                    fp.write(ln)
    except Exception:
        raise

def make_ligand_itp(work_dir: str, log_path: str) -> Optional[tuple]:
    lig_pdb = os.path.join(work_dir, "lig.pdb")
    lig_mol2 = os.path.join(work_dir, "lig.mol2")
    lig_frcmod = os.path.join(work_dir, "lig.frcmod")
    lig_prmtop = os.path.join(work_dir, "lig.prmtop")
    lig_inpcrd = os.path.join(work_dir, "lig.inpcrd")
    lig_top = os.path.join(work_dir, "ligand.top")
    lig_itp = os.path.join(work_dir, "ligand.itp")
    try:
        run_cmd(["antechamber", "-i", lig_pdb, "-fi", "pdb", "-o", lig_mol2, "-fo", "mol2", "-at", "gaff2", "-c", "bcc", "-pf", "y"], work_dir, log_path)
        run_cmd(["parmchk2", "-i", lig_mol2, "-f", "mol2", "-o", lig_frcmod], work_dir, log_path)
        leap_in = os.path.join(work_dir, "tleap.in")
        txt = "\n".join([
            "source leaprc.gaff2",
            f"LIG = loadmol2 {os.path.basename(lig_mol2)}",
            f"loadamberparams {os.path.basename(lig_frcmod)}",
            f"saveamberparm LIG {os.path.basename(lig_prmtop)} {os.path.basename(lig_inpcrd)}",
            f"savepdb LIG lig_amber.pdb",
            "quit",
        ])
        write_file(leap_in, txt)
        run_cmd(["tleap", "-f", "tleap.in"], work_dir, log_path)
        py = os.path.join(work_dir, "parmed_convert.py")
        pytxt = "\n".join([
            "import parmed as pmd",
            f"parm = pmd.load_file('{lig_prmtop}', xyz='{lig_inpcrd}')",
            "ok = False",
            "try:",
            f"    parm.save('{lig_top}', format='gromacs')",
            "    ok = True",
            "except Exception:",
            "    try:",
            "        from parmed.gromacs import GromacsTopologyFile",
            "        gt = GromacsTopologyFile.from_structure(parm)",
            f"        gt.write('{lig_top}')",
            "        ok = True",
            "    except Exception:",
            "        raise",
        ])
        write_file(py, pytxt)
        run_cmd(["python", os.path.basename(py)], work_dir, log_path)
        try:
            with open(lig_top, "r", encoding="utf-8", errors="ignore") as f:
                lines = f.read().splitlines()
            keep_sections = {
                "moleculetype", "atoms", "bonds", "pairs", "angles", "dihedrals", "exclusions",
                "constraints", "settles", "virtual_sites2", "virtual_sites3", "virtual_sites4", "position_restraints"
            }
            cur = None
            out_lines = []
            atomtypes_lines = []
            for ln in lines:
                s = ln.strip()
                if s.startswith("#include"):
                    continue
                m = re.match(r"\s*\[\s*(\w+)\s*\]", s.lower())
                if m:
                    name = m.group(1)
                    cur = name if name in keep_sections else ("atomtypes" if name == "atomtypes" else None)
                    if cur == "atomtypes":
                        atomtypes_lines.append(ln)
                    elif cur:
                        out_lines.append(ln)
                    continue
                if cur == "atomtypes":
                    atomtypes_lines.append(ln)
                elif cur:
                    out_lines.append(ln)
            write_file(lig_itp, "\n".join(out_lines))
            atomtypes_itp = os.path.join(work_dir, "ligand_atomtypes.itp")
            if atomtypes_lines:
                write_file(atomtypes_itp, "\n".join(atomtypes_lines))
            else:
                atomtypes_itp = None
            safe_remove(lig_top)
            safe_remove(leap_in)
            safe_remove(py)
            safe_remove(lig_mol2)
            safe_remove(lig_frcmod)
            safe_remove(lig_prmtop)
            safe_remove(lig_inpcrd)
        except Exception:
            return None
        return (lig_itp, atomtypes_itp) if os.path.exists(lig_itp) else None
    except HTTPException:
        return None
    except Exception:
        append_log(log_path, "error:ambertools_missing")
        return None

def get_moleculetype_name(itp_path: str) -> Optional[str]:
    try:
        with open(itp_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.read().splitlines()
        cur = None
        for i, ln in enumerate(lines):
            s = ln.strip()
            m = re.match(r"\s*\[\s*moleculetype\s*\]", s.lower())
            if m:
                cur = "moleculetype"
                continue
            if cur == "moleculetype":
                if not s or s.startswith(";") or s.startswith("#"):
                    continue
                parts = re.split(r"\s+", s)
                if parts and parts[0] and not parts[0].startswith((";","#")):
                    return parts[0]
                break
        return None
    except Exception:
        return None

def merge_top_with_ligand(prot_top: str, lig_itp: str, out_top: str, atomtypes_itp: Optional[str] = None):
    try:
        with open(prot_top, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.read().splitlines()
        inserted = False
        new_lines = []
        for ln in lines:
            new_lines.append(ln)
            if (not inserted) and ln.strip().startswith("#include") and "forcefield" in ln:
                if atomtypes_itp:
                    new_lines.append(f"#include \"{os.path.basename(atomtypes_itp)}\"")
                new_lines.append(f"#include \"{os.path.basename(lig_itp)}\"")
                inserted = True
        if not inserted:
            if atomtypes_itp:
                new_lines.insert(0, f"#include \"{os.path.basename(atomtypes_itp)}\"")
            new_lines.insert(0, f"#include \"{os.path.basename(lig_itp)}\"")
        idx = None
        for i, ln in enumerate(new_lines):
            if ln.strip().lower().startswith("[ molecules ]"):
                idx = i
                break
        molname = get_moleculetype_name(lig_itp) or "LIG"
        if idx is not None:
            start = idx + 1
            end = start
            has_entry = False
            while end < len(new_lines):
                s = new_lines[end].strip()
                if s.startswith("[") and "]" in s:
                    break
                if s and not s.startswith(";") and not s.startswith("#"):
                    first = re.split(r"\s+", s)[0]
                    if first == molname:
                        has_entry = True
                        break
                end += 1
            if not has_entry:
                new_lines.insert(end, f"{molname}\t1")
                new_lines.insert(end + 1, "")
        else:
            new_lines.append("")
            new_lines.append("[ molecules ]")
            new_lines.append(f"{molname}\t1")
            new_lines.append("")
        text = "\n".join(new_lines)
        if not text.endswith("\n"):
            text += "\n"
        write_file(out_top, text)
        return True
    except Exception:
        return False

def parse_energy_index(edr_path: str, term: str, work_dir: str, log_path: str) -> Optional[int]:
    try:
        p = subprocess.run(["gmx", "energy", "-f", edr_path, "-xvg", "none"], cwd=work_dir, input="0\n", stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        append_log(log_path, p.stdout or "")
        lines = (p.stdout or "").splitlines()
        idx = None
        for ln in lines:
            m = re.match(r"\s*(\d+)\s+(.+?)\s*$", ln)
            if m:
                i = int(m.group(1))
                name = m.group(2).strip()
                if term.lower() in name.lower():
                    idx = i
                    break
        return idx
    except Exception:
        return None

def extract_energy_series(edr_path: str, idx: int, out_path: str, work_dir: str, log_path: str) -> Optional[str]:
    try:
        inp = f"{idx}\n0\n"
        p = subprocess.run(["gmx", "energy", "-f", edr_path, "-o", out_path, "-xvg", "none"], cwd=work_dir, input=inp, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        append_log(log_path, p.stdout or "")
        if p.returncode != 0:
            return None
        return os.path.join(work_dir, out_path)
    except Exception:
        return None

def load_xy(path: str) -> Optional[list]:
    try:
        xs, ys = [], []
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                ln = ln.strip()
                if not ln or ln.startswith(('#','@')):
                    continue
                parts = re.split(r"\s+", ln)
                if len(parts) >= 2:
                    xs.append(float(parts[0]))
                    ys.append(float(parts[1]))
        return list(zip(xs, ys))
    except Exception:
        return None

def plot_xy(data: list, title: str, xlabel: str, ylabel: str, png_path: str):
    if not data or np is None:
        return False
    try:
        # Publication Quality Style Settings
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']
        plt.rcParams['axes.linewidth'] = 1.2
        plt.rcParams['xtick.major.width'] = 1.2
        plt.rcParams['ytick.major.width'] = 1.2
        plt.rcParams['xtick.direction'] = 'in'
        plt.rcParams['ytick.direction'] = 'in'
        
        xs = [d[0] for d in data]
        ys = [d[1] for d in data]
        
        # Increase figure size and DPI for high quality
        fig, ax = plt.subplots(figsize=(6, 4.5), dpi=300)
        
        # Plot with a professional color and style
        ax.plot(xs, ys, lw=1.5, color='#2878B5', alpha=0.9)  # Science Blue
        
        # Titles and Labels
        ax.set_title(title, fontsize=14, fontweight='bold', pad=15)
        ax.set_xlabel(xlabel, fontsize=12, fontweight='medium')
        ax.set_ylabel(ylabel, fontsize=12, fontweight='medium')
        
        # Ticks styling
        ax.tick_params(axis='both', which='major', labelsize=10, pad=6)
        
        # Subtle grid
        ax.grid(True, linestyle='--', alpha=0.4, color='gray', linewidth=0.8)
        
        # Clean spines (Top/Right off)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # Scientific notation if > 3 digits (e.g. 1000 -> 10^3)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(-3, 2), useMathText=True)
        ax.yaxis.get_offset_text().set_fontsize(10)
        
        plt.tight_layout()
        plt.savefig(png_path, bbox_inches='tight', dpi=300)
        plt.close()
        return True
    except Exception:
        return False

def mean_std(vals: list) -> tuple:
    if not vals:
        return (None, None)
    m = float(np.mean(vals)) if np is not None else sum(vals)/len(vals)
    if np is not None:
        s = float(np.std(vals))
    else:
        mu = m
        s = math.sqrt(sum((v-mu)**2 for v in vals)/len(vals))
    return (m, s)

def limit_pdb_models(src: str, dst: str, limit: int):
    cnt = 0
    with open(src, "r", encoding="utf-8") as fi, open(dst, "w", encoding="utf-8") as fo:
        model_open = False
        for line in fi:
            if line.startswith("MODEL"):
                if cnt >= limit:
                    break
                cnt += 1
                model_open = True
            fo.write(line)
        fo.flush()

def count_pdb_models(path: str) -> int:
    try:
        cnt = 0
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("MODEL"):
                    cnt += 1
        return cnt if cnt > 0 else 1
    except Exception:
        return 0

def count_xtc_frames(xtc_path: str, work_dir: str, log_path: str) -> Optional[int]:
    try:
        p = subprocess.run(["gmx", "check", "-f", xtc_path], cwd=work_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        out = p.stdout or ""
        append_log(log_path, out)
        m = re.search(r"Read\s+(\d+)\s+frames", out)
        if not m:
            m = re.search(r"Found\s+(\d+)\s+frames", out)
        return int(m.group(1)) if m else None
    except Exception:
        return None

@app.post("/run")
def run(payload: MdRunInput):
    job_id = payload.job_id
    job_dir, inputs_dir, outputs_dir = ensure_job_dirs(job_id, payload.uid)
    with open(os.path.join(inputs_dir, "inputs.json"), "w", encoding="utf-8") as f:
        json.dump(payload.model_dump(), f, ensure_ascii=False)
    log_path = os.path.join(outputs_dir, "log.txt")
    write_file(log_path, "md_start")
    prot_path = ensure_data_path(payload.protein.path)
    ff = payload.forcefield
    water = payload.water
    steps = payload.steps
    work_dir = os.path.join(outputs_dir, "work")
    os.makedirs(work_dir, exist_ok=True)
    processed_gro = os.path.join(work_dir, "processed.gro")
    topol_top = os.path.join(work_dir, "topol.top")
    if os.path.exists(processed_gro) and os.path.exists(topol_top):
        append_log(log_path, "reuse_topology")
    else:
        src = find_prepared_ff(payload.uid, prot_path)
        if src and os.path.isdir(src.get("work")):
            sw = src["work"]
            ok1 = copy_file(os.path.join(sw, "processed.gro"), processed_gro)
            ok2 = copy_file(os.path.join(sw, "topol.top"), topol_top)
            lig_itp_src = os.path.join(sw, "ligand.itp")
            lig_types_src = os.path.join(sw, "ligand_atomtypes.itp")
            if os.path.isfile(lig_itp_src):
                copy_file(lig_itp_src, os.path.join(work_dir, "ligand.itp"))
            if os.path.isfile(lig_types_src):
                copy_file(lig_types_src, os.path.join(work_dir, "ligand_atomtypes.itp"))
            if ok1 and ok2:
                append_log(log_path, f"reuse_topology_from:{src['job']}")
            else:
                append_log(log_path, "warn:copy_prepare_failed")
        if not (os.path.exists(processed_gro) and os.path.exists(topol_top)):
            try:
                _prep = MdPrepareInput(job_id=payload.job_id, uid=payload.uid, protein=payload.protein, forcefield=ff, water=water)
                _ = prepare(_prep)
                append_log(log_path, "prepare_auto_done")
            except HTTPException:
                append_log("Failed MdPrepareInput!")
        if not (os.path.exists(processed_gro) and os.path.exists(topol_top)):
            append_log(log_path, "Failed")
            raise ValueError("Failed to prepare topology")
    # strip pre-existing SOL entries before solvate to avoid duplication
    sanitize_top_remove_sol(topol_top)
    box_gro = os.path.join(work_dir, "newbox.gro")
    append_log(log_path, "editconf_start")
    run_cmd(["gmx", "editconf", "-f", processed_gro, "-o", box_gro, "-c", "-d", str(steps.box_dist), "-bt", steps.box_type], work_dir, log_path)
    append_log(log_path, "editconf_done")
    solv_gro = os.path.join(work_dir, "solv.gro")
    append_log(log_path, "solvate_start")
    run_cmd(["gmx", "solvate", "-cp", box_gro, "-cs", "spc216.gro", "-o", solv_gro, "-p", topol_top], work_dir, log_path)
    append_log(log_path, "solvate_done")
    ions_mdp = os.path.join(work_dir, "ions.mdp")
    write_file(ions_mdp, "integrator = steep\nemtol = %g\nemstep = 0.01\nnsteps = 2000\nconstraints = none\nnstlog = 100\nnstenergy = 100\n" % steps.em_tol)
    ions_tpr = os.path.join(work_dir, "ions.tpr")
    append_log(log_path, "ions_grompp_start")
    run_cmd(["gmx", "grompp", "-f", ions_mdp, "-c", solv_gro, "-p", topol_top, "-o", ions_tpr], work_dir, log_path)
    append_log(log_path, "ions_grompp_done")
    solv_ions_gro = os.path.join(work_dir, "solv_ions.gro")
    if steps.ion_neutral:
        ok = True
        try:
            append_log(log_path, "genion_start_SOL")
            run_cmd(["gmx", "genion", "-s", ions_tpr, "-o", solv_ions_gro, "-p", topol_top, "-pname", steps.ion_pname, "-nname", steps.ion_nname, "-neutral"], work_dir, log_path, inp="SOL\nSOL\n")
            append_log(log_path, "genion_done_SOL")
        except HTTPException:
            ok = False
        if not ok:
            append_log(log_path, "genion_start_Water")
            run_cmd(["gmx", "genion", "-s", ions_tpr, "-o", solv_ions_gro, "-p", topol_top, "-pname", steps.ion_pname, "-nname", steps.ion_nname, "-neutral"], work_dir, log_path, inp="Water\nWater\n")
            append_log(log_path, "genion_done_Water")
    else:
        solv_ions_gro = solv_gro
    em_mdp = os.path.join(work_dir, "em.mdp")
    npt_mdp = os.path.join(work_dir, "npt.mdp")
    md_mdp = os.path.join(work_dir, "md.mdp")
    write_file(em_mdp, "integrator = steep\nemtol = %g\nemstep = 0.01\nnsteps = %d\nconstraints = none\nnstlog = 100\nnstenergy = 100\n" % (steps.em_tol, steps.em_steps))
    write_file(npt_mdp, "integrator = md\ncontinuation = no\ndt = 0.002\nnsteps = %d\nnstxout-compressed = 100\nnstlog = 100\nnstenergy = 100\ntcoupl = V-rescale\ntc-grps = System\ntau_t = 0.1\nref_t = %g\npcoupl = Berendsen\npcoupltype = isotropic\ntau_p = 2.0\nref_p = 1.0\ncompressibility = 4.5e-5\nconstraints = h-bonds\n" % (steps.npt_steps, steps.npt_temp))
    prod_nstx = int(max(100, math.ceil(steps.prod_steps / 100)))
    if steps.prod_press and steps.prod_press > 0.0:
        md_mdp_txt = "integrator = md\ncontinuation = yes\ndt = 0.002\nnsteps = %d\nnstxout-compressed = %d\nnstlog = 100\nnstenergy = 100\ntcoupl = V-rescale\ntc-grps = System\ntau_t = 0.1\nref_t = %g\npcoupl = Parrinello-Rahman\npcoupltype = isotropic\ntau_p = 2.0\nref_p = %g\ncompressibility = 4.5e-5\nconstraints = h-bonds\n" % (steps.prod_steps, prod_nstx, steps.prod_temp, steps.prod_press)
    else:
        md_mdp_txt = "integrator = md\ncontinuation = yes\ndt = 0.002\nnsteps = %d\nnstxout-compressed = %d\nnstlog = 100\nnstenergy = 100\ntcoupl = V-rescale\ntc-grps = System\ntau_t = 0.1\nref_t = %g\npcoupl = no\nconstraints = h-bonds\n" % (steps.prod_steps, prod_nstx, steps.prod_temp)
    write_file(md_mdp, md_mdp_txt)
    em_tpr = os.path.join(work_dir, "em.tpr")
    em_detail = os.path.join(outputs_dir, "em.log")
    append_log(log_path, "step_start:EM log=em.log")
    append_log(log_path, "em_grompp_start")
    run_cmd(["gmx", "grompp", "-f", em_mdp, "-c", solv_ions_gro, "-p", topol_top, "-o", em_tpr], work_dir, log_path, extra_log=em_detail)
    append_log(log_path, "em_grompp_done")
    append_log(log_path, "em_mdrun_start")
    run_cmd(["gmx", "mdrun", "-deffnm", "em"], work_dir, log_path, extra_log=em_detail)
    append_log(log_path, "em_mdrun_done")
    append_log(log_path, "step_done:EM")
    npt_tpr = os.path.join(work_dir, "npt.tpr")
    npt_detail = os.path.join(outputs_dir, "npt.log")
    append_log(log_path, "step_start:NPT log=npt.log")
    append_log(log_path, "npt_grompp_start")
    run_cmd(["gmx", "grompp", "-f", npt_mdp, "-c", os.path.join(work_dir, "em.gro"), "-p", topol_top, "-o", npt_tpr], work_dir, log_path, extra_log=npt_detail)
    append_log(log_path, "npt_grompp_done")
    append_log(log_path, "npt_mdrun_start")
    run_cmd(["gmx", "mdrun", "-deffnm", "npt"], work_dir, log_path, extra_log=npt_detail)
    append_log(log_path, "npt_mdrun_done")
    append_log(log_path, "step_done:NPT")
    md_tpr = os.path.join(work_dir, "md.tpr")
    prod_detail = os.path.join(outputs_dir, "prod.log")
    append_log(log_path, "step_start:Production log=prod.log")
    append_log(log_path, "prod_grompp_start")
    run_cmd(["gmx", "grompp", "-f", md_mdp, "-c", os.path.join(work_dir, "npt.gro"), "-p", topol_top, "-o", md_tpr], work_dir, log_path, extra_log=prod_detail)
    append_log(log_path, "prod_grompp_done")
    append_log(log_path, "prod_mdrun_start")
    run_cmd(["gmx", "mdrun", "-deffnm", "md"], work_dir, log_path, extra_log=prod_detail)
    append_log(log_path, "prod_mdrun_done")
    append_log(log_path, "step_done:Production")
    em_pdb = os.path.join(outputs_dir, "em.pdb")
    npt_pdb = os.path.join(outputs_dir, "npt.pdb")
    md_pdb = os.path.join(outputs_dir, "md.pdb")
    try:
        run_cmd(["gmx", "trjconv", "-f", os.path.join(work_dir, "em.trr"), "-s", em_tpr, "-o", em_pdb], work_dir, log_path, inp="non-Water\n")
    except HTTPException:
        pass
    try:
        run_cmd(["gmx", "trjconv", "-f", os.path.join(work_dir, "npt.xtc"), "-s", npt_tpr, "-o", npt_pdb], work_dir, log_path, inp="non-Water\n")
    except HTTPException:
        pass
    run_cmd(["gmx", "trjconv", "-f", os.path.join(work_dir, "md.xtc"), "-s", md_tpr, "-o", md_pdb], work_dir, log_path, inp="non-Water\n")
    # use full production frames; nstxout-compressed is set to cap frames to <=100
    # Analysis: EM potential vs step
    em_png = os.path.join(outputs_dir, "em_potential.png")
    em_potential_val = None
    try:
        em_idx = parse_energy_index(os.path.join(work_dir, "em.edr"), "Potential", work_dir, log_path)
        if em_idx is not None:
            em_xvg = extract_energy_series(os.path.join(work_dir, "em.edr"), em_idx, "em_potential.xvg", work_dir, log_path)
            em_data = load_xy(em_xvg) if em_xvg else None
            if em_data:
                plot_xy(em_data, "EM Potential vs Step", "Step", "Potential (kJ/mol)", em_png)
                try:
                    em_potential_val = em_data[-1][1]
                except Exception:
                    em_potential_val = None
                if os.path.exists(em_png):
                    append_log(log_path, f"plot_saved: {os.path.basename(em_png)}")
    except Exception:
        pass
    # NPT: Temperature & Pressure vs time
    npt_temp_png = os.path.join(outputs_dir, "npt_temp.png")
    npt_press_png = os.path.join(outputs_dir, "npt_pressure.png")
    npt_temp_mean, npt_temp_std = (None, None)
    npt_press_mean, npt_press_std = (None, None)
    try:
        npt_temp_idx = parse_energy_index(os.path.join(work_dir, "npt.edr"), "Temperature", work_dir, log_path)
        if npt_temp_idx is not None:
            xvg = extract_energy_series(os.path.join(work_dir, "npt.edr"), npt_temp_idx, "npt_temp.xvg", work_dir, log_path)
            data = load_xy(xvg) if xvg else None
            if data:
                plot_xy(data, "NPT Temperature vs Time", "Time (ps)", "Temperature (K)", npt_temp_png)
                npt_temp_mean, npt_temp_std = mean_std([y for _, y in data])
                if os.path.exists(npt_temp_png):
                    append_log(log_path, f"plot_saved: {os.path.basename(npt_temp_png)}")
        npt_press_idx = parse_energy_index(os.path.join(work_dir, "npt.edr"), "Pressure", work_dir, log_path)
        if npt_press_idx is not None:
            xvg = extract_energy_series(os.path.join(work_dir, "npt.edr"), npt_press_idx, "npt_pressure.xvg", work_dir, log_path)
            data = load_xy(xvg) if xvg else None
            if data:
                plot_xy(data, "NPT Pressure vs Time", "Time (ps)", "Pressure (bar)", npt_press_png)
                npt_press_mean, npt_press_std = mean_std([y for _, y in data])
                if os.path.exists(npt_press_png):
                    append_log(log_path, f"plot_saved: {os.path.basename(npt_press_png)}")
    except Exception:
        pass
    # Production: Temperature & Pressure vs time
    prod_temp_png = os.path.join(outputs_dir, "prod_temp.png")
    prod_press_png = os.path.join(outputs_dir, "prod_pressure.png")
    prod_temp_mean, prod_temp_std = (None, None)
    prod_press_mean, prod_press_std = (None, None)
    try:
        md_temp_idx = parse_energy_index(os.path.join(work_dir, "md.edr"), "Temperature", work_dir, log_path)
        if md_temp_idx is not None:
            xvg = extract_energy_series(os.path.join(work_dir, "md.edr"), md_temp_idx, "md_temp.xvg", work_dir, log_path)
            data = load_xy(xvg) if xvg else None
            if data:
                plot_xy(data, "Production Temperature vs Time", "Time (ps)", "Temperature (K)", prod_temp_png)
                prod_temp_mean, prod_temp_std = mean_std([y for _, y in data])
                if os.path.exists(prod_temp_png):
                    append_log(log_path, f"plot_saved: {os.path.basename(prod_temp_png)}")
        md_press_idx = parse_energy_index(os.path.join(work_dir, "md.edr"), "Pressure", work_dir, log_path)
        if md_press_idx is not None:
            xvg = extract_energy_series(os.path.join(work_dir, "md.edr"), md_press_idx, "md_pressure.xvg", work_dir, log_path)
            data = load_xy(xvg) if xvg else None
            if data:
                plot_xy(data, "Production Pressure vs Time", "Time (ps)", "Pressure (bar)", prod_press_png)
                prod_press_mean, prod_press_std = mean_std([y for _, y in data])
                if os.path.exists(prod_press_png):
                    append_log(log_path, f"plot_saved: {os.path.basename(prod_press_png)}")
    except Exception:
        pass
    # Production: Potential & Kinetic xvg for analysis
    try:
        pot_idx = parse_energy_index(os.path.join(work_dir, "md.edr"), "Potential", work_dir, log_path)
        if pot_idx is not None:
            _x = extract_energy_series(os.path.join(work_dir, "md.edr"), pot_idx, "md_potential.xvg", work_dir, log_path)
            if _x:
                append_log(log_path, "saved_xvg: md_potential.xvg")
        kin_idx = parse_energy_index(os.path.join(work_dir, "md.edr"), "Kinetic", work_dir, log_path)
        if kin_idx is not None:
            _x = extract_energy_series(os.path.join(work_dir, "md.edr"), kin_idx, "md_kinetic.xvg", work_dir, log_path)
            if _x:
                append_log(log_path, "saved_xvg: md_kinetic.xvg")
    except Exception:
        pass
    # Report
    report_csv = os.path.join(outputs_dir, "report.csv")
    lines = ["step,frames,potential,temp_mean,temp_error,press_mean,press_error"]
    em_pot_str = ("%0.2f" % em_potential_val) if em_potential_val is not None else "-"
    lines.append(f"EM,1,{em_pot_str},,,,")
    t_mean = ("%0.2f" % npt_temp_mean) if npt_temp_mean is not None else ""
    t_err = ("%0.2f" % npt_temp_std) if npt_temp_std is not None else ""
    p_mean = ("%0.3f" % npt_press_mean) if npt_press_mean is not None else ""
    p_err = ("%0.3f" % npt_press_std) if npt_press_std is not None else ""
    lines.append(f"NPT,1,-,{t_mean},{t_err},{p_mean},{p_err}")
    pt_mean = ("%0.2f" % prod_temp_mean) if prod_temp_mean is not None else ""
    pt_err = ("%0.2f" % prod_temp_std) if prod_temp_std is not None else ""
    pp_mean = ("%0.3f" % prod_press_mean) if prod_press_mean is not None else ""
    pp_err = ("%0.3f" % prod_press_std) if prod_press_std is not None else ""
    prod_frames = count_xtc_frames(os.path.join(work_dir, "md.xtc"), work_dir, log_path)
    if prod_frames is None:
        prod_frames = count_pdb_models(md_pdb)
    lines.append(f"Production,{prod_frames},-,{pt_mean},{pt_err},{pp_mean},{pp_err}")
    write_file(report_csv, "\n".join(lines))
    outputs_dict = {
        "log": "log.txt",
        "report": "report.csv",
        "traj_pdb": "md.pdb",
        "em_log": ("em.log" if os.path.exists(em_detail) else None),
        "npt_log": ("npt.log" if os.path.exists(npt_detail) else None),
        "prod_log": ("prod.log" if os.path.exists(prod_detail) else None),
        "em_plot": ("em_potential.png" if os.path.exists(em_png) else None),
        "npt_temp_plot": ("npt_temp.png" if os.path.exists(npt_temp_png) else None),
        "npt_press_plot": ("npt_pressure.png" if os.path.exists(npt_press_png) else None),
        "prod_temp_plot": ("prod_temp.png" if os.path.exists(prod_temp_png) else None),
        "prod_press_plot": ("prod_pressure.png" if os.path.exists(prod_press_png) else None),
    }
    try:
        append_log(log_path, "outputs_keys:" + ",".join([k for k,v in outputs_dict.items() if v]))
    except Exception:
        pass
    return {"job_id": job_id, "outputs": outputs_dict}

class MdPrepareInput(BaseModel):
    job_id: str
    uid: str
    protein: MoleculeRef
    forcefield: str = Field(default="amber99sb")
    water: str = Field(default="tip3p")

@app.post("/prepare")
def prepare(payload: MdPrepareInput):
    job_id = payload.job_id
    job_dir, inputs_dir, outputs_dir = ensure_job_dirs(job_id, payload.uid)
    with open(os.path.join(inputs_dir, "inputs.json"), "w", encoding="utf-8") as f:
        json.dump(payload.model_dump(), f, ensure_ascii=False)
    log_path = os.path.join(outputs_dir, "log.txt")
    write_file(log_path, "prepare_start")
    prot_path = ensure_data_path(payload.protein.path)
    ff = payload.forcefield
    water = payload.water
    work_dir = os.path.join(outputs_dir, "work")
    os.makedirs(work_dir, exist_ok=True)
    if has_lig(prot_path):
        append_log(log_path, "detect_ligand:LIG")
        lig_pdb = os.path.join(work_dir, "lig.pdb")
        prot_pdb = os.path.join(work_dir, "protein.pdb")
        split_complex_pdb(prot_path, lig_pdb, prot_pdb)
        processed_gro = os.path.join(work_dir, "processed.gro")
        topol_top_prot = os.path.join(work_dir, "topol_protein.top")
        run_cmd(["gmx", "pdb2gmx", "-f", prot_pdb, "-o", processed_gro, "-p", topol_top_prot, "-ff", ff, "-water", water, "-ignh", "--missing"], work_dir, log_path)
        res = make_ligand_itp(work_dir, log_path)
        if not res:
            raise HTTPException(status_code=500, detail="md_failed")
        lig_itp, atomtypes_itp = res
        topol_top = os.path.join(work_dir, "topol.top")
        ok = merge_top_with_ligand(topol_top_prot, lig_itp, topol_top, atomtypes_itp)
        if not ok:
            raise HTTPException(status_code=500, detail="md_failed")
        try:
            if os.path.exists(topol_top):
                append_log(log_path, "merge_top_saved:" + os.path.basename(topol_top))
        except Exception:
            pass
        # Merge ligand coordinates into processed_gro using Amber-generated PDB
        lig_amber_pdb = os.path.join(work_dir, "lig_amber.pdb")
        ok_coords = merge_gro_with_ligand(processed_gro, lig_amber_pdb, processed_gro, work_dir, log_path)
        if not ok_coords:
            append_log(log_path, "warn:merge_coords_failed")
        else:
            safe_remove(lig_amber_pdb)
        posre_itp = os.path.join(work_dir, "posre.itp")
        return {"job_id": job_id, "outputs": {"log": "log.txt", "topol": "work/topol.top", "gro": "work/processed.gro", "posre_itp": ("work/posre.itp" if os.path.exists(posre_itp) else None), "ligand_itp": (f"work/{os.path.basename(lig_itp)}" if os.path.exists(lig_itp) else None), "ligand_atomtypes_itp": ("work/ligand_atomtypes.itp" if (atomtypes_itp and os.path.exists(atomtypes_itp)) else None)}}
    else:
        processed_gro = os.path.join(work_dir, "processed.gro")
        topol_top = os.path.join(work_dir, "topol.top")
        run_cmd(["gmx", "pdb2gmx", "-f", prot_path, "-o", processed_gro, "-p", topol_top, "-ff", ff, "-water", water, "-ignh", "--missing"], work_dir, log_path)
        posre_itp = os.path.join(work_dir, "posre.itp")
        return {"job_id": job_id, "outputs": {"log": "log.txt", "topol": "work/topol.top", "gro": "work/processed.gro", "posre_itp": ("work/posre.itp" if os.path.exists(posre_itp) else None)}}
