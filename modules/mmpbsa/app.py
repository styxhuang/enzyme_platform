import os
import json
import base64
import datetime
from typing import Optional, Union
from fastapi import FastAPI
from pydantic import BaseModel
import csv
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = os.getenv("ENZYME_DATA_DIR", "/data")

class GroupOpts(BaseModel):
    group: str = "System"
    group_index: Optional[str] = None

class GroupSpec(BaseModel):
    name: str = "System"
    index: Optional[str] = None

class ImplicitGbAdv(BaseModel):
    igb: Optional[str] = None
    salt: Optional[str] = None
    inter_diel: Optional[str] = None
    solv_diel: Optional[str] = None

class ImplicitPbAdv(BaseModel):
    npbopt: Optional[str] = None
    ipb: Optional[str] = None
    ionic_s: Optional[str] = None
    inter_diel: Optional[str] = None
    solv_diel: Optional[str] = None

class RismAdv(BaseModel):
    closure: Optional[str] = None  # e.g., 'kh', 'hnc'
    buffer: Optional[float] = None # e.g., 12.0
    solvbox: Optional[float] = None # e.g., 30.0
    griddim: Optional[int] = None  # e.g., 128
    temperature: Optional[float] = None # Kelvin

class AlanineSettings(BaseModel):
    range: Optional[str] = "manual"
    interface_cutoff: Optional[float] = 10.0
    probe_radius: Optional[float] = 5.0
    residues: Optional[str] = None

class MmpFunctions(BaseModel):
    binding: bool = False
    decomp: bool = False
    alanine: bool = False

class MmpSettings(BaseModel):
    bind_step: Optional[str] = None
    decomp_gran: Optional[str] = None
    # ala_res and scan_range moved to alanine nested model, kept here for backward compatibility if needed, but we will prefer nested
    ala_res: Optional[str] = None
    scan_range: Optional[str] = "manual"
    interface_cutoff: Optional[float] = 10.0
    
    alanine: Optional[AlanineSettings] = None
    functions: Optional[MmpFunctions] = None
    
    implicit_model: Optional[str] = "GB"
    gb_adv: Optional[ImplicitGbAdv] = None
    pb_adv: Optional[ImplicitPbAdv] = None
    rism_adv: Optional[RismAdv] = None
    groups: Optional[dict] = None

class MmpRunInput(BaseModel):
    job_id: str
    uid: str
    source_job_id: str
    group: Union[GroupOpts, str]
    mode: str
    settings: MmpSettings

app = FastAPI()

def ensure_job_dirs(job_id: str, uid: str):
    job_dir = os.path.join(DATA_DIR, uid, "mmpbsa", job_id)
    inputs_dir = os.path.join(job_dir, "inputs")
    outputs_dir = os.path.join(job_dir, "outputs")
    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(outputs_dir, exist_ok=True)
    return job_dir, inputs_dir, outputs_dir

def write_text(path: str, text: str):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)

def append_log(path: str, text: str):
    with open(path, "a", encoding="utf-8") as f:
        f.write(text + "\n")

def which(cmd: str) -> Optional[str]:
    for p in os.environ.get("PATH", "").split(":"):
        fp = os.path.join(p, cmd)
        if os.path.isfile(fp) and os.access(fp, os.X_OK):
            return fp
    return None

def find_gmx_mmpbsa_cmd() -> Optional[list]:
    p = which("gmx_MMPBSA")
    if p:
        return [p]
    for conda in ("/opt/conda/bin/conda", "/root/miniconda3/bin/conda", which("conda")):
        if conda and os.path.isfile(conda):
            return [conda, "run", "-n", os.getenv("MMPBSA_CONDA_ENV", "mmpbsa"), "gmx_MMPBSA"]
    return None

def read_index_groups(ndx_path: str) -> dict:
    groups = {}
    try:
        with open(ndx_path, "r", encoding="utf-8", errors="ignore") as f:
            idx = 0
            for ln in f:
                s = ln.strip()
                if s.startswith("[") and s.endswith("]"):
                    name = s[1:-1].strip()
                    groups[name] = idx
                    idx += 1
        return groups
    except Exception:
        return {}

def _detect_columns(header: list) -> dict:
    idx = {}
    lower = [h.lower() for h in header]
    for name, keys in {
        'vdw': ['vdw','vdwaals','van'],
        'elec': ['elec','eel','coul'],
        'gb': ['gb','egb','gbpolar'],
        'sa': ['sa','surf','esurf','apolar'],
        'total': ['total','tot']
    }.items():
        found = None
        for k in keys:
            for i,h in enumerate(lower):
                if k in h:
                    found = i
                    break
            if found is not None:
                break
        if found is not None:
            idx[name] = found
    # time/frame column
    for i,h in enumerate(lower):
        if 'time' in h or 'frame' in h:
            idx['time'] = i
            break
    return idx

def _write_interaction(outputs_dir: str, rows: list, idxs: dict):
    times = []
    vdw = []
    elec = []
    total = []
    time_idx = idxs.get('time', 0)
    for r in rows:
        try:
            t = r[time_idx]
            times.append(float(t))
        except:
            times.append(len(times))
        try:
            vdw.append(float(r[idxs['vdw']]))
        except:
            vdw.append(float('nan'))
        try:
            elec.append(float(r[idxs['elec']]))
        except:
            elec.append(float('nan'))
        try:
            total.append(float(r[idxs['total']]))
        except:
            total.append(float('nan'))
    csv_path = os.path.join(outputs_dir, 'interaction_PL.csv')
    with open(csv_path,'w',encoding='utf-8') as f:
        f.write('time,vdw,elec,total\n')
        for i in range(len(times)):
            f.write(f"{times[i]},{vdw[i]},{elec[i]},{total[i]}\n")
    def m_s(vals):
        vals2 = [x for x in vals if not math.isnan(x)]
        if not vals2:
            return ('','')
        return (f"{float(np.mean(vals2)):.4f}", f"{float(np.std(vals2)):.4f}")
    mv, sv = m_s(vdw)
    me, se = m_s(elec)
    mt, st = m_s(total)
    sum_path = os.path.join(outputs_dir, 'interaction_PL_summary.csv')
    with open(sum_path,'w',encoding='utf-8') as f:
        f.write('component,mean,error\n')
        f.write(f"P-L VdW. (kcal/mol),{mv},{sv}\n")
        f.write(f"P-L Elec. (kcal/mol),{me},{se}\n")
        f.write(f"P-L Total. (kcal/mol),{mt},{st}\n")
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial','Helvetica','DejaVu Sans']
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rcParams['xtick.major.width'] = 1.2
    plt.rcParams['ytick.major.width'] = 1.2
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    try:
        fig, ax = plt.subplots(figsize=(6,4.5), dpi=300)
        plotted = False
        if any(not math.isnan(x) for x in vdw):
            ax.plot(times, vdw, lw=1.5, label='VdW', color='#2878B5', alpha=0.9)
            plotted = True
        if any(not math.isnan(x) for x in elec):
            ax.plot(times, elec, lw=1.5, label='Elec', color='#D76364', alpha=0.9)
            plotted = True
        if any(not math.isnan(x) for x in total):
            ax.plot(times, total, lw=1.5, label='Total', color='#53A859', alpha=0.9)
            plotted = True
        if not plotted:
            # fallback: plot Total summary value if available
            t_avg, _t_sd = ('', '')
            try:
                valid_total = [x for x in total if not math.isnan(x)]
                t_avg = None if not valid_total else float(np.mean(valid_total))
            except Exception:
                t_avg = None
            if t_avg is None:
                t_avg = 0.0
            xs = times if len(times) > 0 else [1]
            ys = [t_avg] if len(xs) == 1 else [t_avg]*len(xs)
            ax.plot(xs, ys, lw=1.5, label='Total', color='#53A859', alpha=0.9)
        ax.set_title('Interaction Energy', fontsize=14, fontweight='bold', pad=15)
        ax.set_xlabel('Time (ps)', fontsize=12, fontweight='medium')
        ax.set_ylabel('Energy (kcal/mol)', fontsize=12, fontweight='medium')
        ax.grid(True, linestyle='--', alpha=0.4, color='gray', linewidth=0.8)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.ticklabel_format(style='sci', axis='y', scilimits=(-3,2), useMathText=True)
        ax.legend(loc='best', fontsize=10, frameon=False)
        plt.tight_layout()
        png_path = os.path.join(outputs_dir, 'interaction_PL.png')
        plt.savefig(png_path, bbox_inches='tight', dpi=300)
        plt.close()
    except Exception:
        try:
            fig = plt.figure(figsize=(6,4.5), dpi=300)
            plt.plot([1],[0.0], lw=1.5, label='Total', color='#53A859', alpha=0.9)
            plt.legend(loc='best', fontsize=10, frameon=False)
            plt.tight_layout()
            png_path = os.path.join(outputs_dir, 'interaction_PL.png')
            plt.savefig(png_path, bbox_inches='tight', dpi=300)
            plt.close()
        except Exception:
            pass
    return {
        'interaction_csv': 'interaction_PL.csv',
        'interaction_summary': 'interaction_PL_summary.csv',
        'interaction_plot': 'interaction_PL.png'
    }

def write_mmpbsa_in(path: str, startframe: int = 1, endframe: int = 10, model: str = "GB", igb: int = 5, saltcon: float = 0.150, pb_opts: Optional[ImplicitPbAdv] = None, rism_opts: Optional[RismAdv] = None, include_decomp: bool = False, ala_res: Optional[str] = None):
    lines = []
    lines.append("&general")
    lines.append("sys_name=\"Prot-Lig-ST\",")
    lines.append(f"startframe={startframe},")
    lines.append(f"endframe={endframe},")
    if ala_res:
        lines.append("netcdf=0,")
    lines.append("/")
    m = (model or "GB").upper()
    if m == "GB":
        lines.append("&gb")
        lines.append(f"igb={igb}, saltcon={saltcon},")
        lines.append("/")
    elif m == "PB":
        lines.append("&pb")
        if pb_opts and pb_opts.ipb:
            lines.append(f"ipb={pb_opts.ipb},")
        if pb_opts and pb_opts.ionic_s:
            lines.append(f"ionic_s={pb_opts.ionic_s},")
        if pb_opts and pb_opts.inter_diel:
            lines.append(f"inpr={pb_opts.inter_diel},")
        if pb_opts and pb_opts.solv_diel:
            lines.append(f"exdi={pb_opts.solv_diel},")
        lines.append("/")
    elif m in ("RISM", "3D-RISM", "3DRISM"):
        lines.append("&rism")
        closure = (rism_opts.closure if rism_opts and rism_opts.closure else "kh")
        buffer = (rism_opts.buffer if rism_opts and rism_opts.buffer is not None else 12.0)
        solvbox = (rism_opts.solvbox if rism_opts and rism_opts.solvbox is not None else 30.0)
        griddim = (rism_opts.griddim if rism_opts and rism_opts.griddim is not None else 128)
        lines.append(f"closure=\"{closure}\",")
        lines.append(f"buffer={buffer},")
        lines.append(f"solvbox={solvbox},")
        lines.append(f"griddim={griddim},")
        lines.append("/")
    if include_decomp:
        lines.append("&decomp")
        lines.append("idecomp=2,")
        lines.append("dec_verbose=3,")
        lines.append("/")
    if ala_res:
        lines.append("&alanine_scanning")
        lines.append(f"mutant_res='{ala_res}',")
        lines.append("/")
        append_log(os.path.join(os.path.dirname(path), "log.txt"), f"Configuring Alanine Scanning for: {ala_res}")
    else:
        # Standard run, no alanine scanning section
        pass
    write_text(path, "\n".join(lines))

def append_ligand_group(ndx_path: str, gro_path: str) -> bool:
    try:
        atoms = []
        with open(gro_path, "r", encoding="utf-8", errors="ignore") as f:
            lines = f.read().splitlines()
        if len(lines) >= 3:
            for ln in lines[2:-1]:
                ln = ln.strip()
                if not ln:
                    continue
                parts = ln.split()
                if len(parts) >= 4:
                    resname = parts[1]
                    atomnr = parts[3]
                    if str(resname).upper() == "LIG":
                        try:
                            atoms.append(int(atomnr))
                        except Exception:
                            pass
        if atoms:
            os.makedirs(os.path.dirname(ndx_path), exist_ok=True)
            with open(ndx_path, "a", encoding="utf-8") as nf:
                nf.write("\n[ Ligand ]\n")
                buf = []
                for i, a in enumerate(atoms, 1):
                    buf.append(str(a))
                    if i % 15 == 0:
                        nf.write(" ".join(buf) + "\n")
                        buf = []
                if buf:
                    nf.write(" ".join(buf) + "\n")
            return True
    except Exception:
        pass
    return False

def create_index_and_groups(gro: str, outputs_dir: str) -> tuple:
    ndx = os.path.join(outputs_dir, "ana_index.ndx")
    if os.path.exists(gro) and which("gmx"):
        import subprocess
        try:
            subprocess.run(["gmx", "make_ndx", "-f", gro, "-o", ndx], input="q\n", cwd=outputs_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        except Exception:
            pass
        try:
            if os.path.exists(ndx):
                with open(ndx, "r", encoding="utf-8", errors="ignore") as f:
                    lines = f.read().splitlines()
                out = []
                for ln in lines:
                    s = ln.strip()
                    if s.startswith("[") and "LIG" in s:
                        ln = ln.replace("LIG", "Ligand")
                    out.append(ln)
                with open(ndx, "w", encoding="utf-8") as f:
                    f.write("\n".join(out))
        except Exception:
            pass
    if os.path.exists(gro):
        append_ligand_group(ndx, gro)
    groups = read_index_groups(ndx) if os.path.exists(ndx) else {}
    rec_idx = groups.get("Protein")
    lig_idx = groups.get("Ligand")
    return ndx if os.path.exists(ndx) else None, rec_idx, lig_idx


def _plot_complex_grouped_bars(rows: list, outputs_dir: str) -> str:
    labels = [r['Residue'] for r in rows]
    vdw_vals = [r.get('VDW', 0.0) for r in rows]
    eel_vals = [r.get('EEL', 0.0) for r in rows]
    ps_vals = [r.get('PS', 0.0) for r in rows]
    tot_vals = [r.get('TOTAL', 0.0) for r in rows]
    if len(labels) > 20:
        labels = labels[:20]
        vdw_vals = vdw_vals[:20]
        eel_vals = eel_vals[:20]
        ps_vals = ps_vals[:20]
        tot_vals = tot_vals[:20]
    x = np.arange(len(labels))
    width = 0.2
    fig, ax = plt.subplots(figsize=(max(8, len(labels)*0.6), 6), dpi=300)
    ax.bar(x - 1.5*width, vdw_vals, width, label='van der Waals', color='#8E5CE6')
    ax.bar(x - 0.5*width, eel_vals, width, label='Electrostatic', color='#F4C542')
    ax.bar(x + 0.5*width, ps_vals, width, label='Polar', color='#45B769')
    ax.bar(x + 1.5*width, tot_vals, width, label='Total', color='#E53E3E')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.set_ylabel('Energy (kcal/mol)', fontsize=12)
    ax.set_title('Complex Total Energy Decomposition', fontsize=14, fontweight='bold')
    ax.grid(True, axis='y', linestyle='--', alpha=0.4)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(fontsize=10, frameon=False)
    plt.tight_layout()
    png_path = os.path.join(outputs_dir, 'decomposition.png')
    plt.savefig(png_path, bbox_inches='tight', dpi=300)
    plt.close()
    return 'decomposition.png'

def _extract_section_rows_from_final_csv(csv_path: str, section_name: str):
    try:
        with open(csv_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
        idx = -1
        for i, ln in enumerate(lines):
            if ln.strip().lower().startswith(section_name.lower()):
                idx = i
                break
        if idx == -1:
            return None
        j = idx + 1
        while j < len(lines) and not lines[j].strip():
            j += 1
        if j >= len(lines):
            return None
        header_line = lines[j].strip()
        header = [h.strip() for h in header_line.split(',')]
        rows = []
        k = j + 1
        while k < len(lines):
            s = lines[k].strip()
            if not s:
                break
            if s.lower().endswith('energy terms'):
                break
            parts = [p.strip() for p in s.split(',')]
            rows.append(parts)
            k += 1
        if not rows:
            return None
        return (header, rows)
    except Exception:
        return None

def _generate_interaction_from_final_csv(csv_path: str, outputs_dir: str) -> dict:
    out = {}
    try:
        sec = _extract_section_rows_from_final_csv(csv_path, 'Delta Energy Terms')
        if not sec:
            sec = _extract_section_rows_from_final_csv(csv_path, 'Complex Energy Terms')
        if not sec:
            return {}
        header, rows = sec
        idxs = _detect_columns(header)
        if all(k in idxs for k in ('vdw','elec','total')):
            inter = _write_interaction(outputs_dir, rows, idxs)
            out.update(inter)
        return out
    except Exception:
        return {}

def _parse_final_decomp_dat(dat_path: str, outputs_dir: str, log_path: str) -> dict:
    out = {}
    try:
        data = []
        import re
        with open(dat_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()
        # Detect Complex: followed by Total Energy Decomposition:
        complex_line = -1
        for i, ln in enumerate(lines):
            if ln.strip().lower().startswith('deltas:'):
                complex_line = i
                break
        if complex_line != -1:
            j = complex_line + 1
            while j < len(lines) and not lines[j].strip():
                j += 1
            if j < len(lines) and lines[j].strip().lower().startswith('total energy decomposition'):
                j += 1
                rows = []
                k = j + 2
                while k < len(lines):
                    sline = lines[k].strip()
                    if not sline:
                        break
                    parts = [p.strip() for p in sline.split(',')]
                    if len(parts) > 19:
                        break
                    label = parts[0]
                    def _get(idx):
                        try:
                            return float(parts[idx])
                        except Exception:
                            return None
                    # 0-based indices for Avg columns in Complex Total table
                    vdw = _get(4) or 0.0
                    vdw_err = _get(5) or 0.0
                    eel = _get(7) or 0.0
                    eel_err = _get(8) or 0.0
                    ps = _get(10) or 0.0
                    ps_err = _get(11) or 0.0
                    tot = _get(16) or 0.0
                    tot_err = _get(17) or 0.0
                    if tot is None:
                        try:
                            tot = float(parts[-3])
                        except Exception:
                            tot = 0.0
                    rows.append({ 'Residue': label, 'VDW': vdw, 'VDW_ERR': vdw_err, 'EEL': eel, 'EEL_ERR': eel_err, 'PS': ps, 'PS_ERR': ps_err, 'TOTAL': tot, 'TOTAL_ERR': tot_err })
                    k += 1
                if rows:
                    data.extend(rows)
            
        csv_name = 'residues_energy_summary.csv'
        csv_path = os.path.join(outputs_dir, csv_name)
        headers = ['Residue']
        for k in ['VDW','EEL','PS','TOTAL']:
            headers.append(k)
        with open(csv_path, 'w', encoding='utf-8') as f:
            f.write(','.join(headers) + '\n')
            for row in data:
                vals = []
                for h in headers:
                    v = row.get(h, '')
                    if isinstance(v, float):
                        vals.append(f"{v:.4f}")
                    else:
                        vals.append(str(v))
                f.write(','.join(vals) + '\n')
        out['decomposition_csv'] = csv_name
        out['decomposition_plot'] = _plot_complex_grouped_bars(rows, outputs_dir)
        return out
    except Exception:
        return {}

def _get_pdb_chain_res_map(pdb_path: str) -> dict:
    """
    Parses a PDB file and returns a dict mapping atom index (1-based) to 'CHAIN:RESNUM'.
    """
    mapping = {}
    try:
        with open(pdb_path, 'r', encoding='utf-8', errors='ignore') as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        # PDB format:
                        # 7-11 Atom serial number
                        # 22 Chain identifier
                        # 23-26 Residue sequence number
                        serial_str = line[6:11].strip()
                        if not serial_str: continue
                        serial = int(serial_str)
                        
                        chain = line[21:22].strip()
                        if not chain:
                            chain = 'A' # Default to A if empty
                            
                        resnum = line[22:26].strip()
                        mapping[serial] = f"{chain}:{resnum}"
                    except:
                        pass
    except Exception:
        pass
    return mapping

def _get_interface_residues(tpr: str, xtc: str, outputs_dir: str, cutoff: float, log_path: str, ndx_path: str = None) -> list:
    """
    Identifies residues in the interface using gmx select.
    Returns a list of unique 'CHAIN:RESNUM' strings.
    """
    import subprocess
    
    gmx_cmd = ["gmx"]

    # 1. Extract PDB (frame 0) to define chains and residues
    pdb_path = os.path.join(outputs_dir, "ref_struct.pdb")
    if not os.path.exists(pdb_path):
        # echo 0 | gmx trjconv -f tpr -s tpr -o pdb_path -dump 0
        cmd = ' '.join(gmx_cmd + ["trjconv", "-f", xtc, "-s", tpr, "-o", pdb_path, "-dump", "0"])
        try:
            subprocess.run(cmd, input="0\n", cwd=outputs_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True, shell=True)
        except Exception as e:
            if log_path:
                append_log(log_path, f"gmx trjconv failed: {e}")
            return []

    # 2. Map atoms to residues
    atom_map = _get_pdb_chain_res_map(pdb_path)
    if not atom_map:
        if log_path:
            append_log(log_path, "Failed to generate atom map from PDB.")
        return []

    # 3. Run gmx select to find interface atoms
    # We use the generated PDB as structure to ensure atom indices match
    cutoff_nm = cutoff / 10.0
    
    # Use provided index file if available (to ensure Ligand group exists)
    ndx_args = ["-n", ndx_path] if ndx_path and os.path.exists(ndx_path) else []
    
    # Select residues within cutoff of Ligand
    # User feedback: assignment syntax ("Name" = ...) causes errors, but raw selection works.
    select_cmd = f"same residue as (group \"Protein\" and within {cutoff_nm} of group \"Ligand\")"
    
    ndx_out = os.path.join(outputs_dir, "interface.ndx")
    if log_path:
        append_log(log_path, f"Running gmx select with command: {select_cmd}")
    
    # Use shell=True to handle quotes correctly
    cmd = ' '.join(gmx_cmd + ["select", "-s", pdb_path, "-on", ndx_out, "-select", f"'{select_cmd}'"] + ndx_args)
    
    try:
        subprocess.run(cmd, cwd=outputs_dir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, check=True, shell=True)
    except Exception as e:
        if log_path:
            append_log(log_path, f"gmx select failed with command: {cmd}. Error: {str(e)}")
        return []
    
    
    # 4. Parse output ndx to get atoms
    residues = set()
    try:
        current_group = None
        if os.path.exists(ndx_out):
            with open(ndx_out, 'r', encoding='utf-8', errors='ignore') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("["):
                        current_group = line.strip("[] ")
                        continue
                    # We accept any group in the output file since we just generated it for this selection
                    if current_group: 
                        parts = line.split()
                        for p in parts:
                            try:
                                idx = int(p)
                                if idx in atom_map:
                                    residues.add(atom_map[idx])
                            except:
                                pass
    except Exception as e:
        if log_path:
            append_log(log_path, f"Error parsing index file: {e}")
        pass
        
    # Sort residues naturally (e.g. A:2, A:10 instead of A:10, A:2)
    def natural_sort_key(s):
        import re
        parts = re.split(r'(\d+)', s)
        return [int(p) if p.isdigit() else p for p in parts]
        
    return sorted(list(residues), key=natural_sort_key)

def _extract_ddg(dat_path: str):
    try:
        if not os.path.exists(dat_path):
            return 0.0, 0.0
        with open(dat_path, 'r', errors='ignore') as f:
            lines = f.readlines()
        for line in lines:
            if "Delta Delta G binding" in line:
                # Format: Delta Delta G binding =   2.5000 +/-  0.1000
                parts = line.split('=')
                if len(parts) > 1:
                    val_part = parts[1].strip()
                    v_e = val_part.split('+/-')
                    val = float(v_e[0].strip())
                    err = float(v_e[1].strip()) if len(v_e) > 1 else 0.0
                    return val, err
    except:
        pass
    return 0.0, 0.0

def _plot_alanine_scanning(data, outputs_dir):
    if not data: return
    try:
        residues = [d['Residue'] for d in data]
        vals = [d['ddG'] for d in data]
        errs = [d['Error'] for d in data]
        
        plt.figure(figsize=(max(8, len(residues)*0.5), 6), dpi=300)
        plt.bar(residues, vals, yerr=errs, capsize=5, color='#E53E3E', alpha=0.8)
        plt.axhline(0, color='black', linewidth=0.8)
        plt.xlabel('Residue')
        plt.ylabel(r'$\Delta\Delta G_{binding}$ (kcal/mol)')
        plt.title('Alanine Scanning Results')
        plt.xticks(rotation=45, ha='right')
        plt.tight_layout()
        plt.savefig(os.path.join(outputs_dir, "alanine_scanning.png"))
        plt.close()
    except:
        pass

def _image_to_base64(path):
    try:
        if not os.path.exists(path):
            return ""
        with open(path, "rb") as f:
            return base64.b64encode(f.read()).decode('utf-8')
    except Exception:
        return ""

def _generate_html_report(outputs_dir: str, payload: MmpRunInput, outputs: dict) -> Optional[str]:
    job_id = payload.job_id
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Process images to base64
    img_tags = {}
    for k in ['interaction_plot', 'decomposition_plot', 'alanine_scanning_plot']:
        v = outputs.get(k)
        if v:
            p = v if os.path.isabs(v) else os.path.join(outputs_dir, v)
            b64 = _image_to_base64(p)
            if b64:
                img_tags[k] = f'<img src="data:image/png;base64,{b64}" alt="{k}">'
            else:
                img_tags[k] = '<div class="no-data">Image not found</div>'
        else:
            img_tags[k] = ''

    # Extract parameters
    sdict = payload.settings.model_dump() if hasattr(payload.settings, "model_dump") else payload.settings.dict()
    model = (sdict.get("implicit_model") or "GB")
    gb = (sdict.get("gb_adv") or {})
    pb = (sdict.get("pb_adv") or {})
    
    # Read summary CSV for table
    binding_rows = []
    sum_csv = os.path.join(outputs_dir, "interaction_PL_summary.csv")
    if os.path.exists(sum_csv):
        try:
            with open(sum_csv, 'r') as f:
                reader = csv.DictReader(f)
                binding_rows = list(reader)
        except: pass

    # Read Ala Scan CSV
    ala_rows = []
    ala_csv = os.path.join(outputs_dir, "alanine_scanning_results.csv")
    if os.path.exists(ala_csv):
        try:
            with open(ala_csv, 'r') as f:
                reader = csv.DictReader(f)
                ala_rows = list(reader)
        except: pass

    # Read Decomposition CSV
    decomp_rows = []
    decomp_csv = os.path.join(outputs_dir, "residues_energy_summary.csv")
    if os.path.exists(decomp_csv):
        try:
            with open(decomp_csv, 'r') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    try:
                        # Convert to float for sorting
                        row['TOTAL_FLOAT'] = float(row.get('TOTAL', 0))
                        decomp_rows.append(row)
                    except: pass
            # Sort by total energy (ascending, lower is better binding)
            decomp_rows.sort(key=lambda x: x['TOTAL_FLOAT'])
            decomp_rows = decomp_rows[:10]
        except: pass

    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MMPBSA 分析报告 - {job_id}</title>
    <script src="https://polyfill.io/v3/polyfill.min.js?features=es6"></script>
    <script id="MathJax-script" async src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"></script>
    <script>
        window.MathJax = {{
            tex: {{
                inlineMath: [['$', '$'], ['\\\\(', '\\\\)']]
            }}
        }};
    </script>
    <style>
        :root {{ --primary: #2563eb; --text: #1e293b; --bg: #f8fafc; --card: #ffffff; --border: #e2e8f0; }}
        body {{ font-family: 'Inter', system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 40px 20px; line-height: 1.6; }}
        .container {{ max-width: 900px; margin: 0 auto; background: var(--card); padding: 40px; border-radius: 16px; box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1); }}
        header {{ margin-bottom: 40px; border-bottom: 2px solid var(--border); padding-bottom: 20px; }}
        h1 {{ margin: 0 0 10px 0; font-size: 28px; color: var(--primary); font-weight: 700; }}
        .meta {{ color: #64748b; font-size: 14px; }}
        h2 {{ margin: 40px 0 20px 0; font-size: 20px; color: #334155; border-left: 4px solid var(--primary); padding-left: 12px; font-weight: 600; }}
        p {{ margin-bottom: 16px; color: #475569; text-align: justify; }}
        ul {{ margin-bottom: 16px; padding-left: 20px; color: #475569; }}
        table {{ width: 100%; border-collapse: collapse; font-size: 14px; margin: 20px 0; border: 1px solid var(--border); border-radius: 8px; overflow: hidden; }}
        th, td {{ padding: 12px 16px; text-align: left; border-bottom: 1px solid var(--border); }}
        th {{ background: #f8fafc; font-weight: 600; color: #334155; }}
        .chart-card {{ border: 1px solid var(--border); border-radius: 12px; padding: 16px; background: #fff; margin-bottom: 24px; text-align: center; }}
        .chart-card img {{ max-width: 100%; height: auto; border-radius: 8px; }}
        .citation {{ font-size: 12px; color: #64748b; background: #f1f5f9; padding: 10px; border-radius: 6px; margin-top: 10px; }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>MMPBSA 结合自由能分析报告</h1>
            <div class="meta">任务 ID: {job_id} | 生成时间: {now}</div>
        </header>

        <div class="section">
            <h2>1. 方法学 (Methodology)</h2>
            <p>结合自由能 ($\Delta G_{{bind}}$) 使用 {model} 隐式溶剂模型计算。总结合自由能估算公式如下：</p>
            <p style="text-align:center; font-weight:bold;">$$\Delta G_{{bind}} = \Delta H - T\Delta S \\approx \Delta E_{{MM}} + \Delta G_{{sol}} - T\Delta S$$</p>
            <p>其中：</p>
            <ul>
                <li>$\Delta E_{{MM}}$: 气相相互作用能 (范德华力 + 静电)</li>
                <li>$\Delta G_{{sol}}$: 溶剂化自由能 (极性 + 非极性)</li>
                <li>$T\Delta S$: 熵贡献 (如已计算)</li>
            </ul>
        </div>

        <div class="section">
            <h2>2. 模拟参数 (Simulation Parameters)</h2>
            <table>
                <tr><th>参数</th><th>值</th></tr>
                <tr><td>隐式模型</td><td>{model}</td></tr>
                <tr><td>GB 模型 (IGB)</td><td>{gb.get('igb', 'N/A') if model=='GB' else 'N/A'}</td></tr>
                <tr><td>盐浓度</td><td>{gb.get('salt', '0.15') if model=='GB' else pb.get('ionic_s', '0.15')} M</td></tr>
                <tr><td>采样帧</td><td>1 - 10 (间隔: 1)</td></tr>
            </table>
        </div>

        <div class="section">
            <h2>3. 结果：结合自由能 (Binding Free Energy)</h2>
            <p>下表汇总了构成结合自由能的各项能量分量。</p>
            
            <table>
                <thead>
                    <tr>
                        <th>能量分量</th>
                        <th>均值 (kcal/mol)</th>
                        <th>标准差</th>
                    </tr>
                </thead>
                <tbody>
    """
    
    if binding_rows:
        for r in binding_rows:
            html += f"<tr><td>{r.get('component','')}</td><td>{r.get('mean','')}</td><td>{r.get('error','')}</td></tr>"
    else:
        html += "<tr><td colspan='3' style='text-align:center'>无可用数据</td></tr>"

    html += f"""
                </tbody>
            </table>
            
            <div class="chart-card">
                {img_tags.get('interaction_plot', '')}
                <div style="margin-top:8px; font-size:12px; color:#666">图 1. 结合能分量随模拟时间的变化</div>
            </div>
        </div>
    """

    if img_tags.get('decomposition_plot') or decomp_rows:
        html += f"""
        <div class="section">
            <h2>4. 残基能量分解 (Residue Decomposition)</h2>
            <p>通过残基自由能分解识别对结合有关键贡献的残基（能量越低结合越强）。</p>
        """
        
        if decomp_rows:
            html += """
            <table>
                <thead>
                    <tr>
                        <th>残基</th>
                        <th>总能量 (kcal/mol)</th>
                        <th>范德华贡献</th>
                        <th>静电贡献</th>
                    </tr>
                </thead>
                <tbody>
            """
            for r in decomp_rows:
                html += f"<tr><td>{r.get('Residue','')}</td><td>{r.get('TOTAL','')}</td><td>{r.get('VDW','')}</td><td>{r.get('EEL','')}</td></tr>"
            html += """
                </tbody>
            </table>
            """

        html += f"""
            <div class="chart-card">
                {img_tags.get('decomposition_plot', '')}
                <div style="margin-top:8px; font-size:12px; color:#666">图 2. 结合自由能贡献最大的残基</div>
            </div>
        </div>
        """

    if ala_rows:
        html += f"""
        <div class="section">
            <h2>5. 丙氨酸扫描结果 (Alanine Scanning)</h2>
            <p>计算丙氨酸扫描用于识别热点残基。正的 $\Delta\Delta G$ 值表示突变为丙氨酸后复合物稳定性降低（结合亲和力损失），即该残基对结合很重要。</p>
            
            <table>
                <thead>
                    <tr>
                        <th>残基</th>
                        <th>$\Delta\Delta G$ (kcal/mol)</th>
                        <th>误差</th>
                    </tr>
                </thead>
                <tbody>
        """
        for r in ala_rows:
            html += f"<tr><td>{r.get('Residue','')}</td><td>{r.get('ddG','')}</td><td>{r.get('Error','')}</td></tr>"
            
        html += f"""
                </tbody>
            </table>
            
            <div class="chart-card">
                {img_tags.get('alanine_scanning_plot', '')}
                <div style="margin-top:8px; font-size:12px; color:#666">图 3. 丙氨酸扫描 $\Delta\Delta G$ 值</div>
            </div>
        </div>
        """

    html += """
        <div class="section">
            <h2>6. 参考文献 (References)</h2>
            <div class="citation">
                1. Valdés-Tresanco, M.S., et al. "gmx_MMPBSA: A New Tool to Perform End-State Free Energy Calculations with GROMACS." J. Chem. Theory Comput. 2021, 17, 10, 6281–6291.<br>
                2. Miller, B.R., et al. "MMPBSA.py: An Efficient Program for End-State Free Energy Calculations." J. Chem. Theory Comput. 2012, 8, 9, 3314–3321.
            </div>
        </div>
    </div>
</body>
</html>
    """
    
    report_path = os.path.join(outputs_dir, "report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    return "report.html"

def _generate_markdown_report(outputs_dir: str, payload: MmpRunInput) -> Optional[str]:
    lines = []
    lines.append("# MMPBSA Analysis Report")
    lines.append(f"**Job ID:** {payload.job_id}")
    lines.append(f"**Model:** {payload.settings.implicit_model}")
    
    # 1. Binding Energy
    lines.append("\n## 1. Binding Free Energy Components")
    sum_csv = os.path.join(outputs_dir, "interaction_PL_summary.csv")
    if os.path.exists(sum_csv):
        lines.append("| Component | Mean | Std. Dev. |")
        lines.append("| --- | --- | --- |")
        try:
            with open(sum_csv, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    comp = row.get("component", "")
                    mean = row.get("mean", "")
                    err = row.get("error", "")
                    lines.append(f"| {comp} | {mean} | {err} |")
        except:
            lines.append("Error reading summary.")
    else:
        lines.append("No binding energy summary found.")

    # 2. Key Residues
    lines.append("\n## 2. Key Interaction Residues")
    decomp_csv = os.path.join(outputs_dir, "residues_energy_summary.csv")
    if os.path.exists(decomp_csv):
        lines.append("Top 10 residues contributing to binding energy (lower is stronger binding):")
        lines.append("| Residue | Total Energy | VdW | Elec |")
        lines.append("| --- | --- | --- | --- |")
        try:
            data = []
            with open(decomp_csv, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    try:
                        tot = float(row.get("TOTAL", 0))
                        data.append((row, tot))
                    except:
                        pass
            # Sort by total energy (ascending)
            data.sort(key=lambda x: x[1])
            for item in data[:10]:
                r = item[0]
                lines.append(f"| {r.get('Residue')} | {r.get('TOTAL')} | {r.get('VDW')} | {r.get('EEL')} |")
        except:
            lines.append("Error reading decomposition data.")
    else:
        lines.append("No decomposition data found.")

    # 3. Alanine Scanning
    ala_csv = os.path.join(outputs_dir, "alanine_scanning_results.csv")
    if os.path.exists(ala_csv):
        lines.append("\n## 3. Alanine Scanning Results")
        lines.append("Positive ddG indicates that mutation to Alanine destabilizes binding (residue is important).")
        lines.append("| Residue | ddG (kcal/mol) | Error |")
        lines.append("| --- | --- | --- |")
        try:
            with open(ala_csv, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                for row in reader:
                    lines.append(f"| {row.get('Residue')} | {row.get('ddG')} | {row.get('Error')} |")
        except:
            lines.append("Error reading alanine scanning data.")
            
    report_path = os.path.join(outputs_dir, "report.md")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    return "report.md"

@app.post("/run")
def run(payload: MmpRunInput):
    job_id = payload.job_id
    job_dir, inputs_dir, outputs_dir = ensure_job_dirs(job_id, payload.uid)
    with open(os.path.join(inputs_dir, "inputs.json"), "w", encoding="utf-8") as f:
        json.dump(payload.model_dump(), f, ensure_ascii=False)
    log_path = os.path.join(outputs_dir, "log.txt")
    write_text(log_path, "mmpbsa_start")
    src_work = os.path.join(DATA_DIR, payload.uid, "md", payload.source_job_id, "outputs", "work")
    outputs = {}
    tpr = os.path.join(src_work, "md.tpr")
    xtc = os.path.join(src_work, "md.xtc")
    gro = os.path.join(src_work, "processed.gro")
    top = os.path.join(src_work, "topol.top")
    base = find_gmx_mmpbsa_cmd()
    
    if not (base and os.path.exists(tpr) and os.path.exists(xtc) and os.path.exists(top)):
        append_log(log_path, "inputs_or_binary_missing")
        append_log(log_path, "mmpbsa_complete")
        return {"job_id": job_id, "outputs": outputs}

    ndx, rec_idx, lig_idx = create_index_and_groups(gro, outputs_dir)
    if not (ndx and rec_idx is not None and lig_idx is not None):
        append_log(log_path, "index_or_groups_missing")
        append_log(log_path, "mmpbsa_complete")
        return {"job_id": job_id, "outputs": outputs}

    sdict = payload.settings.model_dump() if hasattr(payload.settings, "model_dump") else payload.settings.dict()
    gb = (sdict.get("gb_adv") or {})
    startframe = 1
    endframe = 10
    model = (sdict.get("implicit_model") or "GB")
    igb = int(gb.get("igb") or 5)
    salt = float(gb.get("salt") or 0.150)
    pb_adv = sdict.get("pb_adv")
    rism_adv = sdict.get("rism_adv")
    
    ala_settings = sdict.get("alanine") or {}
    scan_range = ala_settings.get("range") or sdict.get("scan_range") or "manual"
    ala_res_input = ala_settings.get("residues") or sdict.get("ala_res")
    interface_cutoff = float(ala_settings.get("interface_cutoff") or sdict.get("interface_cutoff") or 10.0)
    probe_radius = float(ala_settings.get("probe_radius") or 1.4)
    
    # Check if alanine scanning is actually requested
    # 1. Check functions flag (from frontend checkbox)
    funcs = sdict.get("functions") or {}
    is_alanine_requested = funcs.get("alanine", False)
    
    # 2. Fallback: if functions missing, check payload mode (legacy)
    if "functions" not in sdict and payload.mode == "alanine":
        is_alanine_requested = True
        
    residues_to_scan = []
    if is_alanine_requested:
        if scan_range == "interface":
            append_log(log_path, f"Detecting interface residues (cutoff={interface_cutoff}A)...")
            found_residues = _get_interface_residues(tpr, xtc, outputs_dir, interface_cutoff, log_path=log_path, ndx_path=ndx)
            append_log(log_path, f"Found {len(found_residues)} residues: {', '.join(found_residues)}")
            # Always perform standard analysis (None) first, then scan residues
            residues_to_scan = [None] + found_residues
        elif ala_res_input:
            residues_to_scan = [ala_res_input]

    if not residues_to_scan:
        residues_to_scan = [None] 

    import subprocess
    import shutil

    results_data = []
    
    # If standard run (res is None), only 1 iteration
    for i, res in enumerate(residues_to_scan):
        # Create separate directory for each run to avoid file conflicts
        if res:
            safe_res = res.replace(":", "_")
            run_subdir = os.path.join(outputs_dir, f"ala_{safe_res}")
        else:
            run_subdir = os.path.join(outputs_dir, "reference_calc")
        
        os.makedirs(run_subdir, exist_ok=True)
        
        mmp_in = os.path.join(run_subdir, "mmpbsa.in")
        write_mmpbsa_in(mmp_in, startframe=startframe, endframe=endframe, model=model, igb=igb, saltcon=salt, pb_opts=ImplicitPbAdv(**pb_adv) if pb_adv else None, rism_opts=RismAdv(**rism_adv) if rism_adv else None, include_decomp=True, ala_res=res)
        
        # Note: Input files are absolute paths, so they work from subdir.
        # But 'mmp_in' is now in run_subdir.
        cmd = list(base) + ["-O", "-nogui", "-i", mmp_in, "-cs", tpr, "-ct", xtc, "-ci", ndx, "-cg", str(rec_idx), str(lig_idx), "-cp", top, "-o", "FINAL_RESULTS_MMPBSA.dat", "-eo", "FINAL_RESULTS_MMPBSA.csv"]
        
        if res:
            msg = f"Step {i+1}/{len(residues_to_scan)}: Running Alanine Scanning for residue {res}..."
        else:
            msg = f"Step {i+1}/{len(residues_to_scan)}: Running Standard MMPBSA Analysis (Reference)..."
        
        append_log(log_path, msg)
        
        # Run in subdirectory
        subprocess.run(cmd, cwd=run_subdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        
        dat_path = os.path.join(run_subdir, "FINAL_RESULTS_MMPBSA.dat")
        csv_path = os.path.join(run_subdir, "FINAL_RESULTS_MMPBSA.csv")
        
        if res:
            if os.path.exists(dat_path):
                # Copy result CSV to main outputs dir for easier access/download
                safe_res = res.replace(":", "_")
                if os.path.exists(csv_path):
                    shutil.copy(csv_path, os.path.join(outputs_dir, f"FINAL_RESULTS_{safe_res}.csv"))
                
                ddg_val, ddg_err = _extract_ddg(dat_path)
                results_data.append({"Residue": res, "ddG": ddg_val, "Error": ddg_err})
        else:
            if os.path.exists(csv_path):
                append_log(log_path, f"find {csv_path}")
                # Copy main results to outputs_dir
                shutil.copy(csv_path, os.path.join(outputs_dir, "FINAL_RESULTS_MMPBSA.csv"))
                outputs["results_csv"] = "FINAL_RESULTS_MMPBSA.csv"
                
                inter2 = _generate_interaction_from_final_csv(os.path.join(outputs_dir, "FINAL_RESULTS_MMPBSA.csv"), outputs_dir)
                append_log(log_path, f"inter_from_csv: {inter2}")
                outputs.update(inter2)
            
            decomp_dat = os.path.join(run_subdir, "FINAL_DECOMP_MMPBSA.dat")
            if os.path.exists(decomp_dat):
                 shutil.copy(decomp_dat, os.path.join(outputs_dir, "FINAL_DECOMP_MMPBSA.dat"))
                 decomp_out = _parse_final_decomp_dat(os.path.join(outputs_dir, "FINAL_DECOMP_MMPBSA.dat"), outputs_dir, log_path)
                 append_log(log_path, f"decomp_out: {decomp_out}")
                 outputs.update(decomp_out)

    if results_data:
        sum_csv = os.path.join(outputs_dir, "alanine_scanning_results.csv")
        with open(sum_csv, 'w') as f:
            f.write("Residue,ddG,Error\n")
            for r in results_data:
                f.write(f"{r['Residue']},{r['ddG']},{r['Error']}\n")
        outputs["alanine_scanning_csv"] = "alanine_scanning_results.csv"
        
        _plot_alanine_scanning(results_data, outputs_dir)
        outputs["alanine_scanning_plot"] = "alanine_scanning.png"

    try:
        # Generate both Markdown (legacy) and HTML reports
        _generate_markdown_report(outputs_dir, payload)
        report_file = _generate_html_report(outputs_dir, payload, outputs)
        if report_file:
            outputs["report"] = report_file
    except Exception as e:
        append_log(log_path, f"Report generation failed: {e}")

    append_log(log_path, "mmpbsa_complete")
    return {"job_id": job_id, "outputs": outputs}
