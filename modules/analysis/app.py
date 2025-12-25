import os
import json
import re
import subprocess
import datetime
import base64
from typing import Optional
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

try:
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except Exception:
    np = None

DATA_DIR = os.getenv("ENZYME_DATA_DIR", "/data")

class GeneralOpts(BaseModel):
    temp: bool = False
    press: bool = False
    pot: bool = False
    kin: bool = False

class StructOpts(BaseModel):
    group: str = "System"
    group_index: Optional[str] = None
    rmsd: bool = False
    rg: bool = False
    rmsf: bool = False

class MmpbsaOpts(BaseModel):
    enable: bool = False

class AnaRunInput(BaseModel):
    job_id: str
    uid: str
    source_job_id: str
    general: GeneralOpts
    struct: StructOpts
    mmpbsa: MmpbsaOpts

app = FastAPI()

def ensure_job_dirs(job_id: str, uid: str):
    job_dir = os.path.join(DATA_DIR, uid, "analysis", job_id)
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

def load_xy(path: str) -> Optional[list]:
    try:
        xs, ys = [], []
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            for ln in f:
                ln = ln.strip()
                if not ln or ln.startswith(('#','@')):
                    continue
                parts = ln.split()
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
        
        # Add a light background for contrast if needed, but white is standard for publication
        # fig.patch.set_facecolor('white')
        
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
        s = (sum((v-mu)**2 for v in vals)/len(vals)) ** 0.5
    return (m, s)

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
        append_log(log_path, f"run: gmx energy -f {edr_path} -o {out_path} -xvg none < {idx}")
        p = subprocess.run(["gmx", "energy", "-f", edr_path, "-o", out_path, "-xvg", "none"], cwd=work_dir, input=inp, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        append_log(log_path, p.stdout or "")
        if p.returncode != 0:
            return None
        return os.path.join(work_dir, out_path)
    except Exception:
        return None

def run_gmx_rms(tpr_path: str, xtc_path: str, out_path: str, work_dir: str, log_path: str, sel_idx: int) -> Optional[str]:
    try:
        inp = f"{sel_idx}\n{sel_idx}\n"
        append_log(log_path, f"run: gmx rms -s {tpr_path} -f {xtc_path} -o {out_path} -xvg none < {sel_idx}")
        p = subprocess.run(["gmx", "rms", "-s", tpr_path, "-f", xtc_path, "-o", out_path, "-xvg", "none"], cwd=work_dir, input=inp, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        append_log(log_path, p.stdout or "")
        if p.returncode != 0:
            return None
        return os.path.join(work_dir, out_path)
    except Exception:
        return None

def run_gmx_gyrate(tpr_path: str, xtc_path: str, out_path: str, work_dir: str, log_path: str, sel_idx: int) -> Optional[str]:
    try:
        inp = f"{sel_idx}\n"
        append_log(log_path, f"run: gmx gyrate -s {tpr_path} -f {xtc_path} -o {out_path} -xvg none < {sel_idx}")
        p = subprocess.run(["gmx", "gyrate", "-s", tpr_path, "-f", xtc_path, "-o", out_path, "-xvg", "none"], cwd=work_dir, input=inp, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        append_log(log_path, p.stdout or "")
        if p.returncode != 0:
            return None
        return os.path.join(work_dir, out_path)
    except Exception:
        return None

def run_gmx_rmsf(tpr_path: str, xtc_path: str, out_path: str, work_dir: str, log_path: str, sel_idx: int) -> Optional[str]:
    try:
        inp = f"{sel_idx}\n"
        append_log(log_path, f"run: gmx rmsf -s {tpr_path} -f {xtc_path} -o {out_path} -xvg none < {sel_idx}")
        p = subprocess.run(["gmx", "rmsf", "-s", tpr_path, "-f", xtc_path, "-o", out_path, "-xvg", "none"], cwd=work_dir, input=inp, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        append_log(log_path, p.stdout or "")
        if p.returncode != 0:
            return None
        return os.path.join(work_dir, out_path)
    except Exception:
        return None

def get_md_params(uid, source_job_id):
    try:
        path = os.path.join(DATA_DIR, uid, "md", source_job_id, "inputs", "inputs.json")
        if os.path.exists(path):
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
    except Exception:
        pass
    return {}

def image_to_base64(path):
    try:
        if not os.path.exists(path):
            return ""
        with open(path, "rb") as f:
            return base64.b64encode(f.read()).decode('utf-8')
    except Exception:
        return ""

def generate_html_report(job_id, outputs_dir, gen_rows, struct_rows, outputs, md_params):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    # Process images to base64
    img_tags = {}
    for k, v in outputs.items():
        if k.endswith("_plot") or k.endswith(".png"):
            # If v is absolute path, use it; else join with outputs_dir
            p = v if os.path.isabs(v) else os.path.join(outputs_dir, v)
            b64 = image_to_base64(p)
            if b64:
                img_tags[k] = f'<img src="data:image/png;base64,{b64}" alt="{k}">'
            else:
                img_tags[k] = '<div class="no-data">图片加载失败</div>'

    # Extract MD params
    steps = md_params.get("steps", {})
    ff = md_params.get("forcefield", "amber99sb")
    water = md_params.get("water", "tip3p")
    prod_ns = steps.get('prod_steps', 1000) * 0.002 / 1000
    
    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>模拟分析报告 - {job_id}</title>
    <style>
        :root {{ --primary: #4f46e5; --text: #1e293b; --bg: #f8fafc; --card: #ffffff; --border: #e2e8f0; }}
        body {{ font-family: 'Inter', system-ui, -apple-system, sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 40px 20px; line-height: 1.6; }}
        .container {{ max-width: 900px; margin: 0 auto; background: var(--card); padding: 40px; border-radius: 16px; box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06); }}
        header {{ margin-bottom: 40px; border-bottom: 2px solid var(--border); padding-bottom: 20px; }}
        h1 {{ margin: 0 0 10px 0; font-size: 28px; color: var(--primary); font-weight: 700; }}
        .meta {{ color: #64748b; font-size: 14px; }}
        h2 {{ margin: 40px 0 20px 0; font-size: 20px; color: #334155; border-left: 4px solid var(--primary); padding-left: 12px; font-weight: 600; }}
        h3 {{ margin: 24px 0 12px 0; font-size: 16px; color: #475569; font-weight: 600; }}
        p {{ margin-bottom: 16px; color: #475569; text-align: justify; }}
        ul {{ margin-bottom: 16px; padding-left: 20px; color: #475569; }}
        li {{ margin-bottom: 8px; }}
        .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(350px, 1fr)); gap: 32px; margin-top: 24px; }}
        .chart-card {{ border: 1px solid var(--border); border-radius: 12px; padding: 16px; background: #fff; box-shadow: 0 1px 2px rgba(0,0,0,0.05); }}
        .chart-card img {{ width: 100%; height: auto; display: block; border-radius: 8px; }}
        .chart-title {{ font-size: 14px; font-weight: 600; margin-bottom: 12px; color: #475569; text-align: center; }}
        table {{ width: 100%; border-collapse: collapse; font-size: 14px; margin: 20px 0; border: 1px solid var(--border); border-radius: 8px; overflow: hidden; }}
        th, td {{ padding: 12px 16px; text-align: left; border-bottom: 1px solid var(--border); }}
        th {{ background: #f8fafc; font-weight: 600; color: #334155; }}
        tr:last-child td {{ border-bottom: none; }}
        .no-data {{ color: #94a3b8; font-style: italic; padding: 20px; text-align: center; background: #f8fafc; border-radius: 8px; }}
        .footer {{ margin-top: 60px; pt: 20px; border-top: 1px solid var(--border); font-size: 12px; color: #94a3b8; text-align: center; }}
        strong {{ color: #334155; font-weight: 600; }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>分子动力学模拟分析报告</h1>
            <div class="meta">任务 ID: {job_id} | 生成时间: {now}</div>
        </header>

        <div class="section">
            <h2>1. 方法 (Methods)</h2>
            <p>分子动力学模拟使用 <strong>GROMACS</strong> 软件进行。蛋白质力场选用 <strong>{ff}</strong>，水模型选用 <strong>{water}</strong>。</p>
            <p>模拟流程包括以下步骤：</p>
            <ul>
                <li><strong>能量最小化 (Energy Minimization)</strong>：使用最速下降法 (Steepest Descent)，最大步数 {steps.get('em_steps', 5000)}，收敛容差 {steps.get('em_tol', 1000.0)} kJ/mol/nm，以消除初始结构的不合理接触。</li>
                <li><strong>NPT 平衡</strong>：在恒温恒压系综下进行 {steps.get('npt_steps', 1000) * 0.002:.1f} ps 的平衡模拟。
                    <ul>
                        <li>控温方法：V-rescale，参考温度 {steps.get('npt_temp', 300.0)} K，耦合时间常数 0.1 ps。</li>
                        <li>控压方法：Berendsen，参考压强 1.0 bar，耦合时间常数 2.0 ps。</li>
                    </ul>
                </li>
                <li><strong>生产动力学 (Production MD)</strong>：时长 {prod_ns:.3f} ns ({steps.get('prod_steps', 1000)} 步，步长 2 fs)。
                    <ul>
                        <li>控温方法：V-rescale，参考温度 {steps.get('prod_temp', 300.0)} K。</li>
                        <li>控压方法：{'Parrinello-Rahman' if steps.get('prod_press', 1.0) > 0 else '无 (NVT)'}，参考压强 {steps.get('prod_press', 1.0)} bar。</li>
                        <li>长程静电相互作用：PME (Particle Mesh Ewald)。</li>
                        <li>范德华相互作用：Cut-off (1.0 nm)。</li>
                        <li>约束：所有键均使用 LINCS 算法进行约束。</li>
                    </ul>
                </li>
            </ul>
        </div>

        <div class="section">
            <h2>2. 结果与讨论 (Results)</h2>
            
            <h3>2.1 通用性质 (General Properties)</h3>
            <p>模拟过程中的热力学性质稳定性如下表所示：</p>
            {'<table><thead><tr><th>性质</th><th>平均值</th><th>标准差</th></tr></thead><tbody>' + ''.join([f'<tr><td>{r[0]}</td><td>{r[1]}</td><td>{r[2]}</td></tr>' for r in gen_rows]) + '</tbody></table>' if gen_rows else '<div class="no-data">无数据</div>'}
            
            <div class="grid">
                {''.join([f'<div class="chart-card"><div class="chart-title">温度 (Temperature)</div>{img_tags.get("gen_temp_plot", "")}</div>' if "gen_temp_plot" in outputs else ''])}
                {''.join([f'<div class="chart-card"><div class="chart-title">压强 (Pressure)</div>{img_tags.get("gen_press_plot", "")}</div>' if "gen_press_plot" in outputs else ''])}
                {''.join([f'<div class="chart-card"><div class="chart-title">势能 (Potential Energy)</div>{img_tags.get("gen_pot_plot", "")}</div>' if "gen_pot_plot" in outputs else ''])}
                {''.join([f'<div class="chart-card"><div class="chart-title">动能 (Kinetic Energy)</div>{img_tags.get("gen_kin_plot", "")}</div>' if "gen_kin_plot" in outputs else ''])}
            </div>

            <h3>2.2 结构性质 (Structural Properties)</h3>
            <p>对模拟轨迹的结构分析结果如下：</p>
            {'<table><thead><tr><th>性质</th><th>组 (Group)</th><th>平均值</th><th>标准差</th></tr></thead><tbody>' + ''.join([f'<tr><td>{r[0]}</td><td>{r[1]}</td><td>{r[2]}</td><td>{r[3]}</td></tr>' for r in struct_rows]) + '</tbody></table>' if struct_rows else '<div class="no-data">无数据</div>'}
            
            <div class="grid">
                {''.join([f'<div class="chart-card"><div class="chart-title">RMSD (Root Mean Square Deviation)</div>{img_tags.get("struct_rmsd_plot", "")}</div>' if "struct_rmsd_plot" in outputs else ''])}
                {''.join([f'<div class="chart-card"><div class="chart-title">Radius of Gyration (Rg)</div>{img_tags.get("struct_rg_plot", "")}</div>' if "struct_rg_plot" in outputs else ''])}
                {''.join([f'<div class="chart-card"><div class="chart-title">RMSF (Root Mean Square Fluctuation)</div>{img_tags.get("struct_rmsf_plot", "")}</div>' if "struct_rmsf_plot" in outputs else ''])}
            </div>
        </div>

        <div class="section">
            <h2>3. 引用 (References)</h2>
            <ul>
                <li><strong>GROMACS:</strong> Abraham, M. J., Murtola, T., Schulz, R., Páll, S., Smith, J. C., Hess, B., & Lindahl, E. (2015). GROMACS: High performance molecular simulations through multi-level parallelism from laptops to supercomputers. <em>SoftwareX</em>, 1, 19-25.</li>
                <li><strong>V-rescale Thermostat:</strong> Bussi, G., Donadio, D., & Parrinello, M. (2007). Canonical sampling through velocity rescaling. <em>The Journal of chemical physics</em>, 126(1), 014101.</li>
                <li><strong>Parrinello-Rahman Barostat:</strong> Parrinello, M., & Rahman, A. (1981). Polymorphic transitions in single crystals: A new molecular dynamics method. <em>Journal of Applied physics</em>, 52(12), 7182-7190.</li>
                <li><strong>PME:</strong> Darden, T., York, D., & Pedersen, L. (1993). Particle mesh Ewald: An N⋅log(N) method for Ewald sums in large systems. <em>The Journal of chemical physics</em>, 98(12), 10089-10092.</li>
                <li><strong>LINCS:</strong> Hess, B., Bekker, H., Berendsen, H. J., & Fraaije, J. G. (1997). LINCS: a linear constraint solver for molecular simulations. <em>Journal of computational chemistry</em>, 18(12), 1463-1472.</li>
            </ul>
        </div>

        <div class="footer">
            Generated by Enzyme Platform Analysis Module | {now}
        </div>
    </div>
</body>
</html>
    """
    
    report_path = os.path.join(outputs_dir, "analysis_report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html)
    return "analysis_report.html"

@app.post("/run")
def run(payload: AnaRunInput):
    job_id = payload.job_id
    job_dir, inputs_dir, outputs_dir = ensure_job_dirs(job_id, payload.uid)
    with open(os.path.join(inputs_dir, "inputs.json"), "w", encoding="utf-8") as f:
        json.dump(payload.model_dump(), f, ensure_ascii=False)
    log_path = os.path.join(outputs_dir, "log.txt")
    write_text(log_path, "analysis_start")
    src_work = os.path.join(DATA_DIR, payload.uid, "md", payload.source_job_id, "outputs", "work")
    outputs = {}
    rows = []
    srows = []
    # General properties using gmx energy from md.edr
    try:
        edr_path = os.path.join(src_work, "md.edr")
        if not os.path.exists(edr_path):
            append_log(log_path, "missing_md_edr")
        else:
            if payload.general.temp:
                t_idx = parse_energy_index(edr_path, "Temperature", src_work, log_path)
                if t_idx is not None:
                    xvgp = extract_energy_series(edr_path, t_idx, "ana_temp.xvg", src_work, log_path)
                    data = load_xy(xvgp) if xvgp else None
                    if data:
                        ok = plot_xy(data, "Temperature vs Time", "Time (ps)", "Temperature (K)", os.path.join(outputs_dir, "gen_temp.png"))
                        if ok and os.path.exists(os.path.join(outputs_dir, "gen_temp.png")):
                            outputs["gen_temp_plot"] = "gen_temp.png"
                        vals = [y for _, y in data]
                        m, s = mean_std(vals)
                        rows.append(["Temperature", ("%0.2f" % m) if m is not None else "-", ("%0.2f" % s) if s is not None else "-"])
            if payload.general.press:
                p_idx = parse_energy_index(edr_path, "Pressure", src_work, log_path)
                if p_idx is not None:
                    xvgp = extract_energy_series(edr_path, p_idx, "ana_pressure.xvg", src_work, log_path)
                    data = load_xy(xvgp) if xvgp else None
                    if data:
                        ok = plot_xy(data, "Pressure vs Time", "Time (ps)", "Pressure (bar)", os.path.join(outputs_dir, "gen_press.png"))
                        if ok and os.path.exists(os.path.join(outputs_dir, "gen_press.png")):
                            outputs["gen_press_plot"] = "gen_press.png"
                        vals = [y for _, y in data]
                        m, s = mean_std(vals)
                        rows.append(["Pressure", ("%0.3f" % m) if m is not None else "-", ("%0.3f" % s) if s is not None else "-"])
            if payload.general.pot:
                pot_idx = parse_energy_index(edr_path, "Potential", src_work, log_path)
                if pot_idx is not None:
                    xvgp = extract_energy_series(edr_path, pot_idx, "ana_potential.xvg", src_work, log_path)
                    data = load_xy(xvgp) if xvgp else None
                    if data:
                        ok = plot_xy(data, "Potential vs Time", "Time (ps)", "Potential (kJ/mol)", os.path.join(outputs_dir, "gen_pot.png"))
                        if ok and os.path.exists(os.path.join(outputs_dir, "gen_pot.png")):
                            outputs["gen_pot_plot"] = "gen_pot.png"
                        vals = [y for _, y in data]
                        m, s = mean_std(vals)
                        rows.append(["Potential", ("%0.1f" % m) if m is not None else "-", ("%0.1f" % s) if s is not None else "-"])
            if payload.general.kin:
                kin_idx = parse_energy_index(edr_path, "Kinetic", src_work, log_path)
                if kin_idx is not None:
                    xvgp = extract_energy_series(edr_path, kin_idx, "ana_kinetic.xvg", src_work, log_path)
                    data = load_xy(xvgp) if xvgp else None
                    if data:
                        ok = plot_xy(data, "Kinetic vs Time", "Time (ps)", "Kinetic (kJ/mol)", os.path.join(outputs_dir, "gen_kin.png"))
                        if ok and os.path.exists(os.path.join(outputs_dir, "gen_kin.png")):
                            outputs["gen_kin_plot"] = "gen_kin.png"
                        vals = [y for _, y in data]
                        m, s = mean_std(vals)
                        rows.append(["Kinetic", ("%0.1f" % m) if m is not None else "-", ("%0.1f" % s) if s is not None else "-"])
        gen_csv = os.path.join(outputs_dir, "general_table.csv")
        if rows:
            write_text(gen_csv, "property,mean,error\n" + "\n".join([",".join(r) for r in rows]))
            outputs["general_table"] = "general_table.csv"
        append_log(log_path, "general_done")
    except Exception:
        append_log(log_path, "general_error")
    # Structural properties from md.xtc/md.tpr
    try:
        tpr = os.path.join(src_work, "md.tpr")
        xtc = os.path.join(src_work, "md.xtc")
        if not os.path.exists(tpr) or not os.path.exists(xtc):
            append_log(log_path, "missing_md_tpr_or_xtc")
        else:
            sel_idx = 0
            try:
                if payload.struct.group and payload.struct.group.lower() == "custom" and payload.struct.group_index:
                    sel_idx = int(payload.struct.group_index)
            except Exception:
                sel_idx = 0
            if payload.struct.rmsd:
                xvgp = run_gmx_rms(tpr, xtc, "ana_rmsd.xvg", outputs_dir, log_path, sel_idx)
                data = load_xy(xvgp) if xvgp else None
                if data:
                    ok = plot_xy(data, "RMSD vs Time", "Time (ps)", "RMSD (nm)", os.path.join(outputs_dir, "struct_rmsd.png"))
                    if ok and os.path.exists(os.path.join(outputs_dir, "struct_rmsd.png")):
                        outputs["struct_rmsd_plot"] = "struct_rmsd.png"
                    vals = [y for _, y in data]
                    m, s = mean_std(vals)
                    srows.append(["RMSD", payload.struct.group or "System", ("%0.3f" % m) if m is not None else "-", ("%0.3f" % s) if s is not None else "-"])
                    append_log(log_path, f"rmsd_points={len(vals)}")
                else:
                    append_log(log_path, "rmsd_no_data")
            if payload.struct.rg:
                xvgp = run_gmx_gyrate(tpr, xtc, "ana_gyrate.xvg", outputs_dir, log_path, sel_idx)
                data = load_xy(xvgp) if xvgp else None
                if data:
                    ok = plot_xy(data, "Radius of Gyration vs Time", "Time (ps)", "Rg (nm)", os.path.join(outputs_dir, "struct_rg.png"))
                    if ok and os.path.exists(os.path.join(outputs_dir, "struct_rg.png")):
                        outputs["struct_rg_plot"] = "struct_rg.png"
                    vals = [y for _, y in data]
                    m, s = mean_std(vals)
                    srows.append(["RG", payload.struct.group or "System", ("%0.3f" % m) if m is not None else "-", ("%0.3f" % s) if s is not None else "-"])
                    append_log(log_path, f"rg_points={len(vals)}")
                else:
                    append_log(log_path, "rg_no_data")
            if payload.struct.rmsf:
                xvgp = run_gmx_rmsf(tpr, xtc, "ana_rmsf.xvg", outputs_dir, log_path, sel_idx)
                data = load_xy(xvgp) if xvgp else None
                if data:
                    ok = plot_xy(data, "RMSF", "Index", "RMSF (nm)", os.path.join(outputs_dir, "struct_rmsf.png"))
                    if ok and os.path.exists(os.path.join(outputs_dir, "struct_rmsf.png")):
                        outputs["struct_rmsf_plot"] = "struct_rmsf.png"
                    vals = [y for _, y in data]
                    m, s = mean_std(vals)
                    srows.append(["RMSF", payload.struct.group or "System", ("%0.3f" % m) if m is not None else "-", ("%0.3f" % s) if s is not None else "-"])
                    append_log(log_path, f"rmsf_points={len(vals)}")
                else:
                    append_log(log_path, "rmsf_no_data")
        st_csv = os.path.join(outputs_dir, "struct_table.csv")
        if srows:
            write_text(st_csv, "property,group,mean,error\n" + "\n".join([",".join(r) for r in srows]))
            outputs["struct_table"] = "struct_table.csv"
        append_log(log_path, "struct_done")
    except Exception:
        append_log(log_path, "struct_error")
    # MMPBSA placeholder
    try:
        if payload.mmpbsa.enable:
            plt.figure(figsize=(3,2), dpi=120)
            plt.title("MMPBSA Total")
            plt.tight_layout()
            p1 = os.path.join(outputs_dir, "mmp_total.png")
            plt.savefig(p1)
            plt.close()
            outputs["mmp_total_plot"] = "mmp_total.png"
            plt.figure(figsize=(3,2), dpi=120)
            plt.title("MMPBSA Components")
            plt.tight_layout()
            p2 = os.path.join(outputs_dir, "mmp_comp.png")
            plt.savefig(p2)
            plt.close()
            outputs["mmp_comp_plot"] = "mmp_comp.png"
            mmp_csv = os.path.join(outputs_dir, "mmp_table.csv")
            write_text(mmp_csv, "component,value\n-")
            outputs["mmp_table"] = "mmp_table.csv"
            append_log(log_path, "mmpbsa_done")
    except Exception:
        append_log(log_path, "mmpbsa_error")
    # Report alias
    try:
        md_params = get_md_params(payload.uid, payload.source_job_id)
        rep_html = generate_html_report(job_id, outputs_dir, rows, srows, outputs, md_params)
        outputs["analysis_report"] = rep_html
        
        rep_csv = os.path.join(outputs_dir, "analysis_report.csv")
        write_text(rep_csv, "section,item\nGeneral,See general_table\nStruct,See struct_table\nMMPBSA,See mmp_table")
        outputs["report_csv"] = "analysis_report.csv"
    except Exception:
        pass
    append_log(log_path, "analysis_complete")
    return {"job_id": job_id, "outputs": outputs}
