"""Stage-resolved energy requirements. Version 2.1.0.

Run from any directory: python artifact/cryogenic_model.py
Only NumPy and Matplotlib are required. Inputs marked assumed are not measurements.
This artifact does not estimate an LLM speedup or rank qubit modalities.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
CONFIG_PATH = Path(__file__).with_name("config.json")


def finite(value, name, lower=0.0, upper=None, strict_lower=False):
    if not math.isfinite(value) or value < lower or (strict_lower and value == lower):
        raise ValueError(f"{name} is outside its allowed domain")
    if upper is not None and value > upper:
        raise ValueError(f"{name} is outside its allowed domain")
    return value


def cooling_factor(tc, eta=0.02, th=300.0):
    """Refrigeration work / lifted heat, not total electrical input."""
    finite(th, "hot temperature", strict_lower=True)
    finite(tc, "cold temperature", upper=th, strict_lower=True)
    finite(eta, "Carnot fraction", upper=1.0, strict_lower=True)
    return (th - tc) / (eta * tc)


def electrical_multiplier(tc, eta=0.02, th=300.0):
    """Total thermal amplification factor A (TAF), including direct input.

    Valid for locally supplied electrical energy fully dissipated at tc.
    """
    return 1.0 + cooling_factor(tc, eta, th)


def added_wall_energy(electrical_j, heat_by_stage, fixed_wall_w=0.0,
                      token_rate=1.0, th=300.0):
    """Per-token electrical input + incremental stage cooling + allocated fixed power.

    heat_by_stage contains (temperature_K, eta, incremental_heat_J_per_token).
    Electrical input includes all local supplies and warm supplies, once only.
    Fixed wall power must exclude the incremental costs already accounted for.
    Shared-plant stage factors must be marginal allocations, not repeated totals.
    """
    finite(electrical_j, "electrical energy")
    finite(fixed_wall_w, "fixed wall power")
    finite(token_rate, "token rate", strict_lower=True)
    cooling = 0.0
    for tc, eta, q in heat_by_stage:
        finite(q, "heat")
        cooling += cooling_factor(tc, eta, th) * q
    return electrical_j + cooling + fixed_wall_w / token_rate


def shot_count(epsilon, delta=0.05):
    """Sufficient Hoeffding budget for independent outcomes in [-1, 1]."""
    finite(epsilon, "statistical error", upper=2.0, strict_lower=True)
    finite(delta, "failure probability", upper=1.0, strict_lower=True)
    if delta == 1:
        raise ValueError("failure probability must be less than one")
    return math.ceil(2.0 * math.log(2.0 / delta) / epsilon**2)


def validate_config(c):
    for name in ("qubits", "tokens_per_decision"):
        if isinstance(c[name], bool) or not isinstance(c[name], int) or c[name] < 1:
            raise ValueError(f"{name} must be a positive integer")
    for name in ("controller_power_per_qubit_w", "shot_duration_s", "token_rate_per_s",
                 "available_controller_cooling_w_per_plant"):
        finite(c[name], name, strict_lower=True)
    finite(c["availability"], "availability", upper=1.0, strict_lower=True)
    finite(c["fixed_wall_power_w_per_plant"], "fixed wall power")
    finite(c["allowed_bias"], "allowed bias", upper=2.0)
    finite(c["illustrative_energy_fraction"], "energy fraction", upper=1.0, strict_lower=True)
    finite(c["illustrative_overhead_fraction"], "overhead fraction")
    if c["illustrative_overhead_fraction"] >= c["illustrative_energy_fraction"]:
        raise ValueError("the plotted overhead must be below displaced energy")
    shot_count(c["statistical_error"], c["failure_probability"])
    cooling_factor(c["controller_temperature_k"], c["carnot_fraction"], c["hot_temperature_k"])


def case_result(c, epsilon=None, rate=None, eta=None):
    validate_config(c)
    eps = c["statistical_error"] if epsilon is None else epsilon
    r = c["token_rate_per_s"] if rate is None else rate
    efficiency = c["carnot_fraction"] if eta is None else eta
    finite(r, "token rate", strict_lower=True)
    shots = shot_count(eps, c["failure_probability"])
    duration = shots * c["shot_duration_s"]
    power = c["qubits"] * c["controller_power_per_qubit_w"]
    a = electrical_multiplier(c["controller_temperature_k"], efficiency, c["hot_temperature_k"])
    e_control = a * power * duration
    nu = 1.0 / c["tokens_per_decision"]
    engines = math.ceil(r * nu * duration / c["availability"])
    per_plant = math.floor(c["available_controller_cooling_w_per_plant"] / power)
    if per_plant < 1:
        raise ValueError("one controller exceeds the available controller cooling budget")
    plants = math.ceil(engines / per_plant)
    fixed = plants * c["fixed_wall_power_w_per_plant"] / r
    return {
        "epsilon_stat": eps, "delta": c["failure_probability"], "shots": shots,
        "decision_duration_s": duration, "controller_power_w": power,
        "controller_multiplier": a, "controller_energy_j_per_decision": e_control,
        "controller_energy_j_per_token": nu * e_control,
        "engine_capacity_decisions_per_s": c["availability"] / duration,
        "rate_tokens_per_s": r, "engines_required": engines,
        "engines_per_plant": per_plant, "plants_required": plants,
        "controller_peak_w_all_engines": engines * power,
        "allocated_fixed_j_per_token": fixed,
        "controller_plus_fixed_j_per_token": nu * e_control + fixed,
        "scope": "controller contribution and optional assumed fixed plant allocation; other loads unquantified"
    }


def energy_gain(d, multiplier, f=0.30, h=0.05):
    finite(multiplier, "electrical multiplier", lower=1.0)
    finite(f, "energy fraction", upper=1.0, strict_lower=True)
    finite(h, "overhead fraction")
    d = np.asarray(d, dtype=float)
    if np.any(~np.isfinite(d)) or np.any(d <= 0):
        raise ValueError("displacement intensity must be finite and positive")
    return 1.0 / (1.0 - f + multiplier * f / d + h)


def reference_score(input_angles=None):
    """Synthetic, exactly specified eight-qubit score; no trained scheduler claim.

    Little-endian register. Input RY angles linearly spaced from -0.6 to 0.8.
    Four layers: RY(sin(1+8*l+q)/3), RZ(cos(1+8*l+q)/3) for every q,
    followed by CNOT(q,(q+1) mod 8) in increasing q order. Return <Z_0>.
    """
    n = 8
    angles = np.linspace(-0.6, 0.8, n) if input_angles is None else np.asarray(input_angles, dtype=float)
    if angles.shape != (n,) or not np.all(np.isfinite(angles)):
        raise ValueError("the reference circuit requires eight finite input angles")
    state = np.zeros(2**n, complex)
    state[0] = 1

    def one(gate, q):
        for i in range(2**n):
            if not i & (1 << q):
                j = i | (1 << q)
                state[i], state[j] = gate @ np.array([state[i], state[j]])

    def ry(angle):
        a, b = np.cos(angle/2), np.sin(angle/2)
        return np.array([[a, -b], [b, a]])

    def rz(angle):
        return np.diag([np.exp(-0.5j*angle), np.exp(0.5j*angle)])

    for q, x in enumerate(angles):
        one(ry(x), q)
    for layer in range(4):
        for q in range(n):
            one(ry(np.sin(1+8*layer+q)/3), q)
            one(rz(np.cos(1+8*layer+q)/3), q)
        for q in range(n):
            t = (q+1) % n
            for i in range(2**n):
                if i & (1 << q) and not i & (1 << t):
                    j = i | (1 << t)
                    state[i], state[j] = state[j], state[i]
    probabilities = np.abs(state)**2
    z = float(np.sum(probabilities * np.where(np.arange(2**n) & 1, -1, 1)))
    return {"input_angles_rad": angles.tolist(), "qubits": n, "amplitudes": 2**n, "single_qubit_gates": 72,
            "cnot_gates": 32, "norm": float(probabilities.sum()), "ideal_z0": z,
            "hardware_bias_certified": False, "classical_energy_measured": False}


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def figure_style():
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8,
                         "axes.labelsize": 9, "axes.titlesize": 9,
                         "legend.fontsize": 7, "xtick.labelsize": 8,
                         "ytick.labelsize": 8, "pdf.fonttype": 42,
                         "axes.spines.top": False, "axes.spines.right": False})


def save_figure(fig, out, stem):
    # Publish complete files only; interrupted rendering must not truncate a prior figure.
    pdf_tmp = out / (stem + ".tmp.pdf")
    png_tmp = out / (stem + ".tmp.png")
    fig.savefig(pdf_tmp, bbox_inches="tight", pad_inches=0.04,
                metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(png_tmp, dpi=400, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    pdf_tmp.replace(out / (stem + ".pdf"))
    png_tmp.replace(out / (stem + ".png"))


def draw_stack(out,c):
    fig, ax = plt.subplots(figsize=(7.05, 3.5))
    ax.set(xlim=(0, 10), ylim=(0, 6.3))
    ax.axis("off")
    def box(x, y, w, h, title, body, color):
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle="round,pad=0.04,rounding_size=0.05",
                                   facecolor=color,edgecolor="#52616b",linewidth=0.9))
        ax.text(x+w/2,y+h-0.18,title,ha="center",va="top",fontsize=9.2,fontweight="bold")
        body_y = y + (h/2-0.28 if h > 2 else 0.43)
        ax.text(x+w/2,body_y,body,ha="center",va="center",fontsize=8.0,linespacing=1.4)
    box(0.10,4.45,6.8,1.6,f"{c["hot_temperature_k"]:g} K: classical serving package",
        "Compute chiplets + HBM + I/O\nHost preparation, score verification, queue and fallback", "#f2f0ea")
    box(0.10,2.4,6.8,1.35,f"{c['controller_temperature_k']:g} K: control subsystem",
        f"Control allocation: {c['qubits']} x {c['controller_power_per_qubit_w']*1000:g} mW = {c['qubits']*c['controller_power_per_qubit_w']*1000:g} mW\nLocal electrical input and heat are accounted separately", "#e1ecf1")
    box(0.10,0.35,6.8,1.35,"Millikelvin stage: qubit assembly",
        "Qubit operations, readout absorption and wiring heat\nNo control-stage power is reassigned to this stage", "#cbdce7")
    for y in (3.79,1.75):
        ax.add_patch(FancyArrowPatch((.35,y+0.58),(.35,y+0.02),arrowstyle="-|>",mutation_scale=11,color="#9d3b31"))
    ax.text(0.75,4.07,"Thermally intercepted electrical or optical connections",fontsize=8.0,va="center")
    ax.text(0.75,2.02,"Signal energy and conductive heat: resolve by stage",fontsize=8.0,va="center")
    box(7.3,2.05,2.55,3.9,"Wall energy",
        "Electrical supplies\n+\nIncremental cooling\n+\nAllocated fixed plant\n\nA thermal stage is\nnot a qubit modality", "#fafafa")
    ax.text(0.15,0.02,"Conceptual multi-package system; dimensions and physical integration are not specified.",fontsize=8,va="bottom")
    save_figure(fig,out,"fig1_serving_system")


def draw_gain(out,c):
    fig,ax=plt.subplots(figsize=(3.42,2.85))
    d=np.logspace(-1,8,800)
    f,h=c["illustrative_energy_fraction"],c["illustrative_overhead_fraction"]
    styles=[(c["hot_temperature_k"],"Ambient","#1d6e85","-"),
            (4,"4 K","#287857","--"),
            (1,"1 K","#b57715",":"),
            (.02,"20 mK","#973831","-.")]
    for temp,label,color,ls in styles:
        if temp > c["hot_temperature_k"]:
            continue
        a=electrical_multiplier(temp,c["carnot_fraction"],c["hot_temperature_k"])
        ax.loglog(d,energy_gain(d,a,f,h),label=f"{label}, A = {a:,.6g}",color=color,ls=ls,lw=1.7)
        cross=a/(1-h/f)
        ax.scatter([cross],[1],color=color,s=14,zorder=4)
    ax.axhline(1,color="black",lw=.8)
    if 1-f+h > 0:
        ceiling=1/(1-f+h)
        ax.axhline(ceiling,color="#666666",lw=.8,ls=":")
    ax.set_title(f"f = {f:.2f}; H/Ecls = {h:.2f}",fontsize=8,pad=5)
    ax.set(xlim=(.1,1e8),ylim=(1e-4,2),xlabel="Displacement intensity D",ylabel="Energy gain Ecls / Ehybrid")
    ax.legend(loc="lower right",framealpha=.96)
    ax.grid(True,which="major",alpha=.2)
    fig.tight_layout(pad=.3)
    save_figure(fig,out,"fig2_cryogenic_roofline")


def draw_requirements(out,c):
    fig,axs=plt.subplots(2,1,figsize=(3.42,4.10))
    eps=np.geomspace(.03,.3,150)
    for eta,style,color in [(.01,"--","#973831"),(.02,"-","#1d6e85"),(.1,":","#287857")]:
        vals=[case_result(c,epsilon=x,eta=eta)["controller_energy_j_per_decision"] for x in eps]
        axs[0].loglog(eps,vals,style,color=color,lw=1.5,label=f"eta = {eta:g}")
    axs[0].set(xlabel="Statistical half-width epsilon",ylabel="Controller J/decision",title="(a) Precision sets the shot budget")
    axs[0].set_xticks([.03,.05,.1,.2,.3],labels=["0.03","0.05","0.10","0.20","0.30"])
    axs[0].tick_params(axis="x",which="minor",labelbottom=False)
    axs[0].legend(loc="upper right",fontsize=7)
    rate=np.geomspace(10,10000,1200)
    vals=[case_result(c,rate=x) for x in rate]
    axs[1].loglog(rate,[x["controller_plus_fixed_j_per_token"] for x in vals],color="#1d6e85",lw=1.6,label="Controller + assumed fixed plant")
    axs[1].loglog(rate,[x["controller_energy_j_per_token"] for x in vals],color="#287857",ls="--",lw=1.2,label="Work-proportional controller")
    axs[1].set(xlabel="Completed tokens/s",ylabel="Allocated J/token",title=f"(b) Sharing with a {c["controller_temperature_k"]:g} K capacity limit")
    axs[1].legend(loc="upper right",fontsize=6.5)
    for ax in axs:
        ax.grid(True,which="major",alpha=.2)
    fig.tight_layout(pad=.35,h_pad=1.2)
    save_figure(fig,out,"fig3_precision_capacity")


def build(output_dir=ROOT, config_path=CONFIG_PATH):
    out=Path(output_dir)
    out.mkdir(parents=True,exist_ok=True)
    c=json.loads(Path(config_path).read_text())
    validate_config(c)
    case=case_result(c)
    rows=[case_result(c,epsilon=e) for e in (.2,.1,.05)]
    stage_rows=[{"temperature_k":t,"eta":e,"cooling_factor_chi":cooling_factor(t,e,c["hot_temperature_k"]),
                 "local_electrical_multiplier_A":electrical_multiplier(t,e,c["hot_temperature_k"]),
                 "status":"illustrative efficiency; no complete-system prediction"}
                for t,e in [(c["hot_temperature_k"],c["carnot_fraction"]),(4,c["carnot_fraction"]),
                            (1,c["carnot_fraction"]),(.02,c["carnot_fraction"])] if t <= c["hot_temperature_k"]]
    write_csv(out/"table_cooling_stages.csv",stage_rows)
    write_csv(out/"table_case_inputs.csv",[{**c,"scope":"single controller-allocation case; see parameter_provenance.csv"}])
    write_csv(out/"table_controller_requirements.csv",rows)
    write_csv(out/"table_capacity.csv",[case_result(c,rate=x) for x in (10,100,1000,3000,3200,10000)])
    provenance=[]
    for k,v in c.items():
        status="assumed analysis input"
        source="artifact/config.json; not a measured LLM or complete quantum system"
        if k=="controller_power_per_qubit_w":
            status=("published measurement anchor; linear scaling assumed" if v == .023 and c["controller_temperature_k"] == 4 else "scenario override; not the published operating point")
            source="Underwood et al., PRX Quantum 5, 010326 (2024), abstract; 14 nm, 4 K, active control"
        elif k=="controller_temperature_k":
            status=("temperature of published controller anchor" if v == 4 else "assumed alternative controller temperature; hardware not validated")
            source="Underwood et al., DOI 10.1103/PRXQuantum.5.010326"
        elif k=="artifact_version":
            status="artifact metadata"
        provenance.append({"parameter":k,"value":v,"status":status,"source_or_qualification":source})
    write_csv(out/"parameter_provenance.csv",provenance)
    tex=[]
    macros={"CaseShots":str(case["shots"]),"CasePowerMW":f'{case["controller_power_w"]*1e3:.0f}',
            "CaseEnergy":f'{case["controller_energy_j_per_decision"]:.2f}',
            "CaseTokenEnergy":f'{case["controller_energy_j_per_token"]:.2f}',
            "CaseDurationMS":f'{case["decision_duration_s"]*1e3:.1f}',
            "CaseEngineRate":f'{case["engine_capacity_decisions_per_s"]:.2f}',
            "CaseEngines":str(case["engines_required"]),
            "CasePeakPower":f'{case["controller_peak_w_all_engines"]:.3f}',
            "CaseSubtotal":f'{case["controller_plus_fixed_j_per_token"]:.2f}',
            "CasePlants":str(case["plants_required"]),
            "EnginesPerPlant":str(case["engines_per_plant"])}
    for k,v in macros.items():
        tex.append("\\newcommand{\\"+k+"}{"+v+"}")
    (out/"generated_values.tex").write_text("% Generated by artifact/cryogenic_model.py\n"+"\n".join(tex)+"\n")
    table=[]
    for r in rows:
        table.append(f'{r["epsilon_stat"]:.2f} & {r["shots"]:,} & {r["decision_duration_s"]*1000:.1f} & {r["controller_energy_j_per_decision"]:.2f} & {r["engines_required"]} '+r"\\")
    (out/"table_case_rows.tex").write_text("% Generated by artifact/cryogenic_model.py\n"+"\n".join(table)+"\n")
    table_source = (r"\begin{tabular}{@{}rrrrr@{}}" + "\n" + r"\toprule" + "\n"
                    + r"$\epsilon$ & Shots & ms/decision & J/decision & Engines \\" + "\n"
                    + r"\midrule" + "\n" + "\n".join(table) + "\n"
                    + r"\bottomrule" + "\n" + r"\end{tabular}" + "\n")
    (out/"table_case.tex").write_text(table_source)
    reference=reference_score()
    reference["scope"]="fixed eight-qubit verification circuit, independent of controller sizing overrides"
    (out/"reference_inputs.json").write_text(json.dumps([reference_score(np.linspace(-.6,.8,8)+offset)
                                                       for offset in (-.2,0,.2)],indent=2)+"\n")
    (out/"reference_score.json").write_text(json.dumps(reference,indent=2)+"\n")
    (out/"model_results.json").write_text(json.dumps({"version":c["artifact_version"],"case":case,"reference":reference},indent=2)+"\n")
    figure_style()
    draw_stack(out,c)
    draw_gain(out,c)
    draw_requirements(out,c)
    return case


if __name__=="__main__":
    print(json.dumps(build(),indent=2))
