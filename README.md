# Cryogenic Roofline

Code and numerical data accompanying **The Cryogenic Roofline: Energy Break-Even for Chiplet-Based Quantum-Assisted LLM Serving**, by **Somdip Dey** ([ORCID](https://orcid.org/0000-0001-6161-4637)).

**Artifact version:** 2.1.0  
**Manuscript status:** Unpublished manuscript, 2026. Not accepted for publication yet.

This repository implements a conditional energy and capacity analysis of quantum-assisted LLM serving. It combines stage-resolved refrigeration accounting, precision-dependent measurement budgets, controller energy, and capacity-limited allocation of fixed plant power. A deterministic eight-qubit state-vector calculation provides an ideal reference score.

The numerical results are analytical scenarios, not measurements of an LLM serving system. They do not establish quantum advantage, an application-quality improvement, or a ranking of qubit technologies. The controller-power anchor is 23 mW per qubit under active control at 4 K, reported by [Underwood et al., PRX Quantum 5, 010326 (2024)](https://doi.org/10.1103/PRXQuantum.5.010326). Linear scaling, shot duration, refrigeration efficiency, scheduling, and plant allocations are separate assumptions recorded in the provenance table.

## Repository contents

Keep the following paths. The four files in `artifact/` must use the names shown, without download suffixes such as `(1)`. The notebook expects this directory layout.

| Path | Purpose |
| --- | --- |
| `README.md` | Installation, reproduction, interpretation, and citation instructions. |
| `artifact/cryogenic_model.py` | Shared energy, precision, capacity, state-vector, and output-generation implementation. |
| `artifact/config.json` | Default scenario and artifact version. This is the configuration read by the build. |
| `artifact/requirements.txt` | Pinned NumPy and Matplotlib dependencies. |
| `artifact/validate.py` | Numerical and physical-domain checks, including an independently implemented reference simulator. |
| `cryogenic_roofline_calculator.ipynb` | Interactive explanation and execution of the shared model. |
| `model_results.json` | Default controller/capacity result and ideal reference-circuit result. |
| `parameter_provenance.csv` | Parameter values, evidence status, and source or qualification. |
| `table_case_inputs.csv` | Snapshot of the configuration used to generate the outputs. |
| `table_cooling_stages.csv` | Cooling-work factors and total local electrical multipliers at the tabulated temperatures. |
| `table_controller_requirements.csv` | Precision sweep at statistical half-widths 0.20, 0.10, and 0.05. |
| `table_capacity.csv` | Capacity and allocated energy at 10, 100, 1,000, 3,000, 3,200, and 10,000 tokens/s. |
| `reference_score.json` | Default eight-qubit input angles, ideal score, normalization, and gate counts. |
| `reference_inputs.json` | Three input vectors and their corresponding ideal scores, generated using angular offsets of −0.2, 0, and +0.2 radians. |

The CSV and JSON files are generated numerical outputs and verification fixtures. They are not raw experimental data. The state-vector implementation is `reference_score()` in `artifact/cryogenic_model.py`; the JSON files contain its results, not executable implementations or full state vectors.

## Installation

Use **Python 3.12**. The code uses Python 3.12 f-string syntax and will not parse on Python 3.11 or earlier. The supplied dependency versions are NumPy 2.3.5 and Matplotlib 3.10.8. No quantum SDK, accelerator, or quantum hardware is required.

Download or clone this repository, open a terminal in its root directory, and create an isolated environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r artifact/requirements.txt
```

On Windows, create the environment with `py -3.12 -m venv .venv`, activate it with `.venv\Scripts\Activate.ps1` in PowerShell, and run the same pip command.

## Reproduce the outputs

From the repository root:

```bash
python artifact/cryogenic_model.py
python artifact/validate.py
```

The first command reads `artifact/config.json`, prints the current case as JSON, and regenerates the supplied CSV and JSON outputs at the repository root. **Existing generated files are overwritten.** Keep the default configuration and its outputs in version control before exploring alternatives.

The build also generates the following files locally; they do not need to be present before running it:

- `fig1_serving_system.pdf` and `.png`: system and accounting-boundary schematic.
- `fig2_cryogenic_roofline.pdf` and `.png`: conditional energy-gain curves.
- `fig3_precision_capacity.pdf` and `.png`: precision and capacity requirements.
- `generated_values.tex`, `table_case.tex`, and `table_case_rows.tex`: numerical inputs for manuscript preparation.

No LaTeX installation is needed to generate these files. Matplotlib uses a non-interactive backend, so the script does not require a graphical display.

With the supplied configuration, representative results are:

| Quantity | Default result |
| --- | ---: |
| Sufficient shots per decision | 738 |
| Decision duration | 0.0738 s |
| Active controller power | 0.184 W |
| Local electrical multiplier at 4 K | 3,701 |
| Controller energy per decision, including modeled cooling | 50.2566192 J |
| Controller energy per token | 1.57051935 J |
| Engines required at 1,000 tokens/s | 3 |
| Engines per plant / plants required | 8 / 1 |
| Controller plus allocated fixed plant energy per token | 6.57051935 J |
| Ideal reference score, ⟨Z₀⟩ | approximately 0.507120883394 |

Small floating-point differences may occur across environments. These energy values are partial modeled contributions, not complete system energy or measured savings.

Successful validation prints:

```text
PASS: energy accounting, shot guarantee, crossover, capacity, reference circuit and input domains
```

The checks cover direct-input and cooling accounting, the sufficient shot bound, energy crossover, engine and plant capacity transitions, invalid inputs, and agreement between two state-vector implementations for three input vectors. Validation checks the active configuration and separately uses canonical default values for regression checks. It does not compare every checked-in CSV/JSON file against regenerated output, certify hardware behavior, or establish application-level performance. Run it without Python's `-O` option, which disables assertions.

### Run the notebook

The notebook uses the same model and configuration. To use JupyterLab, install the optional interface in the same environment:

```bash
python -m pip install jupyterlab
python -m jupyterlab cryogenic_roofline_calculator.ipynb
```

Launch JupyterLab from the repository root, select the environment containing the installed dependencies, and run all cells in order. The regeneration cell writes the same outputs as the command-line build. JupyterLab is optional and is not pinned in `artifact/requirements.txt`.

## Model interpretation

For a stage at temperature `Tc`, with hot-side temperature `Th` and efficiency `eta` expressed as a fraction of Carnot coefficient of performance:

```text
chi = (Th - Tc) / (eta * Tc)
A   = 1 + chi
```

`chi` is refrigeration work per joule of heat lifted. `A` is the total thermal amplification factor for locally supplied electrical energy fully dissipated at that stage, including the direct electrical input. At ambient temperature, `A = 1` within this refrigeration boundary; it does not include facility cooling.

The measurement budget is:

```text
shots = ceil(2 * ln(2 / delta) / epsilon^2)
```

This is a sufficient Hoeffding bound for independent outcomes in `[-1, 1]`, with statistical error `epsilon` and failure probability `delta`. It is not a measured minimum shot count. An error bound against the ideal score also requires a separately certified hardware-bias bound.

Controller energy is active controller power multiplied by decision duration and `A`. The model rounds engine counts upward to meet average service demand and limits the number of engines in each plant by the configured controller cooling allocation. Fixed plant power is allocated per completed token and recurs when additional plants are required.

Other host, idle-control, wiring, detection, calibration, error-management, and cryogenic-stage costs remain unquantified. Fixed plant power must exclude contributions already charged through the incremental cooling factor. Capacity sufficiency does not guarantee tail latency or feasibility at every temperature stage.

## Explore alternative scenarios

### Change configuration parameters

Edit `artifact/config.json`, then rerun the model. Change one parameter at a time when examining its effect. The following examples start from the default configuration:

| Parameter | Default → example | Purpose and expected effect |
| --- | --- | --- |
| `statistical_error` | `0.1` → `0.05` | Require tighter statistical precision. Shots increase from 738 to 2,952; decision duration and controller energy increase fourfold. |
| `failure_probability` | `0.05` → `0.01` | Require a lower statistical failure probability. The sufficient shot budget increases. |
| `shot_duration_s` | `0.0001` → `0.00005` | Explore a hypothetical faster complete shot. Controller energy and decision duration halve; at the default rate the required engines fall from three to two. Hardware feasibility must be established separately. |
| `carnot_fraction` | `0.02` → `0.1` | Explore more efficient refrigeration. At 4 K and a 300 K hot side, `A` falls from 3,701 to 741; controller wall energy decreases without changing shots or execution duration. |
| `token_rate_per_s` | `1000.0` → `3200.0` | Examine demand beyond one plant's modeled capacity. Required engines increase to nine and plants to two. Fixed energy per token changes in steps as plants are added. |
| `tokens_per_decision` | `32` → `64` | Assume fewer offload decisions per token. Controller energy per token halves and service demand falls. This changes the workload contract and requires an application justification. |
| `availability` | `0.9` → `0.5` | Reduce usable engine capacity. More engines may be needed; energy per completed decision is unchanged. |
| `controller_power_per_qubit_w` | `0.023` → `0.01` | Explore a hypothetical lower-power controller. Active energy falls and more engines fit within the cooling allocation. The generated provenance marks this as a scenario override. |
| `available_controller_cooling_w_per_plant` | `1.5` → `1.0` | Tighten the controller cooling allocation. Fewer engines fit in each plant, moving replication thresholds. |
| `fixed_wall_power_w_per_plant` | `5000.0` → `2500.0` | Explore a lower separately allocated fixed cost. This halves the fixed-energy contribution without changing controller energy or capacity. |
| `illustrative_energy_fraction` | `0.3` → `0.5` | Increase the assumed fraction of classical energy displaced in the energy-gain plot. This affects Figure 2, not the controller case. |
| `illustrative_overhead_fraction` | `0.05` → `0.1` | Increase normalized additional overhead in Figure 2. Break-even requires greater displacement intensity. This also leaves the controller case unchanged. |

Temperature inputs are `hot_temperature_k` and `controller_temperature_k`. Changing them explores the thermal model, not validated operation of the published controller at a new temperature. `carnot_fraction` must lie in `(0, 1]`, and `0 < controller_temperature_k <= hot_temperature_k`. The controller cooling allocation must support at least one engine. `tokens_per_decision` and `qubits` must be positive integers, `availability` must lie in `(0, 1]`, and the plotted overhead fraction must be smaller than the displaced energy fraction.

Two configuration fields require particular care:

- `allowed_bias` records an assumed bias allowance. It does **not** alter shot counts, energy, or the ideal circuit calculation. Changing it does not certify a hardware-bias bound.
- `qubits` changes linear controller sizing and capacity. The reference simulator remains fixed at eight qubits. Increasing this parameter does not produce or validate a larger reference circuit.

The build contains explicit sweep grids: the precision and rate values listed in the file table above, Figure 3's efficiency curves at 0.01, 0.02, and 0.10, and representative temperatures in the cooling-stage table and Figure 2. Changing a default parameter does not replace these grids. To change the sweep coordinates themselves, edit the corresponding lists in `build()`, `draw_requirements()`, or `draw_gain()` in `artifact/cryogenic_model.py`. Some notebook cells likewise show fixed illustrative examples rather than configuration-driven sweeps.

### Evaluate a scenario without overwriting the supplied data

Run the following in a Python session or notebook started at the repository root:

```python
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path("artifact").resolve()))
import cryogenic_model as model

config = json.loads(Path("artifact/config.json").read_text())
custom = dict(config)
custom["shot_duration_s"] = 50e-6

result = model.case_result(custom)
print(json.dumps(result, indent=2))
```

This evaluates the modified dictionary in memory. It does not edit the configuration or regenerate files. To save a separate scenario and all its outputs:

```python
scenario_dir = Path("results/shorter_shot")
scenario_dir.mkdir(parents=True, exist_ok=True)
scenario_config = scenario_dir / "config.json"
scenario_config.write_text(json.dumps(custom, indent=2) + "\n")

model.build(output_dir=scenario_dir, config_path=scenario_config)
```

The saved scenario configuration records the actual input for this run. Generated provenance uses a generic `artifact/config.json` source label for assumed parameters, so retain the scenario configuration alongside the outputs when using an alternative path. Calling `model.build()` without `config_path` reads the original configuration, not the in-memory `custom` dictionary.

### Change the reference-circuit input

After importing the shared model as above:

```python
import numpy as np

angles = np.linspace(-0.6, 0.8, 8) + 0.2  # radians
reference = model.reference_score(angles)
print(reference["ideal_z0"])  # approximately 0.443847250810
```

`reference_score()` accepts exactly eight finite angles in radians. It applies eight input RY gates, followed by four layers of fixed RY/RZ rotations and a sequential CNOT ring, and returns the ideal expectation of Z on qubit 0. The circuit uses 256 complex amplitudes, 72 single-qubit gates, and 32 CNOT gates, with little-endian qubit indexing.

Changing the input angles changes the score but not the gate count or the current controller energy estimate. The energy model does not derive shot duration or power from circuit simulation. The supplied reference is classically tractable, has no trained parameters, and is a verification task rather than an evaluated LLM routing policy.

## Citation

If you use this code or its numerical data, cite the accompanying unpublished manuscript and identify the repository commit or release used. The following is a provisional citation; it does not imply acceptance or publication in IEEE Micro:

> Somdip Dey. “The Cryogenic Roofline: Energy Break-Even for Chiplet-Based Quantum-Assisted LLM Serving.” Unpublished manuscript, 2026. Accompanying computational artifact, version 2.1.0.

```bibtex
@unpublished{dey2026cryogenic_roofline,
  author = {Dey, Somdip},
  title  = {The Cryogenic Roofline: Energy Break-Even for
            Quantum-Assisted {LLM} Serving},
  year   = {2026},
  note   = {Unpublished manuscript. Accompanying computational
            artifact, version 2.1.0}
}
```

