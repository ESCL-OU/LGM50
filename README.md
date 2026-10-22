# LGM50
DFN | SPMe | ECM Simulation of LGM50

This repository is based on and extends the DFN implementation in [fastDFN](https://github.com/scott-moura/fastDFN) and the SPMe implementation in [SPMeT](https://github.com/scott-moura/SPMeT), with only minor modifications to those original codes. We gratefully acknowledge and thank **Dr. Scott Moura** and collaborators for making these high-quality reference implementations publicly available.

## License

This repository incorporates code from the SPMeT project, which is licensed under the GNU General Public License version 3 (GPL‑3.0). In accordance with the GPL‑3.0’s copyleft requirements, **this project is distributed under the terms of the GNU General Public License v3.0**, as provided in the `LICENSE` file. By using, modifying, or redistributing this code, you agree to the terms of the GPL‑3.0.


## Modifications relative to fastDFN and SPMeT

- **DFN implementation (based on `fastDFN` [3])**
  - **Initialization via `init_cs` (Scott Moura style)**: initial electrode concentrations are computed with `init_cs(p, V)` from `param/params_LGM50.m`.
  - **Added SoC-to-voltage initializer**: this repository adds `calculate_equilibrium_voltage_from_soc.m`, which converts target SoC to equilibrium voltage, then uses `init_cs` for concentration initialization. In the original Scott Moura flow, users typically specify the initial voltage directly.
  - **Thermal dynamics switch**: a parameter flag (e.g., `disable_thermal`) has been added so that the thermal dynamics in the DFN model from [3] can be turned on or off, enabling both isothermal and electrochemical–thermal simulations with a single code base.

- **SPMe implementation (based on `SPMeT` [4])**
  - **Thermal dynamics removed**: only the electrochemical SPMe model from the original SPMeT repository is used; its thermal dynamics are neglected here.
  - **CasADi-based class formulation**: the SPMe model has been refactored into a CasADi-based MATLAB class `SPMe/SPMe.m` that encapsulates state dynamics, single-step integration, voltage outputs, voltage-term decomposition, and optional diagnostic outputs.
  - **Configurable degradation models**: SPMe can run without degradation or with SEI growth models based on Ramadass, Yang, or Safari. Lithium plating and anode porosity-loss add-ons can be enabled independently for SEI-enabled SPMe simulations.

# Install Requirements

We make use of [CasADi](https://web.casadi.org/) for the symbolic variables to allow models to be utilized in MPC and optimal control methods. 
So to run the simulations you will need to install [CasADi](https://web.casadi.org/) framework.
To install [CasADi](https://web.casadi.org/):

1. Download the precomplied binary from the official CasADi website's [get page](https://web.casadi.org/get/)
2. Extract and place the files
Unzip the downloaed file to a permanent location on your computer.
3. Add the directory of the unzipped folder to the MATLAB Path
```
addpath('<yourpath>/casadi-3.7.2-windows64-matlab2018b')
```
Optionally you can also make this change permanent so you don't have to run `addpath` every time you start MATLAB. Go to **Home** tab, click **Set Path**, and add the folder there.

5. Verify Installation

```
import casadi.*
x = MX.sym('x')
disp(jacobian(sin(x),x))
```
The output should display `cos(x)` to confirm that CasADi functions are working correctly.


6. Point the project at your CasADi folder

`startup.m` adds CasADi and `tests/` to the MATLAB path. Edit its hard-coded `casadi_dir` to your CasADi location, then run `startup` from the repository root at the start of each session.

# Running simulations

The main simulation entry points are:

- `sim_profile_LGM50.m`: run one profile segment.
- `sim_profile_sequence_LGM50.m`: run concatenated profile segments continuously, with state handoff between segments.

Both functions use SoC-based initialization:

1. `calculate_equilibrium_voltage_from_soc(init_soc)` computes an equilibrium voltage.
2. `init_cs` computes initial electrode concentrations from that voltage.
3. The selected model builds its state vector from those concentrations, unless an `InitialState` is provided for continuation.

## Supported Models And Profiles

Supported model names are:

- `SPMe`
- `DFN`
- `ECM`

Supported profile definitions are:

- Dynamic drive profiles: `BJDST`, `DST`, `FUDS`, `US06`
- Constant-current profile: `{'CC', crate}`
- Fixed-duration constant-current segment (sequences only): `{'CC', crate, duration_hours}`
- Rest profile: `{'REST', duration_hours}`

Current sign convention:

- `crate > 0`: discharge
- `crate < 0`: charge

Concatenated sequences are passed as token lists:

```matlab
{'REST', 0.25, 'CC', -1, 'CC', 1}
```

## Quick Start

```matlab
% One SPMe dynamic segment
sim = sim_profile_LGM50('SPMe', 0.4, 'DST');

% One DFN 1C charge segment
sim = sim_profile_LGM50('DFN', 0.4, {'CC', -1});

% One ECM dynamic segment, saved to disk
sim = sim_profile_LGM50('ECM', 0.4, 'FUDS', 'SaveOutput', true);

% Repeated SPMe charge/discharge sequence
sim = sim_profile_sequence_LGM50('SPMe', 0.1, {'CC', -1, 'CC', 1}, ...
    'Repeat', 5);

% Sequence with rests and a voltage window
sim = sim_profile_sequence_LGM50('SPMe', 0.1, ...
    {'REST', 0.25, 'CC', -0.5, 'REST', 0.25, 'CC', 0.5}, ...
    'Repeat', 3, ...
    'VMin', 2.5, 'VMax', 4.2);
```

## SPMe Degradation Options

SPMe degradation is configured with `DegradationModel`:

- `none`: no degradation state.
- `ramadass`: Ramadass-style SEI growth.
- `yang`: Yang-style SEI growth.
- `safari`: Safari-style SEI growth with EC diffusion states inside the SEI film.

Lithium plating and anode porosity loss are independent add-ons. They are disabled by default for all SEI models and are ignored when `DegradationModel` is `none`.

```matlab
sim = sim_profile_sequence_LGM50('SPMe', 0.1, {'CC', -1, 'CC', 1}, ...
    'Repeat', 10, ...
    'DegradationModel', 'safari', ...
    'EnableLithiumPlating', true, ...
    'EnableAnodePorosityLoss', true, ...
    'SPMeOutputs', {'delta_sei', 'delta_li', 'delta_film', ...
                    'epsilon_e_n', 'j_sei', 'j_plating', 'j_side'});
```

Useful SPMe optional output names include:

- Degradation states and rates: `delta_sei`, `delta_sei_dot`, `delta_li`, `delta_li_dot`, `delta_film`, `delta_film_dot`
- Porosity and film resistance: `epsilon_e_n`, `epsilon_e_n_dot`, `R_film`
- Side-reaction terms: `j_sei`, `j_plating`, `j_side`, `i0s`
- Side-reaction potentials: `eta_n_side`, `phi_se`, `eta_s`
- Electrolyte and averaged solid diagnostics: `cen_bar`, `ces_bar`, `cep_bar`, `cs_ave_n`, `cs_ave_p`

## Single-Segment Options

`sim_profile_LGM50` accepts one profile definition and returns one simulation segment.

```matlab
sim = sim_profile_LGM50(model_name, init_soc, profile_def, ...
    'Nxn', 20, 'Nxs', 10, 'Nxp', 20, 'Nr', 20, ...
    'DeltaT', 1, 'Tamb', 25, ...
    'DegradationModel', 'none', ...
    'ElectrolyteModel', 'fdm', ...
    'UseElectrolyteActivity', false, ...
    'EnableLithiumPlating', false, ...
    'EnableAnodePorosityLoss', false, ...
    'InitialState', struct(), ...
    'SPMeOutputs', {}, ...
    'SaveOutput', false, 'OutputFile', '');
```

Important options:

- `Nxn`, `Nxs`, `Nxp`, `Nr`: discretization settings for electrochemical models.
- `DeltaT`: fixed simulation timestep `[s]`.
- `Tamb`: ambient temperature `[degC]`.
- `ECMOrder`: number of RC branches for ECM.
- `DegradationModel`: SPMe degradation model.
- `ElectrolyteModel`: SPMe electrolyte dynamics, `fdm` or `fvm`.
- `UseElectrolyteActivity`: include electrolyte activity correction in SPMe electrolyte-polarization voltage. Default is `false`.
- `InitialState`: continue from a previous terminal state, e.g. `struct('x', x_last)` for SPMe.
- `SPMeOutputs`: optional SPMe diagnostics to store in `sim.spme_outputs`.
- `VEnd`: custom voltage termination for a segment.
- `SaveOutput` / `OutputFile`: save `sim` to a `.mat` file.

## Sequence Options

`sim_profile_sequence_LGM50` runs segments continuously by passing the terminal state of each segment into the next segment.

```matlab
sim = sim_profile_sequence_LGM50(model_name, init_soc, profile_sequence, ...
    'Repeat', 1, ...
    'Nxn', 10, 'Nxs', 5, 'Nxp', 10, 'Nr', 15, ...
    'DeltaT', 1, 'Tamb', 25, ...
    'DegradationModel', 'none', ...
    'ElectrolyteModel', 'fdm', ...
    'EnableLithiumPlating', false, ...
    'EnableAnodePorosityLoss', false, ...
    'VMin', [], 'VMax', [], ...
    'ShowProgress', true, ...
    'MemoryLight', false, ...
    'InitialState', struct(), ...
    'CheckpointFile', '', ...
    'SPMeOutputs', {}, ...
    'SaveOutput', false, 'OutputFile', '');
```

Additional sequence-specific options:

- `Repeat`: repeat the full `profile_sequence`.
- `VMin` / `VMax`: truncate each segment if voltage leaves the requested window.
- `ShowProgress`: show a waitbar.
- `SuppressSummary`: suppress the per-segment progress printout.
- `UseElectrolyteActivity`: forwarded to every segment (same meaning as in `sim_profile_LGM50`).
- `MemoryLight`: omit large history fields from the aggregate sequence output while keeping core signals and requested SPMe outputs.
- `CheckpointFile`: save a warm-start checkpoint at the end of each completed repeat. The checkpoint contains the aggregate `sim`, the last segment, terminal state, terminal SOC/voltage, and sequence metadata.

## Warm Starts And Checkpoints

For manual continuation, pass a terminal state through `InitialState`:

```matlab
sim1 = sim_profile_LGM50('SPMe', 0.2, {'CC', -1}, ...
    'DegradationModel', 'safari');

state = struct('x', sim1.x(:, end));

sim2 = sim_profile_LGM50('SPMe', sim1.soc(end), {'CC', 1}, ...
    'DegradationModel', 'safari', ...
    'InitialState', state);
```

For long sequence runs, use `CheckpointFile`:

```matlab
sim = sim_profile_sequence_LGM50('SPMe', 0.1, {'CC', -1, 'CC', 1}, ...
    'Repeat', 100, ...
    'DegradationModel', 'safari', ...
    'EnableLithiumPlating', true, ...
    'EnableAnodePorosityLoss', true, ...
    'CheckpointFile', 'results/SPMe/long_run_checkpoint.mat');
```

The saved checkpoint can be loaded and used to construct an `InitialState` for later continuation.

## Outputs

Simulation functions return a struct `sim`.

Row timing (SPMe, since 2026-07-04): row 1 is the initial state with `current = 0`; row k (k > 1) holds the state, voltage, SoC and optional outputs after integrating one `DeltaT` step with the current stored in that same row. Data produced before this fix stored the current one row ahead of the state it produced.

Common fields:

- `time`, `current`, `voltage`, `soc`
- `model`, `profile`, `p`
- `terminated_early`, `termination_reason`, `t_last`
- `init_soc_input`, `init_voltage_from_soc`, `init_csn0`, `init_csp0`

Sequence-only fields:

- `profile_definition`, `profile_sequence_input`, `repeat_count`
- `segment_index`, `segment_cycle`, `capacity_cycle`
- `segment_log`, `capacity_cycle_log`

Model-specific fields:

- SPMe: `x`, `y`, `c_s_n`, `c_s_p`, `c_e`, and optional degradation fields such as `delta_sei`, `delta_li`, `epsilon_e_n`
- DFN: `x`, `z`, `y`, `c_s_n_pade`, `c_s_p_pade`, `c_e`, `temperature`
- ECM: `x`, `ocv_p`, `ecm_p`

If `SPMeOutputs` is requested, additional named time series are stored under `sim.spme_outputs`.


# PINN-SPMe dataset generation

This repository generates the synthetic datasets used by [PINN-SPMe](https://github.com/yuichikajiura/PINN-SPMe), which estimates electrochemical states and SEI thickness from current/voltage data. The pipeline has three MATLAB steps.

**1. Simulate a degraded cell (hours)**

```matlab
startup
out = test_spme_degradation_addons_sequence_plot([], false, 'fvm');
```

Settings are at the top of `tests/test_spme_degradation_addons_sequence_plot.m`: Safari SEI, anode porosity loss on, Li plating off, start at SoC 0.1, Nr 15, Nxn/Nxs/Nxp 20/5/20. One repeat = one C/5 cycle + 99 1C cycles; 20 repeats. Each completed repeat writes `results/SPMe/spme_degradation_addons_checkpoint_safari_<timestamp>_repeat####_segment#####.mat` plus a `_manifest.mat`. Pass a checkpoint path as the first argument to resume. Plot SEI, porosity, capacity and Li inventory with `tests/plot_spme_degradation_addons_from_data.m`.

The reference degraded cell used for the PINN-SPMe datasets is **repeat 15 of run `20260702_034413`** (about 1,500 cycles, `delta_sei` = 5.05188e-7 m, `epsilon_e_n` = 0.0560).

**2. Generate the datasets**

```matlab
% fresh cell -> results/PINN_SPMe_dataset/
generate_pinn_spme_training_data_LGM50_v5

% degraded cell -> results/PINN_SPMe_dataset/degraded_from_<checkpoint name>/
generate_pinn_spme_training_data_LGM50_v5('results/SPMe/spme_degradation_addons_checkpoint_safari_20260702_034413_repeat0015_segment03030.mat')
```

The generator runs SPMe (FVM electrolyte, Nr 15, Nxn/Nxs/Nxp 10/5/10, 1 s steps) for BJDST, FUDS, US06 and DST drive cycles plus a 1C discharge, all from SoC 0.9 to 2.5 V, and five HPPC sequences (120 s rest, 10 s 2C discharge pulse (1C at SoC 0.3), 120 s rest, 10 s 1C charge pulse, 120 s rest) at SoC 0.9/0.75/0.6/0.45/0.3. Each job writes `<name>.mat` (full `sim`, `data`, options) and `<name>.csv` (time, current, voltage, SoC, surface/average concentrations and the full state in PINN order `[c_s_n, c_e, c_s_p]`), plus `dataset_summary.csv`.

In degraded mode only the degradation states are copied from the checkpoint (SEI thickness, plated Li, anode porosity, Safari EC profile); each job still starts at its own SoC. The Li lost to SEI and plating is removed from the anode baseline (`p.csn0`, `p.n_Li_s`) before `init_cs`, so the degraded cell's cyclable Li is reduced.

Quick look at generated profiles: `plot_pinn_spme_generated_profiles('results/PINN_SPMe_dataset', true)`.

**3. Copy into PINN-SPMe**

Copy the CSVs into `PINN-SPMe/data/` with this rename: fresh `X.csv` -> `spme_X.csv`; degraded `X.csv` -> `spme_X_soh90.csv`; `id_HPPC_*.csv` -> `id_HPPC_*_soh90.csv` ("soh90" is a label only).

## Repository map

| File / folder | Role | Status |
|---|---|---|
| `sim_profile_LGM50.m`, `sim_profile_sequence_LGM50.m` | Simulation entry points | Current |
| `SPMe/SPMe.m`, `SPMe/load_SPMe_params.m`, `param/` | CasADi SPMe class with SEI models; LG M50 parameters, `init_cs` | Current |
| `generate_pinn_spme_training_data_LGM50_v5.m` | PINN-SPMe dataset generator | Current |
| `tests/test_spme_degradation_addons_sequence_plot.m`, `tests/plot_spme_degradation_addons_from_data.m` | Degraded-cell cycling run and its plots | Current |
| `plot_pinn_spme_generated_profiles.m`, `calculate_equilibrium_voltage_from_soc.m`, `startup.m` | Helpers | Current |
| `DFN/`, `sim_DFN_*.m`, `tests/*dfn*`, `run_test_dfn_spatial_batch.m`, `matlab_batch_dfn_spatial.sh` | DFN model and spatial-degradation studies | Separate line of work; not used by PINN-SPMe |
| other `tests/*` | Model validation (PyBaMM comparisons, mass balance, SEI model comparisons) | Diagnostic |
| `generate_pinn_spme_training_data_LGM50.m`, `_v2`, `_v3`, `_v4` | Earlier generators (no HPPC; path and sequence-profile bugs) | Legacy; keep for history only |
| `SPMe/SPMe_legacy.m`, `sim_SPMe_CC.m`, `sim_SPMe_CCCV.m`, `sim_DFN_CC.m`, `sim_DFN_CCCV.m`, `sim_ECM_CC.m`, `ECM 2/` | Older single-purpose scripts | Legacy |
| `matlab_batch.sh`, `run_test_spme_batch.m` | SLURM wrapper for the old Ramadass-vs-Yang test | Legacy |

`results/` is git-ignored. Under `results/SPMe/`, the full degradation runs are `20260625_204743` (FVM, Li plating on) and `20260702_034413` (FVM, plating off, used for PINN-SPMe); the other runs, `*_instability_diagnostics_*` and `*_warmstart_failure_trace_*` files are debugging outputs.


# LG M50 cell description

The models in this repository are parameterized for **LG’s M50** lithium‑ion cell, a high‑energy‑density **21700‑format** cylindrical cell (approx. 21 mm diameter, 70 mm length) with a nominal capacity of about **5 Ah** and nominal voltage around **3.6–3.7 V**.

- **Chemistry and use case**
  - The cell employs a **graphite negative electrode** and a layered oxide **nickel‑rich positive electrode** (NMC‑type), representative of modern high‑energy automotive cells.
  - It is designed for **energy‑focused applications** (e.g., EV packs, stationary storage) where high specific energy is prioritized over extreme power capability.

- **Parameterization in this repository**
  - Geometric, transport, and kinetic parameters are defined in `param/params_LGM50.m`. Those parameters are based on [J. Electrochem. Soc. 167, 110559 (2020)](https://iopscience.iop.org/article/10.1149/1945-7111/ab9050) and [Phys. Chem. Chem. Phys. 24, 8661–8677 (2022)](https://pubs.rsc.org/en/content/articlelanding/2022/cp/d2cp00417h).
  - Initial states are generated from `init_cs` and the parameter set in `param/` (including OCV functions), so that DFN, SPMe, and ECM simulations reflect the voltage, SOC, and dynamic behavior of the LG M50 cell under typical operating conditions.


# References

[1] C.-H. Chen, F. B. Planella, K. O’Regan, D. Gastol, W. D. Widanage, and E. Kendrick, “Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models,” *J. Electrochem. Soc.*, 167(8):080534, 2020. doi:10.1149/1945-7111/ab9050.  
[2] S. E. J. O’Kane, W. Ai, G. Madabattula, D. Alonso-Alvarez, R. Timms, V. Sulzer, J. S. Edge, B. Wu, G. J. Offer, and M. Marinescu, “Lithium-ion battery degradation: how to model it,” *Phys. Chem. Chem. Phys.*, 24:7909–7922, 2022. doi:10.1039/D2CP00417H.  
[3] S. Moura, H. Perez, Z. Gima, S. Park, and D. Zhang, “Open software for electrochemical battery modeling, estimation, and control,” *ECS Meeting Abstracts*, 2018, associated with the `fastDFN` code base.  
[4] S. J. Moura, F. B. Argomedo, R. Klein, A. Mirtabatabaei, and M. Krstic, “Battery State Estimation for a Single Particle Model With Electrolyte Dynamics,” *IEEE Trans. Control Syst. Technol.*, 25(2):453–468, 2017. doi:10.1109/TCST.2016.2571663.
