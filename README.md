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
  - **SEI layer growth neglected**: SEI growth dynamics present in the SPMeT model of [4] are ignored in this implementation.
  - **CasADi-based class formulation**: the SPMe model has been refactored into a CasADi-based MATLAB class `SPMe/SPMe.m` that encapsulates parameter loading, state initialization, single-step integration, and voltage/output functions that were originally distributed across multiple scripts in [4].

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


# Running simulations

Use this entry point for normal simulations:

- `sim_profile_sequence_LGM50.m`

It supports both single and concatenated profiles for all three models:

- **Models**: `DFN`, `SPMe`, `ECM`
- **Profiles**:
  - dynamic: `BJDST`, `DST`, `FUDS`, `US06`
  - constant current: `CC` (direction determined by `Crate` sign)
    - `Crate > 0`: discharge
    - `Crate < 0`: charge

Initialization is SoC-based for all models:

1. `calculate_equilibrium_voltage_from_soc(init_soc)` computes equilibrium voltage.
2. `init_cs` computes initial concentrations from that voltage.
3. Model states are initialized from those concentrations (with model-specific state construction).

## Quick start

```matlab
% SPMe, single dynamic profile
sim = sim_profile_sequence_LGM50('SPMe', 0.4, 'DST');

% DFN, single CC 1C charge from SOC=0.4
sim = sim_profile_sequence_LGM50('DFN', 0.4, {'CC', -1});

% ECM, single dynamic profile with save
sim = sim_profile_sequence_LGM50('ECM', 0.4, 'FUDS', 'SaveOutput', true);

% SPMe, concatenated 1C charge then 1C discharge (one cycle)
sim = sim_profile_sequence_LGM50('SPMe', 0.0, {'CC', -1, 'CC', 1});

% SPMe, repeat the same charge/discharge sequence for 5 cycles
sim = sim_profile_sequence_LGM50('SPMe', 0.0, {'CC', -1, 'CC', 1}, 'Repeat', 5);

% SPMe, repeated FUDS -> 1C charge with custom voltage window
sim = sim_profile_sequence_LGM50('SPMe', 0.4, {'FUDS', 'CC', -1}, ...
    'Repeat', 3, 'VMin', 2.8, 'VMax', 4.2);
```

## Options (`sim_profile_sequence_LGM50`)

```matlab
sim = sim_profile_sequence_LGM50(model_name, init_soc, profile_sequence, ...
    'Nxn', 70, 'Nxs', 35, 'Nxp', 70, 'Nr', 30, ...
    'DeltaT', 1, 'Tamb', 25, ...
    'ECMOrder', 2, ...
    'Crate', -1, ...
    'Repeat', 1, ...
    'VMin', [], 'VMax', [], ...
    'ShowProgress', true, ...
    'VEnd', [], ...
    'SaveOutput', false, 'OutputFile', '');
```

- **`Nxn`, `Nxs`, `Nxp`, `Nr`**: discretization options (DFN/SPMe).
- **`DeltaT`**: simulation timestep `[s]`.
- **`Tamb`**: ambient temperature `[degC]`.
- **`ECMOrder`**: number of RC branches for ECM.
- **`Crate`**: used only if a token is `'CC'` without explicit rate (advanced compatibility).
- **`Repeat`**: repeat entire `profile_sequence` this many times.
- **`VMin` / `VMax`**: optional voltage window truncation for each segment.
- **`ShowProgress`**: show a waitbar while running segments/repeats.
- **`VEnd`** (optional custom stop voltage):
  - dynamic and CC discharge stop when `V <= VEnd`
  - CC charge stops when `V >= VEnd`
  - dynamic profiles automatically repeat until `VEnd` is reached (subject to internal safety repeat cap).
- **`SaveOutput` / `OutputFile`**: save simulation struct to `.mat`.

`profile_sequence` formats:

- single dynamic profile: `'DST'`, `'FUDS'`, `'BJDST'`, `'US06'`
- single CC profile: `{'CC', -1}` or `{'CC', 1}`
- concatenated sequence: `{'FUDS', 'CC', -1}` or `{'CC', -1, 'CC', 1}`

## Outputs

`sim_profile_sequence_LGM50` returns a struct `sim` containing:

- common fields: `time`, `current`, `voltage`, `soc`, `model`, `profile`
- sequence fields: `profile_definition`, `segment_index`, `segment_cycle`, `segment_log`, `repeat_count`
- initialization diagnostics: `init_soc_input`, `init_voltage_from_soc`, `init_csn0`, `init_csp0`
- termination info: `terminated_early`, `termination_reason`, `t_last`
- model-specific states:
  - DFN: `x`, `z`, `y`, `c_s_n_pade`, `c_s_p_pade`, `c_e`, `temperature`
  - SPMe: `x`, `y`, `c_s_n`, `c_s_p`, `c_e`
  - ECM: `x`, plus ECM parameter snapshots


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
