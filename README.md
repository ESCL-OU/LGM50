# LGM50
DFN | SPMe | ECM Simulation of LGM50

This repository is based on and extends the DFN implementation in [fastDFN](https://github.com/scott-moura/fastDFN) and the SPMe implementation in [SPMeT](https://github.com/scott-moura/SPMeT), with only minor modifications to those original codes. We gratefully acknowledge and thank **Dr. Scott Moura** and collaborators for making these high-quality reference implementations publicly available.

## License

This repository incorporates code from the SPMeT project, which is licensed under the GNU General Public License version 3 (GPL‑3.0). In accordance with the GPL‑3.0’s copyleft requirements, **this project is distributed under the terms of the GNU General Public License v3.0**, as provided in the `LICENSE` file. By using, modifying, or redistributing this code, you agree to the terms of the GPL‑3.0.


## Modifications relative to fastDFN and SPMeT

- **DFN implementation (based on `fastDFN` [3])**
  - **Initialization of concentrations**: the original `fastDFN` code initializes electrode and electrolyte concentrations using `init_cs(p, v)`. In this repository, concentrations are instead initialized from:
    - experimentally derived concentration levels corresponding to specific SOC values, as reported in [1] and [2], or
    - steady-state DFN simulations at the desired SOC level.
  - **Thermal dynamics switch**: a parameter flag (e.g., `disable_thermal`) has been added so that the thermal dynamics in the DFN model from [3] can be turned on or off, enabling both isothermal and electrochemical–thermal simulations with a single code base.

- **SPMe implementation (based on `SPMeT` [4])**
  - **Thermal dynamics removed**: only the electrochemical SPMe model from the original SPMeT repository is used; its thermal dynamics are neglected here.
  - **SEI layer growth neglected**: SEI growth dynamics present in the SPMeT model of [4] are ignored in this implementation.
  - **CasADi-based class formulation**: the SPMe model has been refactored into a CasADi-based MATLAB class `SPMe/SPMe.m` that encapsulates parameter loading, state initialization, single-step integration, and voltage/output functions that were originally distributed across multiple scripts in [4].
  - **Initialization of concentrations**: similar to the DFN case, initial solid and electrolyte concentrations are no longer obtained via `init_cs(p, v)`, but are instead set from experimentally based concentration levels derived from [1] and [2] or from simulated profiles at the corresponding SOC level.

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


# To Run
Simply open and run the simulation from MATLAB from the choices of `sim_<model>_<current profile>.m` in home directory


# Using the DFN battery model

The Doyle–Fuller–Newman (DFN) model is implemented in the `DFN` folder and exercised by the example scripts `sim_DFN_CC.m` (constant current) and `sim_DFN_CCCV.m` (CC–CV charging). Refer to [DFN2DAEs.pdf](https://github.com/scott-moura/fastDFN/blob/master/Resources/DFN2DAEs.pdf) from Dr. Moura's repository for numerical implementation. 

- **Parameters and initial conditions**
  - **`load_DFN_params.m`**: builds the parameter struct `p` from `param/params_LGM50.m`, chooses discretization, and constructs all matrices needed by the DFN model.
  - It also creates:
    - **`p.x0`**: initial differential state vector.
    - **`p.z0`**: initial algebraic state vector.
    - **`p.y0`**: initial output vector `[Volt; SOC; …]`.
  - Typical usage:

```matlab
addpath("DFN"); addpath("param");

init_SOC       = 0.3;      % initial SOC (fraction)
t_amb          = 25;       % ambient temperature [degC]
delta_t        = 1.0;      % simulation time step [s]
load_init      = true;     % use pre-saved DFN initial condition
disable_thermal = true;    % ignore temperature dynamics if desired
disable_sei     = true;    % SEI disabled in this implementation

p = load_DFN_params(init_SOC, t_amb, load_init, delta_t, ...
                    disable_thermal, disable_sei);

x = p.x0;   % differential states at current step
z = p.z0;   % algebraic states at current step
y = p.y0;   % outputs at current step
```

### DFN state vectors `x`, `z`, and output `y`

The DFN is formulated as a semi‑explicit DAE:
\[
\dot{x} = f(x,z,I,p), \quad 0 = g(x,z,I,p), \quad y = h(x,z,I,p)
\]

- **Differential states `x`** (vector length `Ncsn + Ncsp + Nce + 1`)
  - `Ncsn = p.PadeOrder * (p.Nxn-1)` : solid-phase concentration modes in the negative electrode.
  - `Ncsp = p.PadeOrder * (p.Nxp-1)` : solid-phase concentration modes in the positive electrode.
  - `Nce = p.Nx - 3` : electrolyte concentration nodes across the cell (excluding boundary points).
  - Layout in `x`:
    - **`c_s_n`** (`Ncsn × 1`): Pade states for solid concentration in the negative electrode.
    - **`c_s_p`** (`Ncsp × 1`): Pade states for solid concentration in the positive electrode.
    - **`c_e`** (`Nce × 1`): electrolyte concentration along `x`.
    - **`T`** (`1 × 1`): cell temperature (if `disable_thermal == false`; otherwise held fixed).

- **Algebraic states `z`** (vector length `3*(Nn+Np) + Nx + 2`)
  - `Nn = p.Nxn-1`, `Np = p.Nxp-1`, `Nx = p.Nx-3`.
  - Layout in `z`:
    - **`phi_s_n`** (`Nn × 1`): solid potential in negative electrode.
    - **`phi_s_p`** (`Np × 1`): solid potential in positive electrode.
    - **`i_en`** (`Nn × 1`): electrolyte current in negative electrode.
    - **`i_ep`** (`Np × 1`): electrolyte current in positive electrode.
    - **`phi_e`** (`(Nx+2) × 1`): electrolyte potential across the full cell including interface points.
    - **`jn`** (`Nn × 1`): interfacial molar flux (negative electrode).
    - **`jp`** (`Np × 1`): interfacial molar flux (positive electrode).

- **Outputs `y`** (stacked diagnostics, one column per time step)
  - For a time history, the scripts store `y` as a matrix with size `N_y × N_T`. The helper `parse_y.m` decodes this stack:
    - **`Volt`** (`1 × N_T`): terminal voltage.
    - **`SOC`** (`1 × N_T`): cell SOC.
    - **`eta_pl`** (`Nn × N_T`): overpotential at the negative current collector.
    - **`c_ss_n`**, **`c_ss_p`**: surface stoichiometries in negative/positive particles.
    - **`c_avg_n`**, **`c_avg_p`**: volume-averaged solid concentrations in negative/positive particles.
    - **`c_ex`** (`(Nx+4) × N_T`): electrolyte concentration including boundary values.
    - **`eta_n`**, **`eta_p`**: Butler–Volmer overpotentials.
    - **`Unref`**: open-circuit potential of the negative electrode at the surface.

The helper functions `parse_x.m`, `parse_z.m`, and `parse_y.m` reshape the stored histories into physically meaningful 2‑D arrays for post‑processing.

### Monitoring lithium plating risk (negative electrode plating potential)

Lithium plating is associated with the local potential of the negative electrode dropping below the Li/Li\(^+\) reference, which is captured in the DFN by the **plating overpotential** `eta_pl` (returned in `y` and unpacked by `parse_y.m`).

- **Extracting plating overpotential near the separator**
  - After a DFN simulation, you can obtain `eta_pl` via:

```matlab
[Volt, SOC, eta_pl, c_ss_n, c_ss_p, c_avg_n, c_avg_p, c_ex, eta_n, eta_p, Unref] = parse_y(y, t_last, p);
```

  - `eta_pl` has size `Nn × t_last`, where each row corresponds to a node in the negative electrode from current collector to separator.
  - The node **closest to the separator** is the last row:

```matlab
eta_pl_sep = eta_pl(end, :);   % plating overpotential at negative electrode / separator interface
```

- **Interpreting `eta_pl_sep`**
  - When **`eta_pl_sep < 0`**, the local negative electrode potential near the separator is below the plating potential, indicating conditions favorable to lithium plating.
  - More negative values correspond to **higher plating driving force**, which is strongly correlated with **accelerated degradation** and, in extreme cases, **internal short circuit risk**.
  - You can monitor this signal during charging to design safer fast‑charge profiles or to enforce constraints in MPC / RL controllers (e.g., keep `eta_pl_sep` above a safety threshold).

### Time stepping with `step_dfn.m`

The low‑level time integration is handled by:

- **`step_dfn.m`**:
  - Signature: `[x_next, z_next, y_next, stats] = step_dfn(x, z, Cur_vec, p, verbosity)`.
  - `Cur_vec = [I_{k-1}, I_k, I_{k+1}]` (in A) is a 3‑point stencil for the applied current used in the Crank–Nicolson time discretization.
  - Internally:
    - Updates solid concentration matrices for the current temperature `T = x(end)`.
    - Calls `cn_dfn.m` to perform a single implicit Crank–Nicolson step on the DAE.
    - Calls `dae_dfn.m` again at the new state to compute the output vector `y_next`.

- **`cn_dfn.m`**:
  - Wraps Newton’s method for the nonlinear implicit step.
  - Solves for consistent algebraic states if the current changes between steps.
  - Uses `dae_dfn.m` and `jac_dfn.m` repeatedly until convergence, and returns iteration statistics in `stats`.

### Using DFN with a custom current‑generating agent

If you want to drive the DFN with a current that is decided at each time step by a custom agent (e.g. an MPC controller or reinforcement learning policy), you can use `step_dfn.m` directly in your loop.

- **Basic loop structure**

```matlab
% 1) Set up parameters and initial state
p = load_DFN_params(init_SOC, t_amb, load_init, delta_t, ...
                    disable_thermal, disable_sei);
x = p.x0;
z = p.z0;
y = p.y0;

Tsim = 3600;                         % total simulation horizon [s]
Nt   = floor(Tsim / p.delta_t);      % number of time steps

I_hist = zeros(1, Nt+1);
I_hist(1) = 0;                       % initial current guess

for k = 1:Nt
    % 2) Agent observes current state/output and decides next current
    Volt = y(1);                     % terminal voltage at time k
    SOC  = y(2);                     % SOC at time k
    obs  = [Volt; SOC];              % (extend as needed)

    I_k   = agent_policy(obs);       % your custom control law [A]
    I_hist(k+1) = I_k;

    % 3) Build Cur_vec = [I_{k-1}, I_k, I_{k+1}]
    if k == 1
        I_prev = I_hist(1);
    else
        I_prev = I_hist(k);
    end

    % If you do not know I_{k+1} yet, you can approximate
    I_next = I_k;                    % zero‑order hold

    Cur_vec = [I_prev, I_k, I_next];

    % 4) Advance DFN one time step
    verbosity = 1;
    [x, z, y, stats] = step_dfn(x, z, Cur_vec, p, verbosity);

    % 5) Check termination conditions (voltage/SOC limits, etc.)
    Volt = y(1);
    SOC  = y(2);
    if Volt < p.volt_min || Volt > p.volt_max
        break;
    end
end
```

- **Notes for custom control**
  - If your agent can compute a full planned sequence of currents, you can populate `Cur_vec` exactly as in `sim_DFN_CC.m`: `Cur_vec = [I(k-1), I(k), I(k+1)]`.
  - If the agent is strictly causal (decides only the current at the present step), using `Cur_vec = [I_{k-1}, I_k, I_k]` is a reasonable zero‑order approximation.
  - You can expose richer observations to the agent (e.g. `SOC`, `c_ss_n`, `c_ss_p`, etc.) by unpacking `y` via `parse_y.m`.

### DFN folder contents and roles

The `DFN` folder contains the core numerical implementation of the DFN model:

- **`step_dfn.m`**: single‑step time advancement wrapper (builds solid concentration matrices at current `T`, calls `cn_dfn`, and computes `y`).
- **`cn_dfn.m`**: Crank–Nicolson nonlinear solver for the DAE using Newton iterations.
- **`dae_dfn.m`**: defines the DFN differential and algebraic equations (`f`, `g`) and diagnostic outputs (`y` and optional detailed outputs).
- **`jac_dfn_pre.m`**: builds state‑independent parts of the Jacobian matrices used in Newton’s method.
- **`jac_dfn.m`**: completes the Jacobian by adding state‑dependent terms (e.g. Butler–Volmer kinetics, temperature coupling, electrolyte conductivity).
- **`load_DFN_params.m`**: assembles the parameter struct `p`, discretization sizes, precomputes matrices (`c_s_mats`, `c_e_mats`, `phi_s_mats`, `i_e_mats`, `phi_e_mats`), initializes `x0`, `z0`, and `y0`.
- **`parse_x.m`, `parse_y.m`, `parse_z.m`**: utilities to unpack the stacked state/output histories into structured physical fields (concentrations, potentials, fluxes, etc.).
- **`potentials.m`**: convenience wrapper around `dae_dfn` to extract detailed potential‑related quantities (overpotentials, reference potentials, IR drops, electrolyte contributions).
- **`c_s_mats.m`**: builds state‑space matrices for solid diffusion (Pade approximation) in both electrodes.
- **`i_e_mats.m`**: builds matrices for electrolyte current dynamics in each region.
- **`phi_s_mats.m`**: assembles matrices for solid‑phase potential in electrodes.
- **`phi_e_mats.m`**: assembles matrices for electrolyte potential (including boundary conditions).
- **`blkdiagFast.m`**: efficient block‑diagonal concatenation helper used when assembling Jacobians and diffusion operators.
- **`symsubsden.m`, `symsubsnum.m`**: helper routines used for symbolic substitutions when generating or manipulating CasADi symbolic expressions (not typically called directly in simulations).

The overall simulation flow is:

1. **Setup**: `sim_DFN_CC.m` or `sim_DFN_CCCV.m` calls `load_DFN_params.m` (and adds `DFN`/`param` to the MATLAB path).
2. **Precomputation**: `load_DFN_params.m` calls the various `*_mats.m` builders and `jac_dfn_pre.m`, and evaluates `dae_dfn.m` once at `(x0, z0)` to create `y0`.
3. **Time stepping**: the sim script builds a current profile `I(k)` and, at each step, calls `step_dfn.m` with the appropriate `Cur_vec`.
4. **Solving DAEs**: `step_dfn.m` invokes `cn_dfn.m`, which repeatedly calls `dae_dfn.m` and `jac_dfn.m` until convergence.
5. **Post‑processing**: after the loop, the sim script uses `parse_x.m`, `parse_z.m`, and `parse_y.m` to unpack the trajectories and saves them into the `results` directory.

With this structure, you can either use the provided `sim_DFN_CC.m` / `sim_DFN_CCCV.m` as templates or directly integrate `step_dfn.m` into your own closed‑loop or agent‑driven simulations.


# Using the SPMe model

The Single Particle Model with electrolyte (SPMe) provides a **reduced‑order approximation** of the DFN model, retaining key physics (solid diffusion in each electrode, 1‑D electrolyte dynamics, voltage and SOC) while being significantly faster and easier to use in control and estimation.

- **Key differences from DFN**
  - Each electrode is represented by a **single representative particle** (rather than many particles along the electrode thickness).
  - The electrolyte is still modeled along the through‑thickness direction as like DFN.
  - The model is implemented in a compact CasADi‑based class `SPMe/SPMe.m` that exposes a discrete‑time **step interface**.
  - Plating overpotential is summarized by a **scalar** `eta_pl` output (see `SPMe.m`), as opposed to the spatially distributed `eta_pl` in the DFN.

- **Basic usage**
  - The script `sim_SPMe_CC.m` shows how to configure and run the SPMe:

```matlab
import casadi.*
addpath("SPMe"); addpath("param");

% Load parameters and initial condition (similar arguments to DFN)
Nr    = 30;
Nxn   = 70; Nxs = 35; Nxp = 70;
delta_t = 1.0;
p = load_SPMe_params(init_SOC, t_amb, load_init, Nr, delta_t, ...
                     Nxn, Nxs, Nxp, disable_thermal, disable_sei);

% Build SPMe system object and initialize state
spme_ss = SPMe(p, 'idas');
x = p.x0;

for k = 1:Nt-1
    I_k = I(k);  % applied current at time step k
    [x, y, c_ss_n, c_ss_p, c_ex] = spme_ss.step(x, I_k);
    V   = y(1);  % terminal voltage
    SOC = y(2);  % SOC
    eta_pl = y(3);  % negative electrode plating overpotential (scalar)
end
```

For more detailed SPMe equations and derivations, refer to Scott Moura’s original SPMeT repository [`SPMeT`](https://github.com/scott-moura/SPMeT).


# Using the ECM (Equivalent Circuit Model)

The ECM is a **lumped, low‑order model** that represents the battery as an open‑circuit voltage source plus one or more RC branches, parameterized as a function of SOC and temperature.
It is intended for applications where **computational efficiency** is critical (e.g., real‑time BMS) and detailed electrochemical states are not required.

- **Model structure**
  - Implemented in `ECM/ECM.m` with parameters loaded by `load_ECM_params.m`.
  - Supports different orders (number of RC branches) via the `n_order` argument in `sim_ECM_CC.m`.
  - Outputs are primarily **terminal voltage** and **SOC**, with internal RC states capturing dynamic polarization.

- **Basic usage**
  - The script `sim_ECM_CC.m` shows a typical workflow:

```matlab
import casadi.*
addpath("ECM"); addpath("param");

n_order = 2;             % number of RC branches
CRate   = -1.0;          % charging/discharging C‑rate
delta_t = 1.0;

% Load parameters and ECM system object
[p, ecm_p, ocv_p] = load_ECM_params(init_SOC, n_order);
ecm_ss = ECM(p, ocv_p, ecm_p, n_order, 'idas');

% Build current profile I(k) (see sim_ECM_CC.m)
for k = 1:Nt-1
    I_k = I(k);
    [x, V] = ecm_ss.step(x, I_k);   % x: [SOC; RC states; ...], V: terminal voltage
end
```


# LG M50 cell description

The models in this repository are parameterized for **LG’s M50** lithium‑ion cell, a high‑energy‑density **21700‑format** cylindrical cell (approx. 21 mm diameter, 70 mm length) with a nominal capacity of about **5 Ah** and nominal voltage around **3.6–3.7 V**.

- **Chemistry and use case**
  - The cell employs a **graphite negative electrode** and a layered oxide **nickel‑rich positive electrode** (NMC‑type), representative of modern high‑energy automotive cells.
  - It is designed for **energy‑focused applications** (e.g., EV packs, stationary storage) where high specific energy is prioritized over extreme power capability.

- **Parameterization in this repository**
  - Geometric, transport, and kinetic parameters are defined in `param/params_LGM50.m`. Those parameters are based on [J. Electrochem. Soc. 167, 110559 (2020)](https://iopscience.iop.org/article/10.1149/1945-7111/ab9050) and [Phys. Chem. Chem. Phys. 24, 8661–8677 (2022)](https://pubs.rsc.org/en/content/articlelanding/2022/cp/d2cp00417h).
  - Initial condition files in `init_models/` and data in `param/` (e.g., OCV curves) are calibrated so that DFN, SPMe, and ECM simulations reflect the voltage, SOC, and dynamic behavior of the LG M50 cell under typical operating conditions.


# References

[1] C.-H. Chen, F. B. Planella, K. O’Regan, D. Gastol, W. D. Widanage, and E. Kendrick, “Development of Experimental Techniques for Parameterization of Multi-scale Lithium-ion Battery Models,” *J. Electrochem. Soc.*, 167(8):080534, 2020. doi:10.1149/1945-7111/ab9050.  
[2] S. E. J. O’Kane, W. Ai, G. Madabattula, D. Alonso-Alvarez, R. Timms, V. Sulzer, J. S. Edge, B. Wu, G. J. Offer, and M. Marinescu, “Lithium-ion battery degradation: how to model it,” *Phys. Chem. Chem. Phys.*, 24:7909–7922, 2022. doi:10.1039/D2CP00417H.  
[3] S. Moura, H. Perez, Z. Gima, S. Park, and D. Zhang, “Open software for electrochemical battery modeling, estimation, and control,” *ECS Meeting Abstracts*, 2018, associated with the `fastDFN` code base.  
[4] S. J. Moura, F. B. Argomedo, R. Klein, A. Mirtabatabaei, and M. Krstic, “Battery State Estimation for a Single Particle Model With Electrolyte Dynamics,” *IEEE Trans. Control Syst. Technol.*, 25(2):453–468, 2017. doi:10.1109/TCST.2016.2571663.
