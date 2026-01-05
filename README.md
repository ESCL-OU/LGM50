# LGM50
DFN | SPMe | ECM Simulation of LGM50

# Install Requirements

We make use of [CasADi](https://web.casadi.org/) for the symbolic variables to be integrated with MPC and optimal control methods. 
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


