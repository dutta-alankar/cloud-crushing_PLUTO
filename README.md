# cloud-crushing_PLUTO

## This is a numerical setup in `PLUTO` that implements the cloud-crushing problem with optically thin radiative cooling

> [!IMPORTANT]
> We are using a custom-modified version of `PLUTO` 4.4 patch 2. The original `PLUTO` code has been written and maintained by Andrea Mignone & collaborators.

### Getting started
- To set up the problem, one needs to run the Python script `prob-prepare.py` located in `python-scripts` directory in the repo (tested on Python version 3.11).
- The following are the free parameters in the problem that can be set in the Python script.
  - Cloud to Wind density contrast $\chi$.
  - Wind to Cloud temperature contrast $\eta$.
  - Mach number of the wind $\mathcal{M}$.
  - Cooling time of mixed gas to cloud crushing time $t_{\rm cool, mix}/t_{\rm cc}$, which sets the cloud size in physical units.
  - Cloud temperature $T_{\rm cl}$.
  - Metallicity with respect to Solar $Z/Z_\odot$.
  - The initial position of the center of the spherical cloud in the simulation domain in the wind direction.
  - Adiabatic index of the gas $\gamma$.
- Additionally, to get quantities in physical units (CGS), one needs to set the cloud density $n_{\rm cl}$, which doesn't directly come up in our problem, but sets the code units to CGS conversion factor.
- Other additional flags are also set in the script. The following are some of the most important ones.
  - `cooling = True/False` Sets optically thin radiative cooling on/off in the simulation.
  - `catalyst = False` Sets whether in-situ `Catalyst` visualization support will be enabled.
  - `auto_compile = True` Sets whether to enable/disable automatic compilation & creation of a separate directory for the compiled binary to run.
  - `boost = True` Sets whether to turn Gallilean frame boost on/off to track the cloud and save computational expense.
  - Change `cooltable_name` in the Python script to use any cooling rate table that must be placed in the `cooltables` directory.
- Edit the sample `local_make` to link `hdf5` library \& `Catalyst` in-situ visualization library (if enabled).
- Sample job scripts in `slurm` is provided in the `jobscripts` directory that can be configured according to the cluster environment where the code will be running.
> [!NOTE]  
> The Python script is self-contained and generates everything necessary for these cloud-crushing simulations. There is **no** need to run `setup.py` as is traditionally done in `PLUTO`.
### Additional info
Like many other numerical problems, cloud-crushing simulations behave best when the code units reflect the dimensions involved in the problem. One such choice is as follows.
  - Code length = Cloud radius
  - Code density = Wind density
  - Code velocity = Wind velocity

Using fluid equations rescaled to these units, it can be shown that the gas pressure/temperature is a degenerate quantity in the absence of radiative cooling.
In these code units, the following are the fluid quantities.
  - Radius of cloud = 1.0
  - Density of cloud = $\chi$
  - Velocity of the cloud = 0
  - Pressure in the cloud = $\frac{\chi/\eta}{\gamma \mathcal{M}^2}$
  - Density of the wind = 1.0
  - Velocity of the wind = 1.0
  - Pressure in the wind = $\frac{1}{\gamma \mathcal{M}^2}$

In these code units, the cloud-crushing time has a trivial value $t_{\rm cc} = \sqrt{\chi}$.

As a demonstration of this code, the following is the evolution of cold mass with radiative cooling at different values of $t_{\rm cool, mix}/t_{\rm cc}$. The other parameters are set to ($\chi$, $\eta$, $\mathcal{M}$, $T_{\rm cl}$, $Z/Z_\odot$) = ($100$, $100$, $1.5$, $4 \times 10^{4}$, $1.0$). In these runs, the Sutherland+93 CIE cooling curve has been used.

![cloud-mass](https://github.com/dutta-alankar/cloud-crushing_PLUTO/assets/39578361/43e8642a-39b2-44df-9a58-9d5590ef7a89)
