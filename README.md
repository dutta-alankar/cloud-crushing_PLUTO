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
  - `boost = True` Sets whether to turn Gallilean frame boost on/off to track the cloud and save computational expense (see [Dutta+19](https://ui.adsabs.harvard.edu/abs/2019RNAAS...3..148D/abstract) for the algorithm).
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

![cloud-mass](https://github.com/dutta-alankar/cloud-crushing_PLUTO/assets/39578361/9353ac16-ad7d-4eec-921c-57b3ea6ef38c)

### Theoretical background

Cloud-crushing simulations are self-similar for different physical properties of the cloud. This means the cloud evolution is identical for different values of these cloud properties when appropriately scaled. Nevertheless, differences can arise in evolution even for degenerate properties, especially in late times, which are purely numerical. Determining this degeneracy across parameters is not only useful in saving computation while making a parameter search of cloud survival/destruction but also in providing a guideline on the robust choice of `code units`. These `code units`, by construction, generate identical initial conditions for different values of degenerate parameters. It is important to note that the set of degenerate parameters is also dependent on the physics involved in the simulations.

We solve the following equations in the presence of optically radiative cooling on PLUTO.

$\partial_t\rho + \nabla . (\rho \textbf{v}) = 0,$

$\partial_t(\rho \textbf{v}) + \nabla . (\rho \textbf{v} \otimes \textbf{v} + p \mathcal{I}) = 0,$

$\partial_t(\frac{1}{2}\rho v^2 + e) + \nabla . [(\frac{1}{2}\rho v^2 + e + p )\textbf{v}] = -n_H^2 \Lambda(T)$.

The symbols here have their usual meaning, and this set of equations is closed by the ideal gas equation of state $e = p/(\gamma -1)$.
We can now de-dimensionalize these equations in units of cloud size $R_{\rm cl}$, cloud density $\rho_{\rm cl}$ and velocity of the wind $v_{\rm wind}$, which can be conveniently expressed in terms of the wind Mach number $\mathcal{M}$. The choice of units that de-dimensionalize the hydrodynamic equations are as follows,

$v_0 = v_{\rm wind} = \mathcal{M} \sqrt{\gamma \eta \frac{k_B T_{\rm cl}}{\mu m_p}}$

$L_0 = R_{\rm cl} = v_{\rm wind} t_{\rm cool, mix}/\left( \sqrt{\chi} \left(t_{\rm cool,mix}/t_{\rm cc}\right) \right),$ 

$\rho_0 = \rho_{\rm wind} = \left(n_{\rm cl}/\chi \right) \mu m_p$. 

> [!NOTE] 
> Apart from the free parameters described earlier, the cloud number density $n_{\rm cl}$ is an additional, degenerate parameter involved in the de-dimensionalization.

Therefore, the de-dimensionalized equations (de-dimenionalized variables are represented by $\sim$ on their top) are as follows,

$\partial_{\rm \tilde{t}}\tilde{\rho} + \tilde{\nabla} . (\tilde{\rho} \tilde{\textbf{v}}) = 0,$

$\partial_{\rm \tilde{t}}(\tilde{\rho} \tilde{\textbf{v}}) + \tilde{\nabla} . (\tilde{\rho} \tilde{\textbf{v}} \otimes \tilde{\textbf{v}} + \tilde{p} \mathcal{I}) = 0,$

$\partial_{\rm \tilde{t}}(\frac{1}{2}\tilde{\rho} \tilde{v}^2 + \tilde{e}) + \tilde{\nabla} . [(\frac{1}{2}\tilde{\rho} \tilde{v}^2 + \tilde{e} + \tilde{p} )\tilde{\textbf{v}}] = - \zeta \tilde{n}_H^2 \tilde{\Lambda}(\tilde{T)}$.

Here $\zeta = \frac{1}{\gamma (\gamma -1)}\mathcal{M}^{-2} \chi^{-3/2} \left(\chi/\eta\right)^{3/2} \left(\frac{t_{\rm cool, mix}}{t_{\rm cc}}\right)^{-1}\left[\frac{\Lambda(T_0)}{\Lambda(\sqrt{\eta} T_{\rm cl})}\right]$ and $T_0 = \mathcal{M}^2 \gamma \eta T_{\rm cl}$. Non-radiative simulations correspond to $\zeta =0$.

With $L_0$, $v_0$ and $\rho_0$ (as presented earlier) as the choice of code units, the initial pressures $\tilde{p}$ in the cloud and the wind in code units are 

$p_{\rm wind}^{\rm (code)} = \frac{n_{\rm wind} k_B T_{\rm wind}}{\rho_0 v_0^2}=\frac{1}{\gamma \mathcal{M}^2}$ and

$p_{\rm cl}^{\rm (code)} = \frac{n_{\rm cl} k_B T_{\rm cl}}{\rho_0 v_0^2}=\frac{(\chi/\eta)}{\gamma \mathcal{M}^2}$ respectively.

By design, $v_{\rm wind}^{\rm (code)} = 1$ and $R_{\rm cl}^{\rm (code)} = 1$.

One can clearly see that any quantity that doesn't appear in either the equations or in the initialization is degenerate. For example, $n_{\rm cl}$ is a degenerate quantity (so are absolute pressure values in the simulation). 

The following cartoon demonstrates the initialization in our code.

![image](https://github.com/dutta-alankar/cloud-crushing_PLUTO/assets/39578361/e11ca732-af1e-4699-bbc4-5f0f1257d2f5)
