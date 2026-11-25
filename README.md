# MATLAB code: kinetic SEIRV model with hierarchical opinion dynamics

This repository contains the MATLAB scripts used for the numerical experiments in

> M. C. Romero Longar, N. Saintier, A. Silva. *Coupling hierarchical opinion formation with epidemic spread: a kinetic SEIRV framework for vaccine hesitancy.*

The model couples an SEIRV epidemic with an opinion dynamics in which agents have a continuous opinion in [-1, 1] and a fixed hierarchy level. Agents adopt the opinion of higher-ranked agents, and a lower-ranked agent influences a higher-ranked one only with probability `p`. The transmission and vaccination rates depend on the opinion.

## Files

`seirv_hierarchy_regular.m` integrates the closed (monokinetic) system for `N` hierarchy levels and reproduces the four scenarios of the paper. The state is `(rho_J^i, M_J^i = rho_J^i m_J^i)`, so no quotients appear and no threshold is needed when a compartment is empty. For each scenario it prints the effective reproduction number at `t = 0` and `t = 50`, the epidemic peak and its time, the recovered and vaccinated fractions at `t = 400`, the consensus opinion `m_inf` and the consensus weights `pi` of Theorem 4.6, and the eigenvalues of `(eta/tau) L`. It plots the compartments, the mean opinion by level, the infectious fraction by level and `R_eff(t)` for the scenario chosen with `plot_case`.

`seirv_hierarchy_extras.m` produces the figures of the sensitivity analysis: heat maps of the recovered fraction and the epidemic peak in the planes `(p, m_S^3(0))` and `(p, tau)`, and a comparison between the closed system and a Monte Carlo simulation of the Boltzmann model with `10^5` agents. The three blocks can be switched off with the flags `do_heat1`, `do_heat2` and `do_mc`.

## Usage

Put both files in the same folder, open MATLAB there and run

```matlab
seirv_hierarchy_regular
seirv_hierarchy_extras
```

No toolboxes are needed. The scenarios are defined in the structure `sc` at the top of `seirv_hierarchy_regular.m`; each one has the fields `bT` (`beta_T`), `bV` (`beta_V`), `p`, `f` (level proportions) and `mS0` (initial opinion of the susceptibles at each level). Common parameters are set in `set_parameters`: `sigma_E = 0.2`, `gamma_I = 0.1`, `eta = 0.05`, `tau = 1` (rates in day^-1, time in days), levels `h = (0, 0.5, 1)`, and initial fractions `rho_S = 0.98`, `rho_I = 0.02` at every level.

## Expected output

For the four scenarios the first script should print the values of Table 1 of the paper: epidemic peaks of 0.0200, 0.0200, 0.0628 and 0.0200, with the interior peak of scenario 3 at `t* = 97.2`, and `rho_R(400)` equal to 0.0494, 0.0207, 0.5989 and 0.0394. The second script writes `heat_p_m3.jpg`, `heat_p_tau.jpg` and `mc_vs_closure.jpg` in the current folder.

## Notes

The closed system is integrated with `ode45` (tolerances `1e-10` and `1e-12`), which is enough for `tau = 1`. For small `tau` the system is stiff, and the second script uses `ode15s` in the `(p, tau)` map. The Monte Carlo block is the slowest part (three initial-opinion modes and four seeds, `10^5` agents, time step 0.1 day); reduce `Na` or `nseed` for a quick test. Monte Carlo results depend on the random seed, so they agree with the paper only up to sampling error.

The equations were also implemented independently in Python (SciPy, Radau method); that version reproduces the values of Table 1 and was used to cross-check the MATLAB code.

Tested with MATLAB R20XX. <!-- fill in the version you used -->

## Citation

If you use this code, please cite the paper above. <!-- add the DOI or arXiv link when available -->
