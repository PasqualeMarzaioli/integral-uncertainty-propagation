# Results — INTEGRAL Uncertainty Propagation

**Author:** Pasquale Marzaioli

Numerical output from a complete run of `integral_uncertainty_propagation.m` (MATLAB R2025b, batch mode). All figures are saved under [`plots/`](plots/).

---

## 1. Run configuration

| Quantity | Value |
|----------|--------|
| Spacecraft | INTEGRAL |
| Reference frame (dynamics) | Earth-centred inertial (ECI), J2000 |
| Gravitational parameter $\mu$ | $398600.4418\,\mathrm{km}^3/\mathrm{s}^2$ |
| Initial epoch | `2025-10-26T01:10:57.769` UTC |
| Number of apsides per type | $5$ (five revolutions) |
| Monte Carlo samples $N_{\mathrm{MC}}$ | $2000$ |
| Random seed | $42$ (Mersenne Twister) |
| Convergence tolerance | $0.05\,\sigma$ (component-wise, relative to final sample standard deviation) |
| UT parameters | $\alpha = 0.1$, $\beta = 2$, $\kappa = 0$ |
| Integrator | `ode113`, relative tolerance $10^{-12}$, absolute tolerance $10^{-13}$ |
| Figure export | PNG, $300\,\mathrm{dpi}$, directory `plots/` |

### Initial mean state

Position in $\mathrm{km}$, velocity in $\mathrm{km}/\mathrm{s}$:

$$
\bar{\mathbf{x}}_0 =
\begin{bmatrix}
7896.37847595 \\
-3152.23575783 \\
149961.25963557 \\
-0.6035374932 \\
0.0159618300 \\
0.0321154658
\end{bmatrix}
$$

### Initial covariance $P_0$

Position block in $\mathrm{km}^2$, velocity block in $(\mathrm{km}/\mathrm{s})^2$; position–velocity cross terms are zero:

$$
P_0 =
\begin{bmatrix}
6.7\times 10^{-3} & -2.1\times 10^{-3} & 7.6\times 10^{-4} & 0 & 0 & 0 \\
-2.1\times 10^{-3} & 9.7\times 10^{-3} & 5.3\times 10^{-4} & 0 & 0 & 0 \\
7.6\times 10^{-4} & 5.3\times 10^{-4} & 2.1\times 10^{-3} & 0 & 0 & 0 \\
0 & 0 & 0 & 4.7\times 10^{-7} & 0 & 0 \\
0 & 0 & 0 & 0 & 3.1\times 10^{-7} & 0 \\
0 & 0 & 0 & 0 & 0 & 1.6\times 10^{-7}
\end{bmatrix}
$$

---

## 2. Detected apsides

Keplerian orbital period derived from the initial two-body energy:

$$
T = 227964.078\,\mathrm{s}
\quad (\approx 63.323\,\mathrm{h}).
$$

### Pericentre epochs and Cartesian states

| # | Epoch (UTC) | $r_x$ [km] | $r_y$ [km] | $r_z$ [km] | $v_x$ [km/s] | $v_y$ [km/s] | $v_z$ [km/s] |
|---|-------------|------------|------------|------------|--------------|--------------|--------------|
| 1 | 2025-10-27T08:50:39.808 | $-584.072060276294$ | $233.161673661399$ | $-11092.1964440523$ | $8.15954199503125$ | $-0.215796406027516$ | $-0.434185959423701$ |
| 2 | 2025-10-30T00:10:03.886 | $-584.072060228324$ | $233.161673659339$ | $-11092.1964440145$ | $8.15954199505057$ | $-0.215796406028542$ | $-0.434185959398452$ |
| 3 | 2025-11-01T15:29:27.964 | $-584.072060150789$ | $233.161673656615$ | $-11092.1964439842$ | $8.15954199505955$ | $-0.215796406029679$ | $-0.434185959352897$ |
| 4 | 2025-11-04T06:48:52.042 | $-584.072060132699$ | $233.161673655352$ | $-11092.1964439451$ | $8.15954199507861$ | $-0.215796406030426$ | $-0.434185959341506$ |
| 5 | 2025-11-06T22:08:16.120 | $-584.072060100554$ | $233.161673653324$ | $-11092.1964438867$ | $8.15954199510238$ | $-0.215796406031761$ | $-0.434185959306767$ |

Pericentre radius magnitude (first event): $\|\mathbf{r}_p\| \approx 11109.4\,\mathrm{km}$.

### Apocentre epochs and Cartesian states

| # | Epoch (UTC) | $r_x$ [km] | $r_y$ [km] | $r_z$ [km] | $v_x$ [km/s] | $v_y$ [km/s] | $v_z$ [km/s] |
|---|-------------|------------|------------|------------|--------------|--------------|--------------|
| 1 | 2025-10-28T16:30:21.847 | $7896.37853221587$ | $-3152.23575930360$ | $149961.259631837$ | $-0.603537493113665$ | $0.0159618299649736$ | $0.0321154674682870$ |
| 2 | 2025-10-31T07:49:45.925 | $7896.37853195445$ | $-3152.23575932070$ | $149961.259633078$ | $-0.603537493107977$ | $0.0159618299648500$ | $0.0321154674666123$ |
| 3 | 2025-11-02T23:09:10.003 | $7896.37853112786$ | $-3152.23575921782$ | $149961.259628982$ | $-0.603537493123303$ | $0.0159618299653023$ | $0.0321154674650249$ |
| 4 | 2025-11-05T14:28:34.081 | $7896.37853095894$ | $-3152.23575924930$ | $149961.259630829$ | $-0.603537493115215$ | $0.0159618299651117$ | $0.0321154674634104$ |
| 5 | 2025-11-08T05:47:58.160 | $7896.37852991737$ | $-3152.23575920206$ | $149961.259629878$ | $-0.603537493117808$ | $0.0159618299652606$ | $0.0321154674594425$ |

Apocentre radius magnitude (first event): $\|\mathbf{r}_a\| \approx 150199.9\,\mathrm{km}$, consistent with the supplied initial state (near apocentre).

Successive pericentre (respectively apocentre) times are spaced by one Keplerian period to within the verification tolerance ($0.1\,\mathrm{s}$).

---

## 3. Uncertainty propagation methods

Three independent representations of the mapped uncertainty were evaluated at every detected apsis:

1. **LinCov** — linearized covariance via the state-transition matrix $\Phi(t)$ integrated with the variational equations.
2. **Unscented Transform (UT)** — $2n+1 = 13$ scaled sigma points propagated nonlinearly and recombined with Julier–Uhlmann weights.
3. **Monte Carlo** — $2000$ independent samples drawn from $\mathcal{N}(\bar{\mathbf{x}}_0, P_0)$ and integrated with the same two-body dynamics.

After propagation, means and covariances were rotated into the local radial–transverse–normal (RTH) frame attached to each nominal apsis state.

---

## 4. Monte Carlo convergence

At the final pericentre and final apocentre, the cumulative sample mean was monitored as a function of sample index. Convergence is declared when, for every subsequent cumulative estimate, the maximum absolute component-wise deviation from the full-sample mean stays within $0.05$ times the corresponding component standard deviation (evaluated jointly over both apsides).

**Result:**

$$
N_{\min} = 834
\quad\text{at tolerance}\quad
0.050\,\sigma.
$$

Thus the Monte Carlo mean estimates stabilize well before the full budget of $2000$ samples; the full ensemble is retained for covariance comparisons and ellipse plots.

Console message:

```text
Monte Carlo means remain within 0.050 sigma after N = 834 samples.
```

---

## 5. Verification

All built-in physical and numerical checks **passed**:

| Check | Criterion | Outcome |
|-------|-----------|---------|
| Pericentre spacing | $\|\Delta t_p - T\| < 0.1\,\mathrm{s}$ for successive events | Passed |
| Apocentre spacing | $\|\Delta t_a - T\| < 0.1\,\mathrm{s}$ for successive events | Passed |
| Specific orbital energy | Relative drift $< 10^{-10}$ vs. initial energy | Passed |
| STM vs. finite differences | Relative Frobenius error of $\Phi$ vs. central differences $< 10^{-5}$ at first target | Passed |
| Covariance symmetry / PSD | LinCov, UT, and Monte Carlo series at every epoch | Passed |

Console message:

```text
Uncertainty-propagation verification passed: apsis timing, energy, symmetry, and PSD checks.
```

---

## 6. Figure inventory

All plots are written to [`plots/`](plots/) at $300\,\mathrm{dpi}$:

| File | Content |
|------|---------|
| [`plots/pericentre_mean_difference.png`](plots/pericentre_mean_difference.png) | UT $-$ LinCov mean difference in RTH at pericentres |
| [`plots/apocentre_mean_difference.png`](plots/apocentre_mean_difference.png) | UT $-$ LinCov mean difference in RTH at apocentres |
| [`plots/pericentre_rth_standard_deviations.png`](plots/pericentre_rth_standard_deviations.png) | LinCov vs. UT $1\sigma$ standard deviations (RTH, pericentre) |
| [`plots/apocentre_rth_standard_deviations.png`](plots/apocentre_rth_standard_deviations.png) | LinCov vs. UT $1\sigma$ standard deviations (RTH, apocentre) |
| [`plots/pericentre_rth_cross_covariances.png`](plots/pericentre_rth_cross_covariances.png) | LinCov vs. UT cross-covariances (RTH, pericentre) |
| [`plots/apocentre_rth_cross_covariances.png`](plots/apocentre_rth_cross_covariances.png) | LinCov vs. UT cross-covariances (RTH, apocentre) |
| [`plots/monte_carlo_mean_convergence.png`](plots/monte_carlo_mean_convergence.png) | Cumulative MC mean convergence; $N_{\min}=834$ marked |
| [`plots/pericentre_monte_carlo_distributions.png`](plots/pericentre_monte_carlo_distributions.png) | MC scatter and $3\sigma$ ellipses in the radial–transverse plane (pericentre) |
| [`plots/apocentre_monte_carlo_distributions.png`](plots/apocentre_monte_carlo_distributions.png) | MC scatter and $3\sigma$ ellipses in the radial–transverse plane (apocentre) |
| [`plots/maximum_covariance_axes.png`](plots/maximum_covariance_axes.png) | Maximum $3\sigma$ position/velocity axes; LinCov/UT differences from Monte Carlo |

Interpretation of these figures is developed in the [README](README.md).

---

## 7. Summary of numerical findings

- The orbit is highly eccentric: pericentre near $\|\mathbf{r}_p\|\approx 1.11\times 10^{4}\,\mathrm{km}$ and apocentre near $\|\mathbf{r}_a\|\approx 1.50\times 10^{5}\,\mathrm{km}$, with period $T\approx 2.28\times 10^{5}\,\mathrm{s}$.
- Two-body apsis timing and energy conservation hold to the requested tolerances over five revolutions.
- Monte Carlo sample means are stable after $N_{\min}=834$ samples at the $0.05\sigma$ criterion; $2000$ samples provide a conservative ensemble for covariance comparison.
- LinCov, UT, and Monte Carlo are compared in the RTH frame at every apsis; diagnostic plots quantify mean bias (UT vs. LinCov), standard deviations, cross-covariances, sample geometry, and maximum principal axes relative to Monte Carlo.

---

*End of results document.*
