# Results — INTEGRAL Uncertainty Propagation

**Author:** Pasquale Marzaioli

Numerical output from a complete run of `integral_uncertainty_propagation.m` (MATLAB R2025b, batch mode). All figures are saved under [`plots/`](plots/). The state tables retain solver digits for deterministic comparison; those digits do not imply equivalent physical orbit accuracy.

---

## 1. Run configuration

| Quantity | Value |
|----------|--------|
| Scenario | INTEGRAL-inspired benchmark; not an operational ephemeris |
| Reference frame (dynamics) | Earth-centred inertial (ECI); source orientation not documented |
| Gravitational parameter $\mu$ | $398600.4418\,\mathrm{km}^3/\mathrm{s}^2$ |
| Initial epoch | `2025-10-26T01:10:57.769` UTC |
| Number of apsides per type | $5$ (five revolutions) |
| Monte Carlo samples $N_{\mathrm{MC}}$ | $2000$ |
| Random seed | $42$ (Mersenne Twister) |
| Mean-stability tolerance | $0.05\,\sigma$ (component-wise, relative to the same run's final sample mean and standard deviation) |
| UT parameters | $\alpha = 0.1$, $\beta = 2$, $\kappa = 0$ |
| Integrator | `ode113`, relative tolerance $10^{-12}$, absolute tolerance $10^{-13}$ |
| Figure export | PNG, $300\,\mathrm{dpi}$, directory `plots/` |

### Initial mean state

The repository does not document the estimation source of the state or covariance. These values define a reproducible illustrative benchmark and must not be interpreted as a traceable flight orbit-determination solution.

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

The same state gives $a=80656.052\,\mathrm{km}$ and $e=0.8622545$, hence analytic two-body radii $r_p=11110.010\,\mathrm{km}$ and $r_a=150202.093\,\mathrm{km}$.

### Pericentre epochs and Cartesian states

| # | Epoch (UTC) | $r_x$ [km] | $r_y$ [km] | $r_z$ [km] | $v_x$ [km/s] | $v_y$ [km/s] | $v_z$ [km/s] |
|---|-------------|------------|------------|------------|--------------|--------------|--------------|
| 1 | 2025-10-27T08:50:39.808 | $-584.072060276294$ | $233.161673661399$ | $-11092.1964440523$ | $8.15954199503125$ | $-0.215796406027516$ | $-0.434185959423701$ |
| 2 | 2025-10-30T00:10:03.886 | $-584.072060228324$ | $233.161673659339$ | $-11092.1964440145$ | $8.15954199505057$ | $-0.215796406028542$ | $-0.434185959398452$ |
| 3 | 2025-11-01T15:29:27.964 | $-584.072060150789$ | $233.161673656615$ | $-11092.1964439842$ | $8.15954199505955$ | $-0.215796406029679$ | $-0.434185959352897$ |
| 4 | 2025-11-04T06:48:52.042 | $-584.072060132699$ | $233.161673655352$ | $-11092.1964439451$ | $8.15954199507861$ | $-0.215796406030426$ | $-0.434185959341506$ |
| 5 | 2025-11-06T22:08:16.120 | $-584.072060100554$ | $233.161673653324$ | $-11092.1964438867$ | $8.15954199510238$ | $-0.215796406031761$ | $-0.434185959306767$ |

Pericentre radius magnitude (first event): $\|\mathbf{r}_p\| = 11110.010\,\mathrm{km}$.

### Apocentre epochs and Cartesian states

| # | Epoch (UTC) | $r_x$ [km] | $r_y$ [km] | $r_z$ [km] | $v_x$ [km/s] | $v_y$ [km/s] | $v_z$ [km/s] |
|---|-------------|------------|------------|------------|--------------|--------------|--------------|
| 1 | 2025-10-28T16:30:21.847 | $7896.37853221587$ | $-3152.23575930360$ | $149961.259631837$ | $-0.603537493113665$ | $0.0159618299649736$ | $0.0321154674682870$ |
| 2 | 2025-10-31T07:49:45.925 | $7896.37853195445$ | $-3152.23575932070$ | $149961.259633078$ | $-0.603537493107977$ | $0.0159618299648500$ | $0.0321154674666123$ |
| 3 | 2025-11-02T23:09:10.003 | $7896.37853112786$ | $-3152.23575921782$ | $149961.259628982$ | $-0.603537493123303$ | $0.0159618299653023$ | $0.0321154674650249$ |
| 4 | 2025-11-05T14:28:34.081 | $7896.37853095894$ | $-3152.23575924930$ | $149961.259630829$ | $-0.603537493115215$ | $0.0159618299651117$ | $0.0321154674634104$ |
| 5 | 2025-11-08T05:47:58.160 | $7896.37852991737$ | $-3152.23575920206$ | $149961.259629878$ | $-0.603537493117808$ | $0.0159618299652606$ | $0.0321154674594425$ |

Apocentre radius magnitude (first event): $\|\mathbf{r}_a\| = 150202.093\,\mathrm{km}$, consistent with the supplied initial state.

Successive pericentre (respectively apocentre) times are spaced by one Keplerian period to within the verification tolerance ($0.1\,\mathrm{s}$).

---

## 3. Uncertainty propagation methods

Three independent representations of the mapped uncertainty were evaluated at every detected apsis:

1. **LinCov** — linearized covariance via the state-transition matrix $\Phi(t)$ integrated with the variational equations.
2. **Unscented Transform (UT)** — $2n+1 = 13$ scaled sigma points propagated nonlinearly and recombined with Julier–Uhlmann weights.
3. **Monte Carlo** — $2000$ independent samples drawn from $\mathcal{N}(\bar{\mathbf{x}}_0, P_0)$ and integrated with the same two-body dynamics.

After propagation, means and covariances were rotated into the local radial–transverse–normal (RTH) frame attached to each nominal apsis state. All methods are evaluated at the fixed nominal event times; sigma points and Monte Carlo samples do not use sample-specific apsis events.

---

## 4. Monte Carlo mean stability

At the final pericentre and final apocentre, the cumulative sample mean was monitored as a function of sample index. The reported index is the first point after which every subsequent cumulative estimate remains within $0.05$ times the corresponding component standard deviation of the same run's full-sample mean, evaluated jointly over both apsides.

**Result:**

$$
N_{\min} = 834
\quad\text{at tolerance}\quad
0.050\,\sigma.
$$

Thus the cumulative means are retrospectively stable well before the full budget of $2000$ samples. This is not an independent convergence proof because the final sample defines the reference and the condition necessarily holds at $N=2000$. The full ensemble is retained as a finite-sample covariance reference; for Gaussian data, a sample variance at this size has approximate relative standard error $\sqrt{2/(N-1)}=3.16\%$.

Console message:

```text
Relative to the full sample, Monte Carlo means remain within 0.050 sigma after N = 834 samples.
```

---

## 5. Quantitative method comparison

For each epoch, the largest position and velocity $3\sigma$ axes are $3\sqrt{\lambda_{\max}}$ of the corresponding $3\times3$ covariance block. Endpoint values and signed discrepancies from the Monte Carlo estimate are:

| Epoch | MC position axis [km] | LinCov difference | UT difference | MC velocity axis [km/s] | LinCov difference | UT difference |
|-------|-----------------------|-------------------|---------------|---------------------------|-------------------|---------------|
| Pericentre 1 | $922.127$ | $+0.314\%$ | $+0.314\%$ | $0.363320$ | $+0.372\%$ | $+0.371\%$ |
| Pericentre 5 | $6187.547$ | $+2.105\%$ | $+2.086\%$ | $2.387676$ | $+4.534\%$ | $+4.465\%$ |
| Apocentre 1 | $102.872$ | $+1.029\%$ | $+1.029\%$ | $0.00370241$ | $+1.149\%$ | $+1.149\%$ |
| Apocentre 5 | $514.340$ | $+1.032\%$ | $+1.033\%$ | $0.0151297$ | $+1.120\%$ | $+1.120\%$ |

The dominant axis can hide errors in smaller RTH components. Final-epoch component standard deviations and relative discrepancies are:

| RTH component | Pericentre 5 MC | LinCov difference | UT difference | Apocentre 5 MC | LinCov difference | UT difference |
|---------------|-----------------|-------------------|---------------|----------------|-------------------|---------------|
| $r_r$ [km] | $148.273$ | $-81.770\%$ | $+4.956\%$ | $1.03209$ | $-95.478\%$ | $-0.236\%$ |
| $r_t$ [km] | $2062.413$ | $+2.102\%$ | $+2.083\%$ | $171.447$ | $+1.032\%$ | $+1.033\%$ |
| $r_h$ [km] | $1.96485$ | $-99.631\%$ | $-97.915\%$ | $0.189279$ | $-48.233\%$ | $-48.207\%$ |
| $v_r$ [km/s] | $0.795777$ | $+4.541\%$ | $+4.472\%$ | $0.00499781$ | $+1.121\%$ | $+1.121\%$ |
| $v_t$ [km/s] | $0.104039$ | $-89.784\%$ | $+8.828\%$ | $0.000677449$ | $+1.092\%$ | $+1.094\%$ |
| $v_h$ [km/s] | $0.00741044$ | $+1.588\%$ | $+1.588\%$ | $0.000553572$ | $+0.589\%$ | $+0.589\%$ |

The maximum component-wise mean discrepancy, standardized by the Monte Carlo standard error $\sqrt{P_{\mathrm{MC},jj}/N}$, is $28.96$ for LinCov and $2.25$ for UT over all epochs. This supports the nonlinear UT mean correction; it does not make the Monte Carlo mean exact.

Sampling noise remains material for covariance comparisons. Splitting the ensemble into two independent halves of $1000$ changes the final maximum position axis by $1.31\%$ at pericentre and $1.25\%$ at apocentre. Differences of similar size should not be ranked as method accuracy.

The propagated distribution also departs from Gaussian shape at pericentre. At the fifth pericentre, the largest absolute component-wise sample skewness is $2.84$ and the largest absolute excess kurtosis is $12.43$; the corresponding fifth-apocentre values are $0.515$ and $0.518$. Covariance ellipses therefore summarize scale and correlation but do not fully describe the late-pericentre distribution.

---

## 6. Verification

All built-in physical and numerical checks **passed**:

| Check | Criterion | Outcome |
|-------|-----------|---------|
| Pericentre spacing | $\|\Delta t_p - T\| < 0.1\,\mathrm{s}$ for successive events | Passed |
| Apocentre spacing | $\|\Delta t_a - T\| < 0.1\,\mathrm{s}$ for successive events | Passed |
| Analytic apsis radii | Numerical radii within $10^{-5}\,\mathrm{km}$ of $a(1\mp e)$ | Passed |
| Specific orbital energy | Relative drift $< 10^{-10}$ vs. initial energy | Passed |
| Specific angular momentum | Relative vector drift $< 10^{-10}$ vs. initial value | Passed |
| STM vs. finite differences | Relative Frobenius error of $\Phi$ vs. central differences $< 10^{-5}$ at first target | Passed |
| UT input moments | Dimensionless mean/covariance reconstruction errors below $10^{-12}$ / $10^{-8}$ | Passed |
| Monte Carlo states | All propagated states finite | Passed |
| Covariance integrity | Finite, symmetric, and PSD after correlation normalization at every epoch | Passed |

Console message:

```text
Verification passed: max apsis-spacing error 8e-06 s, relative energy drift
2.34e-11, relative angular-momentum drift 1.02e-11, STM error 5.11e-09.
```

---

## 7. Figure inventory

All plots are written to [`plots/`](plots/) at $300\,\mathrm{dpi}$:

| File | Content |
|------|---------|
| [`plots/pericentre_mean_difference.png`](plots/pericentre_mean_difference.png) | UT $-$ LinCov mean difference in RTH at pericentres |
| [`plots/apocentre_mean_difference.png`](plots/apocentre_mean_difference.png) | UT $-$ LinCov mean difference in RTH at apocentres |
| [`plots/pericentre_rth_standard_deviations.png`](plots/pericentre_rth_standard_deviations.png) | LinCov, UT, and MC $1\sigma$ standard deviations (RTH, pericentre) |
| [`plots/apocentre_rth_standard_deviations.png`](plots/apocentre_rth_standard_deviations.png) | LinCov, UT, and MC $1\sigma$ standard deviations (RTH, apocentre) |
| [`plots/pericentre_rth_cross_covariances.png`](plots/pericentre_rth_cross_covariances.png) | LinCov, UT, and MC cross-covariances (RTH, pericentre) |
| [`plots/apocentre_rth_cross_covariances.png`](plots/apocentre_rth_cross_covariances.png) | LinCov, UT, and MC cross-covariances (RTH, apocentre) |
| [`plots/monte_carlo_mean_convergence.png`](plots/monte_carlo_mean_convergence.png) | Retrospective cumulative MC mean stability; $N_{\min}=834$ marked |
| [`plots/pericentre_monte_carlo_distributions.png`](plots/pericentre_monte_carlo_distributions.png) | MC scatter and $3\sigma$ ellipses in the radial–transverse plane (pericentre) |
| [`plots/apocentre_monte_carlo_distributions.png`](plots/apocentre_monte_carlo_distributions.png) | MC scatter and $3\sigma$ ellipses in the radial–transverse plane (apocentre) |
| [`plots/maximum_covariance_axes.png`](plots/maximum_covariance_axes.png) | Maximum $3\sigma$ position/velocity axes; LinCov/UT differences from Monte Carlo |

Interpretation of these figures is developed in the [README](README.md).

---

## 8. Summary of numerical findings

- The orbit is highly eccentric: pericentre near $\|\mathbf{r}_p\|\approx 1.11\times 10^{4}\,\mathrm{km}$ and apocentre near $\|\mathbf{r}_a\|\approx 1.50\times 10^{5}\,\mathrm{km}$, with period $T\approx 2.28\times 10^{5}\,\mathrm{s}$.
- Two-body apsis timing and energy conservation hold to the requested tolerances over five revolutions.
- Monte Carlo cumulative means remain within the retrospective $0.05\sigma$ band after $N_{\min}=834$ samples; this is a stability diagnostic, not proof of population-moment convergence.
- UT materially improves the propagated mean relative to LinCov; its worst standardized mean discrepancy from Monte Carlo is $2.25$, compared with $28.96$ for LinCov.
- Maximum covariance-axis discrepancies are $0.3\%$–$4.5\%$, while Monte Carlo covariance sampling noise is itself material at $N=2000$.
- At pericentre 5, UT captures radial-position and transverse-velocity dispersion far better than LinCov, but both deterministic methods miss at least $97.9\%$ of the Monte Carlo normal-position standard deviation.
- Strong late-pericentre skewness and excess kurtosis show that a Gaussian covariance ellipse is not a complete uncertainty description there.
- The two-body-only result isolates Keplerian nonlinearity and is not suitable for operational prediction of the actual spacecraft.

---

*End of results document.*
