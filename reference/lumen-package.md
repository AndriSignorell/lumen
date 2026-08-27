# lumen: Tests and Distributions for DescToolsX

The **lumen** package provides a comprehensive collection of statistical
hypothesis tests and probability distributions designed to complement
the functionality of DescToolsX.

The package focuses on a consistent user interface, extended
methodological coverage, and seamless integration of inferential
procedures and distribution functions within a unified framework.

## Hypothesis Tests

The package implements a wide range of statistical tests:

**Goodness-of-Fit Tests**

- `andersonDarlingTest`, `cramerVonMisesTest`

- `lillieTest`, `jarqueBeraTest`, `shapiroFranciaTest`

- `pearsonTest`

**Nonparametric Tests**

- `signTest`, `jonckheereTerpstraTest`, `pageTest`

- `mosesTest`, `siegelTukeyTest`, `vanWaerdenTest`

**Post-hoc Procedures**

- `dunnTest`, `conoverTest`, `nemenyiTest`

- `scheffeTest`, `postHoc`

**Parametric Tests**

- `tTestA`, `yuenTTest`, `zTest`, `varTest`

- `hotellingsT2Test`, `leveneTest`

**Contingency Table Tests**

- `barnardTest`, `bhapkarTest`, `breslowDayTest`

- `cochranArmitageTest`, `cochranQTest`

- `mantelTrendTest`, `woolfTest`, `stuartMaxwellTest`

**Time Series Tests**

- `adfTest`, `kpssTest`

- `durbinWatsonTest`, `breuschGodfreyTest`

**Randomness and Independence Tests**

- `runsTest`, `BartelsRankTest`, `vonNeumannTest`

## Probability Distributions

The package provides density (`d*`), distribution (`p*`), quantile
(`q*`) and random generation (`r*`) functions for a range of
distributions, following the standard R conventions.

**Extreme Value Distributions**

- Generalized Extreme Value: `dgev`, `pgev`, `qgev`, `rgev`

- Generalized Pareto: `dgpd`, `pgpd`, `qgpd`, `rgpd`

- Gumbel and extended Gumbel: `dgumbel`, `dgumbelx`, ...

- Frechet and reverse Weibull: `dfrechet`, `drweibull`, ...

- Maxima/minima distributions: `dextreme`, `pextreme`, ...

**Special Distributions**

- Benford distribution: `dbenford`, `pbenford`, `qbenford`, `rbenford`

- Order and triangular distributions: `dorder`, `dtri`, ...

## Utilities

- `scores` – Score generation for ordinal contingency tables

- `corTest` – Fast correlation testing for matrices

## Design Principles

- Unified interface across hypothesis tests and distributions

- Standard `d/p/q/r` naming for distributions

- Clear classification via `@family` and `@concept`

- Separation of statistical procedures and data transformation utilities
