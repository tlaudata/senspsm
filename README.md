# Overview

This package is a wrapper to conduct sensitivity to unmeasured confounding in 1:1 propensity score matching. You provide your assumptions on the strength of the potential unobserved confounding and confidence interval of ATT is estimated using bootstrap method.

# Setup

To install the package, start by installing devtools package, if not already installed:

```
install.packages("devtools")
```

Then install senspsm package as follows:

```
library(devtools)
install_github("tlaudata/senspsm",build_vignettes=TRUE)
```

Vignette can be accessed using the following commanc

```
vignette("senspsm")
```