# Identify-sub-network-whose-distribution-altered
This R code tries to find the sub-network whose distribution alters under different conditions.

Suppose we have lots of variables, some of which may alter their distributions when some observable condition is changed. The most commonly used metric to reflex the change in distribution is the expectation, and F-test could be used to detect that. However, when the number of variables is large, F-test will yield a high FDR. If we have another piece of information, say the dependence structure of theses variables, we could make a better identification.

In the example illustrated by main.R, the network represents the potential dependence structure of variables, which are represented by vertices. The red and yellow vertices are variables whose joint distribution is altered the condition is changed. Variables identified by this R code is marked by green.

To run the R code, simply run main.R.
