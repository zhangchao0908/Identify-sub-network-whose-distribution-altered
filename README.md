# Identify-sub-network-whose-distribution-altered
This R code tries to find the sub-network whose distribution alters under different conditions.

Suppose we have lots of variables, some of which may alter their distributions when some observable condition is changed. The most commonly used metric to reflex the change is the expectation. If the condition takes two values, say 0 and 1, F-test could be used to detect the change in the expectation. However, when the number of variables is large, F-test will yield a high FDR. If we have another piece of information, say the dependence structure of theses variables, we could make better identification.
