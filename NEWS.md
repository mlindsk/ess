# molic v1.1.1.2 (2021-05-31)

 * Small changes in default parameter values

# molic v1.1.1 (2021-05-28)

 * Dependency on package =Matrix= removed
 * The =fit_graph= is now lightyears faster than the previous version (redundant copying was removed)
 * If the number of observations is small, the underlying probability tables tends to be sparse. In these cases the usual penalty term of AIC and BIC is often too restrictive. One can now set =sparse_qic= to =TRUE=. This penality is computed according to a sparse criteria. The criteria resembles the usual penalty as the number of observations grow.
 * The =q= parameter (0 for AIC and 1 for BIC) can now be any positive number. The higher, the more penalty and the fewer edges will be added.

# molic v1.0.0 (2020-05-24)

 * First release.
