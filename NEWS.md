# lorps 0.1.0

- Provides `lorps_mi()` for fixed-rank or LPPL-selected matrix imputation.
- Supports five disjoint validation subsets by default and portable PSOCK CV.
- Preserves joint missing-entry draws in a compact representation.
- Adds `complete_lorps()` and `lorps_suggest_K()`.
- Validates MCMC settings against the compiled sampler's exact saving rule.
