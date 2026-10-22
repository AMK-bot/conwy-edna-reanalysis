# Smoke-test receipts

Auditable PASS lines written by `R/smoke_test_*.R` (latest run per engine).

- **Stan** (cmdstanr 0.9.0 / CmdStan 2.39.0): PASS ✓ — Bernoulli compile+sample, theta=0.749, Rhat 1.014, 10.5s
- **JAGS** (occumb 1.3.0 / jagsUI 1.6.3 / rjags 4.17 / JAGS 4.3.2): PASS ✓ — occumb chain runs; jagsUI occupancy converges (psi_hat=0.66, max Rhat 1.001), 1.4s
