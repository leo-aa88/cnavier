## Summary

<!-- What changes and why. -->

## Checks

CI runs formatting, static analysis, the CPU tests, the command-line and
regression tests, and the memory checks. It has no GPU, so:

- [ ] If this touches the numerics (`fluiddyn.c`, `poisson.c`, `finitediff.c`,
      `cudasolver.cu`, ...), `make CUDA=1 test` passes on a machine with an
      NVIDIA GPU (exit status 0, no GPU tests skipped)
- [ ] Same for `make OPENMP=1 CUDA=1 test` if it touches the Gauss-Seidel/SOR solvers
- [ ] If results change on purpose, `tests/reference/` is regenerated and the
      reason is given here
