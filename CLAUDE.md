# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

EM ("Embedded MESA") is a Fortran library that drives MESA stellar evolution in-process and calls GYRE at each step to compute oscillation frequencies, stopping evolution where model frequencies best match observed ones. It is used as the fitness function for asteroseismic model fitting (the AMP pipeline, a parallel genetic algorithm in `src/amp/`).

## Build

Requires the MESA SDK to be initialized (`$MESASDK_ROOT` set; minimum version enforced by `build/check_sdk_version`), plus `MESA_DIR` and `GYRE_DIR` pointing at builds whose versions exactly match `src/pkgconfig/mesa.pc.in` and `gyre.pc.in` (checked by `scripts/configure_*.sh`). When bumping MESA/GYRE, update the `Version:` in those `.pc.in` files (and `Requires.private` in `src/pkgconfig/em.pc`).

```
make                  # build + install into bin/, lib/, include/, lib/pkgconfig/
make DEBUG=yes        # -Og, -fcheck=all, -Wall, backtraces
make clean
make check_src        # verify every source file listed in build/Makefile exists
make -C build linter  # fortitude lint (config: build/fortitude.toml) on the generated .f90 files in build/
```

Options in the top-level `Makefile` (`FRONTENDS`, `DEBUG`, `SHARED`, `OMP`, `FPE`, `CRMATH`, `IEEE`) are passed as `make VAR=yes|no`.

Build mechanics: all compilation happens in `build/`. `.fypp` sources are found via `vpath` across `src/*`, preprocessed with fypp (from the MESA SDK's Python) into `build/*.f90`, then compiled with gfortran `-std=f2018`. Dependencies are generated automatically (`fypp_deps`, `makedepf08.awk`). **Adding a new source file requires adding it to the relevant `*_SRCS` list in `build/Makefile`** (and a new directory to `SRC_DIRS` in the top-level `Makefile`). `src/amp/` is built separately by its own legacy `Makefile.f77` (MPI), not by the main build.

## Testing

There is no unit test suite. The two frontends serve as end-to-end checks (both run full MESA+GYRE evolutions, so they are slow):

- `bin/em_test` — hard-coded observed frequencies; creates a star, evolves to ZAMS, then runs `evolve_star_seismic` and prints the best model and its frequencies.
- `bin/testffamp2 <params...>` — calls the AMP fitness function `userff` (`src/frontend/ff_amp2.fypp`) with scaled parameters from the command line and prints `chisq_r`. Reads `obs.dat` from the current directory (header line: `n iwflag isflag`; `n` lines of `l nu dnu ...`; then up to 5 non-seismic constraint lines). `isflag == 0` selects the frequency-ratio chi-squared; otherwise surface-corrected frequencies are used.

## Architecture

Module layering (each `.fypp` defines one module of the same name):

- `em_m` (`src/lib`) — public facade. Re-exports selected procedures from the modules below, holds best-model scalars (Teff, logg, R, L, age, ...), and provides a legacy `evolve_star_seismic` overload that stores these via a callback. External callers (AMP) `use em_m` only.
- `em_mesa_m` (`src/mesa`) — wraps MESA `star_lib`: `create_star` sets up a star entirely in code (empty inlist filenames — no inlist files are read), `evolve_star*` drive the evolution loop. `evolve_star_seismic_process_` hooks `s%extras_check_model` with a contained state machine (`START`/`EVOLVE`/`BRACKET`/...) that starts calling GYRE once the model large separation is within `f_enter` of observed, brackets the sign change in the lowest radial-mode frequency offset, refines the timestep to locate the best match, and stops when Δν has dropped `f_exit` past observed.
- `em_gyre_m` (`src/gyre`) — builds a GYRE `evol_model_t` directly from the in-memory MESA star and runs GYRE's mode finder for a given `l`, scanning a frequency window around the observed modes. Returns a `freqs_t`.
- `freqs_m` (`src/freq`) — `freqs_t` type (nu, dnu, n_pg, E_norm; one harmonic degree per instance, arrays of them indexed `0:L_MAX` with `L_MAX = 3`), module-level storage of observed/model frequencies, mode matching and surface corrections (`apply_cubic_correction`, `apply_combined_correction`).
  - Radial orders matter: ratios pair modes by `n_pg`. Observed data must be built with the order-assigning constructor — `freqs_t(nu, dnu, E_norm, 0)` for radial modes first, then `freqs_t(nu, dnu, E_norm, l, fr_0)` for l > 0 — which assigns asymptotic orders (gap-aware, following AMP 1.3). Model `n_pg` are GYRE's own orders and are used as-is; they follow the same convention for pure p modes, but not for mixed modes. `check_mod_orders` (called from `eval_chisq_ratios`) only prints a warning when they disagree with the asymptotic relation.
- `ratios_m` (`src/freq`) — `ratios_t`, constructed as `ratios_t(fr, 'R010'|'R02'|'R13', rt_ref)` from a `freqs_t(0:L_MAX)` array (Roxburgh & Vorontsov 2003, see `doc/`); with `rt_ref` the ratios are interpolated onto its frequencies (used for model ratios). `eval_ratios_deriv` gives finite-difference derivatives for the covariance.
- `chisq_m` (`src/freq`) — chi-squared evaluators: individual frequencies (`eval_chisq_radial`, `eval_chisq_all`) and ratios (`eval_chisq_ratios`, full covariance via LAPACK95). Its `DEBUG` parameter prints ratios in AMP 1.3's output format (`ratio nu_ref rt_obs rt_err rt_mod`) for comparison with AMP.
- `interp_m`, `deriv_m` (`src/math`) — cubic interpolation/derivative helpers (adapted and simplified from GYRE).
- `utils_m` (`src/lib`) — small utilities (e.g. `linear_fit`).

Dependencies come from MESA (`star_lib`, etc.), GYRE (`gyre_m`), and ForUM (`forum_m`: `RD` kind, `ASSERT` macros via `forum.inc`, pulled in by `src/include/em.inc`).

## Conventions

- Every module/program starts with `#:include 'em.inc'` and uses `implicit none (type, external)`; reals are `real(RD)`.
- Private helper procedures carry a trailing underscore (`run_gyre` vs `make_gyre_model_`); module variables carry `_m`.
- Error checks use the fypp macro `@:ASSERT(cond, 'message')`.
- Max line length 132 (fortitude). `ff_amp2` is an external function (not in a module) for compatibility with the F77 AMP driver, hence its lint exemptions.
