# Changelog

## Unreleased

### Added
- Regression coverage for scalar uniform spline values and derivatives across degrees 1–5,
  shifted grids, knot boundaries, the first grid interval and extrapolation policies.

### Fixed
- Correct interval selection for scalar uniform spline values, including multi-interval
  linear and clamped cubic splines.
- Prevent index underflow when evaluating scalar cubic derivatives in the first grid interval.
- Correct coefficient indexing for scalar highest-order derivatives, including the first
  derivative of a linear spline.

## Rusty-Fitpack 0.1.3 (2026-09-23)

### Added
- `Extrapolation` enum (`Clamp`, `Zero`, `SmoothDecay(width)`) and the functions
  `splev_uniform_ext`, `splder_uniform_ext` and `spline_support` to choose how a spline with
  uniform knots is continued outside its support. Values and derivatives are always consistent.
  `SmoothDecay` continues the spline beyond its last knot with a quintic polynomial that matches
  value, first and second derivative and decays to zero (C2) within the given width, e.g. for
  tabulated integrals that should vanish beyond the end of their grid.

### Changed

### Fixed

### Documentation
- Documented that the legacy `splev_uniform`/`splder_uniform` are not consistent outside the
  support (clamped value, but one-sided derivative at the nearest end). Their behaviour is
  unchanged for backward compatibility; `splev_uniform_ext`/`splder_uniform_ext` provide
  consistent continuations.


### Internal

## Rusty-Fitpack 0.1.2 (2023-11-18)
[2a1ad8f...HEAD](https://github.com/mitric-lab/Rusty-FITPACK/compare/2a1ad8f...HEAD)

### Added

### Changed

### Fixed
Fixed all warnings

### Documentation

### Internal

## Rusty-Fitpack 0.1.1 (2023-11-17)
[6790a512...HEAD](https://github.com/mitric-lab/Rusty-FITPACK/compare/6790a512...HEAD)

### Added
- Added metadata to Cargo.toml
- First registering/publishing to crates.io

### Changed
- made fpbspl public [#1](https://github.com/mitric-lab/Rusty-FITPACK/pull/1)

### Fixed
Fixed all warnings

### Documentation

### Internal

