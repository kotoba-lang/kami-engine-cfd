# kami-engine-cfd

The Rust `kami-cfd` prototype has been retired from this repository.

Current ownership:

- `kotoba-lang/cae-solver` owns the portable `cae.solver/solve` dispatch
  contract, including the `:lbm` solver kind.
- `kotoba-lang/aero` owns the reduced-order aerodynamic drag model and case/datom
  interface.
- High-fidelity CFD execution is now a host-adapter concern. A host may provide
  an `:lbm` implementation behind the same EDN case contract, but this repository
  no longer owns Cargo, Rust source, or binary build steps.

`docs/adr/0001-lbm-solver.md` remains as a migration record for the removed
prototype and its calibration notes.
