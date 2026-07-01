# kami-engine-cfd

`kami-engine-cfd` is the Kotoba/CLJC contract layer for CFD solver requests and
results used by Kami Engine and CAE workflows.

The repository no longer contains a Rust LBM runtime. Native, GPU, or remote CFD
backends should live in adapter repositories and conform to the EDN contracts
defined here.

## Contract

| path | role |
|---|---|
| `src/kami_engine/cfd/contract.cljc` | CLJC request/result validators |
| `test/kami_engine/cfd/contract_test.cljc` | contract and fixture conformance tests |
| `docs/adr/0001-lbm-solver.md` | historical solver design and migration notes |

## EDN Shapes

Requests describe solver intent:

```edn
{:kami.cfd/solver :lbm
 :kami.cfd/dim 3
 :kami.cfd/geometry :ahmed
 :kami.cfd/reynolds 3000.0
 :kami.cfd/steps 1500}
```

Results describe adapter output:

```edn
{:kami.cfd/solver :lbm
 :kami.cfd/dim 3
 :kami.cfd/geometry :ahmed
 :kami.cfd/vehicle-cd 0.29
 :kami.cfd/steps 1500
 :kami.cfd/status :ok}
```

## Development

```bash
clojure -M:test
bb test:cljc
```

The default repository path should not contain `Cargo.toml`, `Cargo.lock`, `.rs`,
or Rust toolchain files.
