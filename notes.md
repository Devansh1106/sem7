## MUMPS.jl

- `solve()` function assume that A and rhs are available at all nodes.
    - It does all combined initialize / analyze / factorize / solve.  

- `solve!()` requires problem definition on just host node. Then `factorize!()` and `solve!()` should be called by all the process.
    - `associate!(rhs)` and `associate(A)` can be called on host process only.
