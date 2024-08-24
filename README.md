# linicrypt-solver

## Setup the enviroment
I use rye for managing the dependencies: [install rye](https://rye.astral.sh/).

Then set up the venv with this command in the root directory of this project
```bash
  rye sync
```
If you didn't enable the rye shim to shadow the local python installation you
need to activate the venv with `source ./.venv/bin/activate` in bash.

## Results from the galois implementation of the CR corollary

The implementation in `src/linicrypt_solver` consists mainly of these modules:
- `field.py`: Definition of the field that is used. I am using the prime order field $F_29$
- `__init__.py`: Defines the interface/trait `Constraint`
- `random_oracle.py`: Constraint is implemented for the Random Oracle constraints
- `ideal_cipher.py`: Constraint is implemented for Ideal Cipher constraints
- `solvable.py`: Implements an API for a set of constraints,
  in particular the algorithms for finding out if it is solvable fixing a subspace of the dual space,
  outside of a subspace of the state space.
- `algebraic_representation.py`: Here the CR and 2PR corollaries are implemented with the high level functions provided by lower modules.
- `merkle_damgard.py`: Implements helper methods to construct the PGV compression functions and their MD constructions.
- `cli.py`: Examples are implemented here

The most important method on a `Constraint` is `map`.
This enables:
- Forcing variables to collapse
- Checking CR for repeated nonces szenarios by solving the set of constraints in which the outputs are collapsed
- Easy construction of Merkle-Damgard construction by collapsing the input and output variables of the chained compression functions

The implementation of CR and 2PR reproduces the categorization from PGV and BRS.
That means every secure function is identified as secure,
and for every insecure function an attack is found in form of a
set of solvable linicrypt constraints and a linear function mapping its solutions to collisions.
The corollary goes further: this code finds every possible "fast" attack on CR ("slow" meaning guessing with negligible probability)

Setting the log level to `INFO` produces this output when testing CR for all PGV compression functions in $H_2$.
Example output:
```
E(1h + 1m, 0h + 0m) + 1h + 0m PGV: ('B', 59), BRS: g
2024-08-24 15:04:31 | INFO
H_n:
I:
[[1 0 0 0 0 0 0]
 [0 1 0 0 0 0 0]
 [0 0 0 1 0 0 0]
 [0 0 0 0 0 1 0]]
O:
[[1 0 0 0 0 0 0]
 [0 0 0 0 0 0 1]]
C:
[0 0 0 0 0 0 0] <- [1 1 0 0 0 0 0] -> [28  0  1  0  0  0  0]
[0 0 0 0 0 0 0] <- [0 0 1 1 0 0 0] -> [ 0  0 28  0  1  0  0]
[0 0 0 0 0 0 0] <- [0 0 0 0 1 1 0] -> [ 0  0  0  0 28  0  1]
  0%|                                                                  | 0/203 [00:00<?, ?it/s]2024-08-24 15:04:31 | INFO
Found solvable subspace:
2024-08-24 15:04:31 | INFO
Partition of the constraints is [[2, 3], [0, 1, 4, 5]]
2024-08-24 15:04:31 | INFO
Solvable subspace of F^(2d) is
2024-08-24 15:04:31 | INFO
[[ 1  0  0  0  0]
 [ 0  1  0  0  0]
 [ 0  0  1  0  0]
 [ 1  1 28  0  0]
 [28  0  2  0  0]
 [ 0  0  0  1  0]
 [ 0  0  0  0  1]
 [ 1  0  0  0  0]
 [27  0  2  1  0]
 [ 2  0 27  0  1]
 [28  1  2  0 28]
 [ 1  0 28  0  1]
 [ 0  1  1  0 28]
 [ 0  0  0  0  1]]
2024-08-24 15:04:31 | INFO
Solvable constraints in that subspace are
2024-08-24 15:04:31 | INFO
[0 0 0 0 0] <- [1 1 0 0 0] -> [28  0  1  0  0]
[0 0 0 0 0] <- [28  0  2  1  0] -> [ 1  0 27  0  1]
  5%|███                                                     | 11/203 [00:00<00:01, 139.97it/s]
Collision resistant: False

...

E(0h + 1m, 0h + 1m) + 0h + 1m PGV: ('-', 5), BRS: a
  0%|                                                                  | 0/203 [00:00<?, ?it/s]2024-08-24 15:04:50 | INFO
Found solvable subspace:
2024-08-24 15:04:50 | INFO
Partition of the constraints is [[0], [1, 2, 3, 4, 5]]
2024-08-24 15:04:50 | INFO
Solvable subspace of F^(2d) is
2024-08-24 15:04:50 | INFO
[[1 0 0 0 0]
 [0 1 0 0 0]
 [0 0 1 0 0]
 [0 0 0 1 0]
 [0 0 0 0 1]
 [0 0 0 1 0]
 [0 0 0 0 1]
 [1 0 0 0 0]
 [0 0 0 1 0]
 [0 0 0 0 1]
 [0 0 0 1 0]
 [0 0 0 0 1]
 [0 0 0 1 0]
 [0 0 0 0 1]]
2024-08-24 15:04:50 | INFO
Solvable constraints in that subspace are
2024-08-24 15:04:50 | INFO
[0 1 0 0 0] <- [0 1 0 0 0] -> [ 0 28  1  0  0]
[0 0 0 1 0] <- [0 0 0 1 0] -> [ 0  0  0 28  1]
  0%|▎                                                         | 1/203 [00:00<00:02, 81.24it/s]
Collision resistant: False
E(0h + 1m, 1h + 0m) + 0h + 1m	PGV: ('B', 21), BRS: c
100%|████████████████████████████████████████████████████████| 203/203 [00:02<00:00, 69.30it/s]
Collision resistant: True
```
It also confirms the suspicion of the authors of 'Characterizing Collision and Second-Preimage Resistance in Linicrypt'.
In the section "5.2 Why the Restriction to Distinct Nonces?" they look at the example $P(x,y) = H(H(x)) - H(y)$.

This this example is not collision resistant, but the authors claim that it might by second preimage resistant.
The two corollaries and this code confirm this claim. The output of `cli.py:test_cr()` is

```
The program
I:
[[1 0 0 0 0]
 [0 1 0 0 0]]
O:
[[0 0 0 1 1]]
C:
[1 0 0 0 0] |-> [0 0 1 0 0]
[0 0 1 0 0] |-> [0 0 0 1 0]
[0 1 0 0 0] |-> [0 0 0 0 1]

 15%|████████▌                                               | 31/203 [00:00<00:00, 306.10it/s]2024-08-24 16:18:01 | INFO
Found solvable subspace:
2024-08-24 16:18:01 | INFO
Partition of the constraints is [[0], [1, 2, 3], [4, 5]]
2024-08-24 16:18:01 | INFO
Solvable subspace of F^(2d) is
2024-08-24 16:18:01 | INFO
[[1 0 0 0]
 [0 1 0 0]
 [0 1 0 0]
 [0 0 1 0]
 [0 0 1 0]
 [0 1 0 0]
 [0 0 1 0]
 [0 0 1 0]
 [0 0 0 1]
 [0 0 0 1]]
2024-08-24 16:18:01 | INFO
Solvable constraints in that subspace are
2024-08-24 16:18:01 | INFO
[1 0 0 0] |-> [0 1 0 0]
[0 1 0 0] |-> [0 0 1 0]
[0 0 1 0] |-> [0 0 0 1]
 17%|█████████▋                                              | 35/203 [00:00<00:00, 255.55it/s]
is CR: False
100%|███████████████████████████████████████████████████████| 203/203 [00:01<00:00, 142.87it/s]
is 2PR: True
```

## Results from the sympy implementation for the first conjecture

This only finds basis change matrices, collapses are not implemented yet.

```
  Permutation attack c1 <-> c2:
-----------------------------
f(0h + 0m, 0h + 0m) + 0h + 1m	PGV: ('-', 56), BRS: a
⎡ 1    0    0      0        0   ⎤
⎢                               ⎥
⎢B₁₀  B₁₁  B₁₂  B₁₃ + 1  B₁₄ - 1⎥
⎢                               ⎥
⎢B₁₀  B₁₁  B₁₂    B₁₃      B₁₄  ⎥
⎢                               ⎥
⎢ 0    1   -1      0        1   ⎥
⎢                               ⎥
⎣ 0    0    0      0        1   ⎦

f(0h + 0m, 0h + 0m) + 1h + 0m	PGV: ('-', 60), BRS: a
⎡ 1    0    0    0    0 ⎤
⎢                       ⎥
⎢B₅   B₆   B₇   B₈   B₉ ⎥
⎢                       ⎥
⎢ 1    0   -1    0    1 ⎥
⎢                       ⎥
⎢B₁₅  B₁₆  B₁₇  B₁₈  B₁₉⎥
⎢                       ⎥
⎣ 0    0    0    0    1 ⎦

f(0h + 0m, 0h + 1m) + 1h + 0m	PGV: ('D', 12), BRS: b
⎡1  0  0   0  0⎤
⎢              ⎥
⎢0  0  0   1  0⎥
⎢              ⎥
⎢1  0  -1  0  1⎥
⎢              ⎥
⎢0  1  0   0  0⎥
⎢              ⎥
⎣0  0  0   0  1⎦

f(0h + 0m, 1h + 1m) + 1h + 0m	PGV: ('D', 44), BRS: b
⎡1   0  0   0  0 ⎤
⎢                ⎥
⎢-1  0  1   1  0 ⎥
⎢                ⎥
⎢1   0  -1  0  1 ⎥
⎢                ⎥
⎢0   1  1   0  -1⎥
⎢                ⎥
⎣0   0  0   0  1 ⎦

f(0h + 1m, 0h + 0m) + 1h + 0m	PGV: ('P', 57), BRS: f
⎡1  0  0   0  0⎤
⎢              ⎥
⎢0  0  0   1  0⎥
⎢              ⎥
⎢1  0  -1  0  1⎥
⎢              ⎥
⎢0  1  0   0  0⎥
⎢              ⎥
⎣0  0  0   0  1⎦

f(0h + 1m, 0h + 1m) + 1h + 0m	PGV: ('P', 9), BRS: f
⎡1  0  0   0  0⎤
⎢              ⎥
⎢0  0  0   1  0⎥
⎢              ⎥
⎢1  0  -1  0  1⎥
⎢              ⎥
⎢0  1  0   0  0⎥
⎢              ⎥
⎣0  0  0   0  1⎦

f(1h + 1m, 0h + 0m) + 1h + 0m	PGV: ('B', 59), BRS: g
⎡1   0  0   0  0 ⎤
⎢                ⎥
⎢-1  0  1   1  0 ⎥
⎢                ⎥
⎢1   0  -1  0  1 ⎥
⎢                ⎥
⎢0   1  1   0  -1⎥
⎢                ⎥
⎣0   0  0   0  1 ⎦

f(1h + 1m, 1h + 1m) + 1h + 0m	PGV: ('B', 43), BRS: g
⎡1   0  0   0  0 ⎤
⎢                ⎥
⎢-1  0  1   1  0 ⎥
⎢                ⎥
⎢1   0  -1  0  1 ⎥
⎢                ⎥
⎢0   1  1   0  -1⎥
⎢                ⎥
⎣0   0  0   0  1 ⎦

f(0h + 0m, 0h + 0m) + 1h + 1m	PGV: ('D', 64), BRS: b
⎡   1        0        0         0        0   ⎤
⎢                                            ⎥
⎢ -B₁₅    1 - B₁₆    -B₁₇    1 - B₁₈   -B₁₉  ⎥
⎢                                            ⎥
⎢1 - B₁₅  1 - B₁₆  -B₁₇ - 1   -B₁₈    1 - B₁₉⎥
⎢                                            ⎥
⎢  B₁₅      B₁₆      B₁₇       B₁₈      B₁₉  ⎥
⎢                                            ⎥
⎣   0        0        0         0        1   ⎦

f(0h + 0m, 0h + 1m) + 1h + 1m	PGV: ('P', 16), BRS: f
⎡1  0  0   0  0⎤
⎢              ⎥
⎢0  0  0   1  0⎥
⎢              ⎥
⎢1  0  -1  0  1⎥
⎢              ⎥
⎢0  1  0   0  0⎥
⎢              ⎥
⎣0  0  0   0  1⎦

f(0h + 1m, 0h + 0m) + 1h + 1m	PGV: ('P', 61), BRS: f
⎡1  0  0   0  0⎤
⎢              ⎥
⎢0  0  0   1  0⎥
⎢              ⎥
⎢1  0  -1  0  1⎥
⎢              ⎥
⎢0  1  0   0  0⎥
⎢              ⎥
⎣0  0  0   0  1⎦

f(0h + 1m, 0h + 1m) + 1h + 1m	PGV: ('P', 13), BRS: f
⎡1  0  0   0  0⎤
⎢              ⎥
⎢0  0  0   1  0⎥
⎢              ⎥
⎢1  0  -1  0  1⎥
⎢              ⎥
⎢0  1  0   0  0⎥
⎢              ⎥
⎣0  0  0   0  1⎦

12
```
