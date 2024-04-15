# linicrypt-solver

Describe your project here.

## Results
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
