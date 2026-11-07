# canonical root の Casimir 補足：raw contraction の全線形必要条件

2026-11-07。主義務 G。新しい手証明の producer note、独立 review 前。
凍結済 `CUBIC_SCROLL_NORMAL_ROOT_DIRECTRIX_REDUCTION.md` は変更しない。

## 1. 扱う全源と未処理の補集合

smooth cubic scroll `S⊂P4`、実既約被約 quintic `D`、一般 rank 二 restriction の仮定を使う。
前 note の残る rank 二 root pair 五組のうち、`A0=C0+4f` を含む四組、及び rank 一 root `A0` を扱う。
この saturated root は、`q=q2=x2x4−x3²` の normal derivative の generic line に等しい。

**必要条件。** `Pi` を全150 homogeneous coefficients の canonical divergence-zero lift とすると

`Pi#(dq)|S = 0`                                           (1)

が raw polynomial vector の全五成分について成立する。
projective tangency の `Pi#dq∈C·Euler` より強い条件である。
これを残る `A=(2,5)` rank 一及び `(2,5)+(2,5)` rank 二へ拡張しない。
全 G の上界又は空性をこの補題だけから主張しない。

## 2. canonical lift の strict Casimir

`E=Σxi∂i` とし、`Pi` の coefficient degree は二、`div Pi=0` とする。
projective Jacobi は `E∧[Pi,Pi]=0` である。
BV divergence は `div[Pi,Pi]=0` を与え、`[E,[Pi,Pi]]=0`。
従って `div(E∧[Pi,Pi])=5[Pi,Pi]`（全体の符号 convention は零判定に影響しない）となり、
characteristic 零から strict `[Pi,Pi]=0`。

標準 volume に `Pi∧Pi` を contraction した polynomial one-form を `omega` とする。
`div(Pi∧Pi)=2[Pi,Pi]=0` により `d omega=0`。
`omega` の係数は degree 四なので Cartan formula は

`d(i_E omega)=5 omega`。

従って実 homogeneous quintic `D=i_E omega` は、固定非零定数を除いて Euler-contracted actual Pfaffian である。
rank 四 bivector の linear algebra から `Pi#omega=0` なので

`Pi#dD=0`。                                              (2)

`D` の正規化定数は以下の零恒等式へ影響しない。

## 3. projective root-line preserving から raw 零へ

`X={q,−}_Pi` とする。parallel normal root の必要条件は、chart `x0≠0` で

`{q/x0²,xi/x0} ∈ (q/x0²)+I_S²`。

分母を払う numerator は正確に

`H_i=x0{q,xi}−xi{q,x0}−2q{x0,xi}`。                     (3)

`dH_i∧dq|S=0` と projective tangency はこの必要条件を polynomial に表す。
`Q={q=0}` の generic smooth point of `S` では、`S` は Cartier divisor である。
`A0` root は generic `dq` line なので、normal quadratic の `z²` 係数が零となり、
`D|Q` の `S` に沿う order は `m≥3`。
`D` は既約かつ quadric `q` の倍でないため `D|Q≠0`、従って `m` は有限である。

(3) を `Q` に制限すると、ある regular rational function `Lambda=X(x0)/x0` に対して

`X = Lambda E + Y`, `Y∈I_(S/Q)²·T_Q`                    (4)

が `S` の generic point 近傍で成立する。
ここでは `Q` の singular locus `C0` を除いており、全 polynomial 恒等式は最後に dense extension する。

`E(D)=5D` と (2) から

`0=X(D)|Q=5 Lambda(D|Q)+Y(D|Q)`。

局所 DVR の uniformizer を `z` とすれば、`D|Q=z^m u`、`u` は unit。
`Y` の全係数は `z²` の倍なので `Y(D|Q)∈(z^(m+1))`。
従って `Lambda∈(z)` である。(4) より `X|S=0` が従う。
全成分は polynomial なので generic `S\C0` 上の零が integral `S` 全体へ延び、(1) を得る。

この議論では canonical lift を勝手に Euler gauge へ変更していない。
非 canonical lift では (2) がそのまま成り立たないので、(1) を全 raw lift に課してはいけない。
また `D|Q` を零へ置いて有限 order を仮定してもいない。

## 4. 最小の次の有限目標

凍結済 full21 canonical source へ、次の三組の線形必要条件を一括して戻す。

1. 全源の `R⊃C0`、すなわち `t18=t10=0`。
2. 全 intrinsic root-line equations `(dH_i∧dq)|S=0`。
3. 本補題の全五 raw equations `(Pi#dq)|S=0`。

complete linear kernel の affine dimension が七以下なら、この root branch の projective PGL 像は
`24+dim(kernel)−6−1≤24` となる。
dimension が大きければ、その完全 kernel を残し、実 `D` の因子又は whole nonlinear caseproof が必要である。
generic minor の増設、GB/jet の反復又は他の二 root 型の成功への読み替えはしない。

先行する root-line 条件だけの全620×150 matrix は source 上 rank 五/kernel十六であり、
`CLASSIFICATION_CUBIC_SCROLL_RANK_TWO_CANONICAL_ROOTLINE_20261107.json` に未閉結果として凍結した。
この新補題はその minor の取り替えではなく、strict Casimir の finite-order argument による別の必要条件である。
