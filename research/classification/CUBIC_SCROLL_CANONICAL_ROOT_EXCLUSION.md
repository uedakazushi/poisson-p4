# smooth cubic scroll の canonical q2-root 全空：Hartogs と Casimir・divergence の矛盾

2026-11-07。主義務 G。純粋な幾何・多項式証明。追加 CAS は不要。
producer 証明であり、別作者の独立 review を受けるまで全分類の登録結果にはしない。
既存主原稿、旧 report、旧負の dimension 結果を変更しない。

## 1. 正確な仮定と新しい結論

`S=S(1,2)≅F1⊂P4`、`H=C0+2f` とし、ideal を

`q0=x0x3−x1x2`, `q1=x0x4−x1x3`, `q2=x2x4−x3²`

とする。rank 四 Poisson 構造の実退化因子 `D` は既約被約 quintic とし、
`S⊂Sing D`、`pi|S` は一般 rank 二とする。
既済の nonplane secant argument から `D` の一般 normal multiplicity は二であり、
normal quadratic `G` は零でない。intrinsic conormal は可換で、`S\R` の conormal connection は flat。

**命題。** `G` が saturated normal root `ell0:N→O(C0+4f)` を持つ全 source は空である。
ここで `ell0` は `dq2` の generic line の saturation である。
normal quadratic が rank 一でも二でもよく、`R=div(pi|S)` の可約・非被約、root collision、
codimension 二の root basepoints、vertical 成分を全て許す。

凍結済 `CUBIC_SCROLL_NORMAL_ROOT_DIRECTRIX_REDUCTION.md` の全候補表と結合すると、
smooth scroll の未処理 source は次の **二型だけ**へ還元される。

* rank 一 normal quadratic：`A=2C0+5f`, `V=C0`。
* rank 二 normal quadratic：`A1=A2=2C0+5f`, `V=C0`, `Delta=C0+3f`。

本命題はこの二型、twisted cubic cone、全 G の非優越又は全 E 帰属を証明しない。
新しい成分数は生じない。

## 2. canonical lift は strict Poisson で、実 D は Casimir

`Pi` を全150 homogeneous coefficients における divergence-zero canonical lift とする。
`Pi` の係数は degree 二、`E=Σxi∂i` とする。
projective Jacobi は `E∧[Pi,Pi]=0`。
BV divergence identity と `div Pi=0` から `div[Pi,Pi]=0`。
`[Pi,Pi]` は degree 三 coefficient の三-vector なので Euler weight は零であり、
`div(E∧[Pi,Pi])=5[Pi,Pi]`（全体の符号 convention は零判定に影響しない）。
従って strict `[Pi,Pi]=0`。

標準五次元 volume と `Pi∧Pi` を contraction した polynomial one-form を `omega` とする。
`div(Pi∧Pi)=0` より `d omega=0`。
係数 degree 四なので `d(i_E omega)=5 omega`。
実 anticanonical quintic `D` は、非零定数を除いて `i_E omega` そのものである。
rank 四 bivector の linear algebra から `Pi#omega=0`。従って

`Pi#dD=0`。                                               (1)

また全 polynomial `f` に対する Hamiltonian vector field `X_f={f,−}_Pi` は
`div X_f=0`。これは `div Pi=0` と Hessian の対称性から直接も分かる。
以下では `X=X_q`, `q=q2` を使い、

`X(D)=0`, `X(q)=0`, `div X=0`                              (2)

の三恒等式を同時に保持する。勝手な Euler gauge への変更はしない。

## 3. rank 三 quadric の有限 cover と D の全許容形

`Q={q=0}` の affine cone に有限 cover を取る：

`x2=z²`, `x3=zw`, `x4=w²`。

その coordinate ring は
`A=C[x0,x1,z,w]`、`Q` の coordinate ring は `A^(±1)`。
involution は `(z,w)→(−z,−w)` である。
projective homogeneous degree を二倍し、weights を

`wt(x0,x1,z,w)=(2,2,1,1)`

とする。`Q` 上の `S` の inverse image は irreducible hypersurface

`T=x0*w−x1*z=0`, `wt(T)=3`。

cover の branch locus は `z=w=0`（vertex line の inverse image）であり、
generic `S` 上では étale である。

`dq` が normal quadratic の一根であるため、generic smooth `S⊂Q` で
`D|Q` の multiplicity は三以上である。
実際、ambient normal coordinates `(q,u)` で `D∈(q,u)²`、normal quadratic の `u²` 項が零だから
`D∈(q)+(q,u)³`。
従って cover へ pullback した `Dbar` は `T³` の倍である。
`D` が既約 quintic なので `q∤D`、従って `Dbar≠0`。
`wt(Dbar)=10`、`wt(T³)=9` から全許容形は正確に

`Dbar=T³ ell(z,w)`, `ell=a*z+b*w`, `(a,b)≠(0,0)`。         (3)

これは全係数での形であり、generic section や一つの projection の検査ではない。
`D|Q` の `S` に沿う order が四以上なら weight 上不可能となり、既に `q|D` に戻る。

## 4. parallel root と Hartogs：X の唯一の全形

flat conormal connection は parallel な quadratic `G/sigma` を保つので、
その generic root line `C·dq` を保つ。rank 一の場合も唯一の root line が保たれる。
chart `x0≠0` の intrinsic bracket numerator は

`H_i=x0{q,xi}−xi{q,x0}−2q{x0,xi}`。

projective tangency と root-line preserving は、generic `S` で

`H_i∈(q)+I_S²`。                                        (4)

従って `Q` の cover 上で `x0 X(xi)−xi X(x0)` は `T²` の倍となる。
`T` と `x0` は互いに素なので、全二つの indices に対して

`xi X(xj)−xj X(xi)∈(T²)`                                (5)

が polynomial identity として成立する。
vertex や chart の補集合を除いたままにしてはいない。

ここで `B=A/(T²)` を取る。これは dimension 三の hypersurface ring なので Cohen–Macaulay、
origin における depth は三である。
五つの `xi=(x0,x1,z²,zw,w²)` の同時零は origin だけ。
(5) により、punctured `Spec B` では `X(xi)/xi` が重なり上で一致し、一つの regular function を定める。
local cohomology の `H^0_origin(B)=H^1_origin(B)=0`、すなわち CM Hartogs によりこれは `B` 全体の元である。
nonreduced ring `B` に対してもこの depth argument は有効である。

この元の weight は四。`wt(T²)=6>4` なので `A` の weight 四 polynomial `Lambda` へ一意に lift する。
`X(xi)` は weight 六、`Lambda xi` も六だから、(5) の全形は

`X(xi)=Lambda xi+T² c_i` in `A`, `c_i∈C`。                (6)

`X(q)=0` を (6) に代入すると
`c2 w²−2c3 zw+c4 z²=0`、従って `c2=c3=c4=0`。
よって cover 上の derivation の全形は

`X=Lambda Etilde+T² V`,

`Etilde=x0∂x0+x1∂x1+(z/2)∂z+(w/2)∂w`,

`V=c0∂x0+c1∂x1`。                                      (7)

cover への derivation の lift は generic étale locus で一意であり、(6) の `z²,zw,w²` に対する
式は `X(z)=Lambda z/2`, `X(w)=Lambda w/2` を与えるので、(7) は branch locus を越えて regular。
`Lambda` の lift、全 constant vector、vertex 二方向以外の零を全て説明した。

## 5. Casimir と divergence の係数 3 対 2

`Etilde(Dbar)=5Dbar`、`V(ell)=0` である。
(2),(3),(7) の `X(Dbar)=0` は integral ring `A` で

`5 Lambda+3T V(T)=0`。                                  (8)

一方 `Q` の canonical residue volume は cover 上の
`dx0∧dx1∧dz∧dw` の非零定数倍である。
例えば `x2≠0` で

`dx0∧dx1∧dx2∧dx3 / (∂q/∂x4) = 2 dx0∧dx1∧dz∧dw`。

ambient `X(q)=0` かつ `div X=0` なので、この volume に対する induced divergence も零。
`div Etilde=3`、`Etilde(Lambda)=2Lambda`、`div V=0` より

`0=div X=5 Lambda+2T V(T)`。                             (9)

(8)−(9) から `T V(T)=0`。
`V(T)=c0 w−c1 z` なので `c0=c1=0`、更に `Lambda=0`。
従って (7) は `X|Q=0` を与える。

すなわち `Q` の ideal は projective Poisson ideal であり、`Q` は Poisson hypersurface。
rank 四 symplectic point に Poisson hypersurface は存在しないため、`Q⊂D` でなければならない。
これは `D` の既約性と次数五に矛盾する。命題が従う。

## 6. 保存する範囲と残余

この証明は frozen directrix note の canonical-root pair 四組と rank 一一組を全て空にする。
全 linear source の16次元という未閉上界を改善したという計算上の主張ではなく、
実 `D` の全形と canonical Hamiltonian の二つの恒等式による反証である。
旧二つの coarse dimension failures、三枝 rank2/2/3、後続 rootline matrix kernel16 は成功に読み替えない。

残りの `A=2C0+5f` の二型、及び cubic cone は別義務である。
全 source の非優越・分類上界・新成分・全 E 帰属を本命題から結論しない。
