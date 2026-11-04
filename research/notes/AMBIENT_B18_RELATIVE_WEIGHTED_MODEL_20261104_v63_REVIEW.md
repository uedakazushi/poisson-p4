# B18 相対 weight 模型の独立レビュー

2026-11-04。原著 `AMBIENT_B18_RELATIVE_WEIGHTED_MODEL_20261104_v63` を検証する。
実際の四項 `P0=e_y∧(x²∂h+hx∂z+(h²+z²)∂w)` を固定する torus
`(diag(1,1,λ⁻¹,1,1), tensor λ)` だけを使う。

独立した AST/Fraction sparse 演算で、以前に全物理係数が独立受理された
CA/CB の core10 源多項式を**全150成分・全次数**読み直した。
`(p,k0,…,k8)` のうち `k2` だけが character1、他は0で、全出力が対応する物理
character に homogeneous である。CB の ratio と binary f の積を含む全多項式を保持する。
定数、全一次係数、pure p²係数を以前の源証明書と照合し、全20源係数と源速度も照合した。
kernel Lie matrix は torus と可換、源速度の positive fiber 成分は零である。

物理 target と flag 微分を別に再構成した augmented source block の階数は双方29。
normalizer29 の characters は `0^24,1^5`、共通 physical19 群切断の単元 minor は `−32`、
切断 characters は `0^15,1^4`。従って残りの regular normalizer10 と共通 even base10 は
ともに `0^9,1^1`。actual odd 最小10生成元は `0^8,1^2` である。

## 全次数の群積と切断の正当化

選んだ GL 方向は torus 共役の eigenvector であり、その順序付き指数積と scalar を
formal 群 section に取れる。19個の physical 切断自体が torus-stable であり、ambient
作用写像の微分は単元なので formal inverse は存在し、一意性から equivariant である。
同じ積分解を actual CA/CB、proper source の normalizer、単射と cokernel に同時に適用する。
従って群19変数は cokernel 上にも正則であり、切断後にも target→normalizer の単射が残る。
任意の finite map の injectivity を無条件に base change した議論ではない。

切断後の実 even base を `R=R0[[w]]`、`R0=C[[x9]]` と書く。
tracefree free rank-one summand には character0 の generator τ を選べる。これは
linearly reductive torus の equivariant Nakayama によるもので、源 kernel p の character0
と整合する。従って `τ²=ν∈R0`、actual image は `A=R⊕Iτ`。
8本の weight0 generator と2本の weight1 generator を取れば
`I=I0 R+w J0 R` となる。これは actual common base 上の主張であり、二つの abstract99
normal form の相対座標を同一視していない。

次数d≥1の homogeneous first syzygy は最初の8成分が `w^d c_j`、最後の2成分が
`w^(d−1)d_j` の形となる。regular domain で `w^(d−1)` を除くと degree1 syzygy。
degree0 syzygy と併せ、first syzygies は degree≤1 で生成される。
全代数 kernel はこれらと `Z_i Z_j−ν b_i b_j` の55積関係で生成される。
積関係で任意の式を `R+ΣR Z_i` に還元し、basis `(1,τ)` を使えば残余関係が
ちょうど first syzygy となるため、有限 jet に限った主張ではない。
したがって各成分の最小20変数 presentation は positive degree≤2 でよい。
受理済みの最小99関係数上界は formal smooth 群因子の除去と equivariant Nakayama で残る。

## 全 D では87 graph差を保持する

canonical126 の character 数を独立に数え直すと `0^36,1^60,2^24,3^6`。
共通 actual tangent39 は `0^32,1^7`。従って最小39座標外の87個の物理 ambient 座標は
weight0/1/2/3 にそれぞれ4/53/24/6個である。
最小座標を含む equivariant Cohen presentation を取り、余剰座標の homogeneous lift を
選べる。零 weight 変数の形式的函数を係数にすると、weight d の lift は positive 変数に
ついて次数dの多項式である。87 unit linear relations の formal IFT による消去とも一致する。
lift の kernel modulo の非一意性は、この存在と次数上界に影響しない。

CA の graph で余剰87変数を消した全 D には、両 component kernel と**全87 graph差**を残す。
次数≤3、生成元数≤`99+99+87=285` が安全な全次数上界である。
正重み最小変数は3個、零重みは17個なので、係数環 `C[[X17]]` 上で
`1,3,6,10` 種類の monomial に付く formal 函数へ問題が有限化される。
可能な三次 graph差6本を落として degree≤2 とすることは未証明である。

## 固定部分の範囲

非負 grading により、positive3変数を零にした fixed quotient は degree0 subring。
その cotangent dimension は17。既受理 nil-cotangent の2 class
`[u46+u95]`、`[u63]` はそれぞれ weight1、0である。
完備な torus-stable nilradical の weight0 射影は同じ非零 linear class を保つ。
`(m_D²)_0=m_fixed²` と subring の nilradical の一致から、fixed 部分の nil-cotangent 像は
正確に1次元で `[u63]` が基底となる。生の `u63` 自身の nilpotence は主張しない。

零17変数の実係数函数、相対 graph差、全 D の nil ideal・socle・depth19/20 は未計算。
positive3変数の正則列性・flatness を仮定せず、fixed quotient の情報から全 D の depth を
移すこともしない。新点、源 jet の追加計算、GB、消去、原著 producer 再実行は用いない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_relative_weighted_model_review_20261104_v63.py --output work/v62-review/B18-relative-weighted-independent
```
