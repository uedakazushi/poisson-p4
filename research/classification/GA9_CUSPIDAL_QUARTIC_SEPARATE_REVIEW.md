# Ga9 power1/power3 を含む新24成分：独立監査

2026-11-07、PRIMARY E。**全項目 PASS**。
全七 parameter の実族閉包は24次元の全 Poisson 成分を与え、凍結した prior102 と相異なる。
Ga9 c=0 の power1/power3 の全源閉包はこの成分に含まれる。
local uniqueness、generic smoothness、Ga9 linear-w又は型6の全支持被覆は主張しない。
最終 claim ID と全体の計数は root が行う。

## 1. 全七源と Saito frame の独立再構成

checker は `code/classification_ga9_cuspidal_quartic_review.py`。
producer、SymPy、CAS を import せず、別作者の univariate Fraction/gcd による QQ(a) と
疎多項式で全係数を扱った。base=(h,y,w)、fiber=(x,z)。

\[
 F_1=(3hw+2y^2)^2-16y(h+ay)^3,\quad
 F_3=(hw+6y^2)^2-144y^3(h+ay),\qquad D_i=hF_i.
\]

全18 quadratic gradient-syzygy 係数の行列を gradient(Di) から別に構成した。
行数は22と20、a=0で rank16、保存16小行列の a=0 値とも一致した。
producer の大きな rational K を再利用せず、matched U,V と commutator から

`λU=divU+β`, `λV=divV−α`,
`SU=U−λU E/5`, `SV=V−λV E/5`

という別の二本の gradient syzygy を構成した。全 gradient との積が零、二本が独立であり、
generic に全 syzygy kernel を尽くす。
元の全18係数を保持した rank16＋kernel2 であり、選んだ二導分だけを ansatz としていない。

全 parameter のまま

`Uh=Vh=0`, `[U,V]=αU+βV`,
`det(U,V)/F1=12/(2a³−3)`, `det(U,V)/F3=252/(18a−7)`

を直接検査した。従って Euler を加えた determinant は unit*hF。
reduced D の chart では Saito criterion が全 logarithmic basis を与える。
探索時の他の16minor因子を反転する必要はなく、最終 source の分母が実際に
`2a³−3` 又は `18a−7` の因子だけであることを QQ[a] gcd で確認した。

任意の六係数を残した
`R=r0h²+r1hy+r2hw+r3y²+r4yw+r5w²` に対し

\[
 P=\partial_x\wedge U+\partial_z\wedge V
       +(\alpha x+\beta z+R)\partial_x\wedge\partial_z
\]

の全150係数が producer と一致した。
全10 strict Jacobi が QQ(a)[r0,…,r5,h,x,y,z,w] 上で零なので、全350 projective Jacobi も零。
六磁気を特殊化した後の確認ではない。

## 2. 全源の cofactor 構造と、openness 文言の訂正

**N²=0 自体は開条件ではない。** producer の「冪階数2,0とh条件が一般parameterで開」
という一文は、次の構造的恒等式と真の開条件に分けて用いる。

上の全 P で発散 Z は fiber ∂x,∂z 成分だけを持つ。
各成分は base h,y,w の linear form であり、Z は base を殺す。
従って全7parameter の5×5 cofactor行列は base→fiber block のみ、**恒等的に N²=0**。
これは独立実装でも全有理係数のまま検証した。
a=1 の点で rankN=2、h∉im(Nᵀ) を再生した。この二条件だけが真の開条件である。
base bracket {y,w}=0、h-Casimir も全源の構造に含まれる。
従って一般源は既済の型6 nil22 image・plane-quartic G 上界≤24の仮定を満たす。
端点 a=0 で distinguished-H 条件が変わる可能性を捨てていない。

以下の真切片の局所上界だけでも成分性は証明出来るので、この G 上界は追加の整合確認である。

## 3. 全磁気端点と二実族の full150 transfer

両実族の分母は a=0 で非零。r5=8又は72、他の r=0 とした a=0 の全150係数は、
公開 Ga9 normalform の power1/power3 canonical tensor と厳密に一致した。
a=0 で任意の r を戻すと、差はそのまま `(R−Rspecial)∂x∧∂z`。
従って端点でも全六磁気 affine fiber を保持しており、磁気零も同じ実族に入る。

a≠0 上の base 変更

`H=h`, `Y=−(h+ay)/a³`, `W=9w/a⁴−12y/a⁵−6h/a⁶`

の determinant は −9/a6。
`F3(a³)(H,Y,W)=9F1(a)(h,y,w)/a8` を QQ(a) 上で検査した。
さらに報告された2×2 anchor matrix T と fiber matrix が厳密に T⁻ᵀ であることを確認し、
完成した5×5 PGL matrixによる **全150括弧係数**の pushforward を別に展開した。
結果は F3(a³) の source に保存された六磁気変換を代入したものと完全一致する。

六磁気の線形 map は加法定数項を持たず、全六係数を保ち、その determinant は正確に

\[
 \frac{117649a^{12}(2a^3-3)^6}{9(18a^3-7)^6}\ne0.
\]

a=0 を transfer の式で割ってはいない。
generic transfer により F1 の稠密 parameter 開部分が F3 source閉包に入り、
その閉性から F1 の全 regular parameter と a=0 にも延びる。
F3 は自身の実族で a=0 まで正則である。
Ga9の全十一源の飽和閉包が canonical23軌道閉包そのものであることは、
受理済み `NIL211_SHORT_TRANSLATION_REVIEW.md` の全源上界と actual orbit23 を依存とする。
そのため二 canonical点の帰属から、両 power の全 parameter 源と磁気零境界まで被覆出来る。
零磁気が非零磁気と同じ orbit だとはしていない。

## 4. 複素絶対既約性と一般 member

w の二次 discriminant は独立に

`576h²y(h+ay)³`, `576h²y³(h+ay)`

と一致した。先頭係数は h²の定数倍で、h=0 の定数項は
`4(1−4a³)y4`、`36(1−4a)y4`。
各 primitive unit が非零なら係数 gcd=1。
discriminant は C(h,y) 上で平方でなく、a=0 でも h5y 又は h3y3 の奇数次数因子を持つ。
Gauss により F は複素数体上で既約。h は F の因子でないので一般 hF は被約。

F3 の normalization は
`[h:y:w]=[(s0²−as1²)²:(s0²−as1²)s1²:12s0s1³−6s1⁴]`。
代入で F3=0、a≠0,1/4 で basepoint がなく、
`t=s0/s1=((w/h)/(y/h)²+6)/12` から generic birational である。
t=±√a の二点は `(0:0:1)` に写る。両枝は smooth で共通接線 h=0、
h/y² の値は ±12√a−6 と異なり A3 tacnode。
s1=0 の affine orders は2,3で A2 cusp。
delta2+1が quartic の算術種数3を使い切るため他の特異点はない。
Hは二 tacnode枝に各二重に接する。
これは一般 member の記述であり、除外された parameter の配置を同一視しない。

## 5. 真の49切片と非零二次余核による成分性

a=1、六 r=(1,2,3,5,7,11) の点とその integer scale11を全係数で再生した。
全150変数/350行 D、全 GL5/vertical/scalar M と全7 source微分を加えた A は

`rankD=100`, `rankM=48`, `rankA=49`, `DA=0`。

a 微分は QQ(a) の numerator/denominator を正確に微分したもので、数値差分ではない。
全350行の積零、保存100/48/49小行列を確認した。
実源微分49方向から source/PGL 飽和閉包の射影次元は24以上。

独立49 source列を F とし、その49 coordinate行を cut とする。
保存 R は全49cutで零、DR=0。
F⊕R は50次元で、Dの100可逆 coordinate方向 U100 と合わせて全150空間になる。
独立 checker は U100 の補 coordinate50行上の det(F⊕R)≠0 を直接確認し、
保存 D100 minor と合わせて真正150直和を証明した。
U100 が最初からcutに入ると仮定せず、必要なら
`U'=U−F*(F_cut)⁻¹U_cut` と project する。DF=0 なのでDU'=DU、可逆100minorは変わらない。
従って101次元の真切片は U'100 と一つの R で張られる。

R の全二次 projective Jacobi Q を直接構成し、保存350係数と一致させた。
`[D|Q]` の非零101小行列を再生したので、Qは全 imD の外。
selected100式の形式陰関数で U' を消すと、一変数 t の残式の二次係数はこの非零余核 class。
従って残式の少なくとも一つは `t²*unit`、切片の被約支持は原点である。
高次を求める必要はなく、正重みの完全多項式 IFT 又は nil指数計算は使わない。

49本のlinear cutに対するKrull不等式から全raw局所次元≤49、vertical24とscalar1を除いて≤24。
実 source は irreducible な七parameter域と connected PGL の像であり、その閉包 Y は既約。
rank49の下界と合わせdimY=24。
Yを含む全成分は上の点を通りdim≤24なので、その成分はY自身である。
source49には非群作用方向もあるため、局所成分の一意性やsmoothnessまでは主張しない。

一般 actual PGL orbit の次元は正確に23である。
全 strict Poisson source に対し BV divergence identity は
`[divΠ,Π]=0` を与える。上で確認した一般の divΠ は非零 nilpotent rank2 の linear field で、
projective scalar fieldでない。従って一般 stabilizer は少なくとも一次元、orbit≤23。
anchor の実行列 M48 は orbit≥23 の開条件を与えるので一般 orbit=23。
新24成分の一般 orbit codimension は一、Rosenlicht quotient の次元も一である。
parameter a がその quotient のbirational座標であることや族写像の次数までは主張しない。

## 6. 全 prior102 と相異なること

commit `3fdd4295f818338f23079fa6201bfae6f291e390` の COMPONENTS blob
`f588ae1d60364ca771f746a46cc23fa5e11e6cf9` を read-only で取得し、SHAと全102行が
producer の immutable snapshot に一致することを確認した。

新 source の D は全 parameter で三 base変数にしか依存しない。
actual D の方向微分 map の rank≤3は閉条件であり、PGLと閉包で保たれる。
同rank・同次元の MPS八件 T02,T03,T16,T17,T21,T36,T37,T38について、
公開全sourceの edge を零にし、全 period を残したactual tensorを別に構成した。
各 source の全 logcanonical係数と signed Pfaffian を一致させ、
`D=κ(p)x0x1x2x3x4`, κ≠0 を確認。
五つの missing-coordinate monomialを行に取る方向微分 minor は正確に κ5≠0。
各named sourceがその全MPS成分に入るという既済識別を用いると、
全成分が新しい閉rank≤3 locusと一致することは不可能である。
一般因子の未証明推測や、同次元でない MPS23監査の流用ではない。

残る三同次元成分は全体が閉 `1+2+2`、非被約、`1+1+1+2` locusに入るという既済事実で区別する。
新 source の一般 D は単純 H と複素既約 quarticなので、そのどの閉locusにも入らない。
他の全項目は全成分次元又は閉rank2条件で区別される。
これで全102件に理由を割り当てた。

## 7. 証拠と未処理の範囲

- 有限契約：`reports/CLASSIFICATION_GA9_CUSPIDAL_QUARTIC_REVIEW_JOB_20261107.json`
- 独立結果：`reports/CLASSIFICATION_GA9_CUSPIDAL_QUARTIC_REVIEW_20261107.json`
- 全入力hashと独立実装全文を結果に保持。

gate PASS後1.6291秒、最大RSS22,456KiB、全再構成PASS。
sourceのSaito chart分母、primitive unit、generic transferのaと18a³−7を全て記録した。
excluded chartを型6全支持から削った証明ではない。
Ga9 linear-w と他の型6 plane-quartic源は未処理。
producer・既存証明書・共有台帳・Git・原稿は編集していない。
