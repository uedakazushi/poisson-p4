# Spin(4) 二形式切片と混合判別式の独立監査

2026-11-03。対象は `Q_FIRSTJET_SPIN4_MIXED_DISCRIMINANT_20261103_v63`。
先行する Spin(5) 構造とその独立監査を受理済み入力とする。
今回受理するのは、実12次元切片の同定、混合判別式 D の既約性・二重次数(8,8)、
line normalizer に対する偶性、p0!=0 上での regular 不変式への降下である。
極の位数 k、24次 polynomial lift、NEW25 との対応は受理しない。

## 1. 物理的切片を別実装で再構成した

author module は import していない。元の24係数の順序を2形式の組と線形変数から作り直し、
対称な A-spinor / B-spinor 係数を1/2で戻して全24×12埋込みを Fraction 演算で再構成した。
閉条件、primitive 条件、全 gamma trace を直接確認した。
cross 部分16係数への gamma map の階数は4、構成した埋込みは階数12である。
従ってこれは gamma-traceless な e0 行零切片の全体である。
受理済み24×16基底への埋込みも全係数で照合した。

元の六つの block diagonal sp2⊕sp2 生成元に対し、物理的な2形式の Lie 微分を作り、
U(2,1)、V(1,2) の各単項式を変数微分する別の係数計算と比較した。
全6×24×12係数が一致する。元の pullback は反作用なので、Lie 作用を通常の左作用の
符号へ黙って変更しない。不変式への結論は連結な群の作用について述べる。
有限 normalizer g も元24係数へ実際に作用させ、

`(U(a,b),V(a,b)) -> (−V(b,a),U(b,a))`

の全係数を確認した。cross4×4 determinant delta の全係数を独立な subset determinant で
作り、g*delta=−delta を検査した。

## 2. 制限行列と判別式回路

V=a0*v0(b)+a1*v1(b) に対して、v0,v1 が共通零点を持たない開 R=Res2(v0,v1)!=0 では
V=0 は b を底とする graph a=(-v1,v0) である。
U の制限 F5=u0*v1²−u1*v0*v1+u2*v0² を convolution から作り、六つの全係数を比較した。
6×6制限行列の determinant は、resultant API を使わず subset determinant によって
全係数で R³ と確認した。従ってこの開で U の六係数と全 binary quintic の六係数は同型である。
幾何的にも restriction sequence の核・余核は O(1,−1) の H0,H1 であり、いずれも0である。

一般 quintic の判別式59項も author の discriminant 計算を実行していない。
F と F' の9×9 Sylvester determinant を整数係数の subset DP で計算し、leading coefficient
で割って全59係数を新たに得た。保存公式と全項が一致する。
170-node 回路は六つの F5 係数と R の節点で切り分け、下側は12係数の次数4以下の辞書、
上側は F0,...,F5,R の七記号の Laurent 辞書で全恒等式を検証した。上側の出力は正確に
上記59項とその R^-2 倍であり、元入力から cut を迂回する経路も無い。
従って二つの fixture だけから回路式を外挿していない。
全節点、二重次数、保存済み二つの normalization fixture も Fraction で検査した。
12変数の合成判別式は展開せず、新しい first-jet 点も評価していない。

## 3. R² の整除と正確な付値

quadratic resultant R は既約である。共通零点を P1 に記した incidence は既約で、その一般点は
唯一の共通単根を持つ。したがって像は既約 hypersurface である。一般の横断変形で R は一次の
uniformizer になるため、resultant がその hypersurface の高い冪であることもない。

R の一般点の DVR を取り、共通根を t=0 とする。二つの quadratic は
v0=t*a(t)+epsilon*p(t)、v1=t*b(t)+epsilon*q(t) と書けるため、制限 quintic F5 は
局所 ideal (t,epsilon)² に属する。一般点では epsilon=0 の F5 は t=0 にちょうど二重根を持ち、
残り三根は単純である。二次 Weierstrass 因子は

`t² + epsilon*B(epsilon)*t + epsilon²*C(epsilon)`

となり、その判別式は epsilon² で割れる。残りの因子・cross resultant は unit である。
よって Disc(F5) の R-adic valuation は少なくとも2である。
全係数環は UFD なので R² は多項式として整除する。

独立 code は保存された小 fixture
v0=t、v1=t²+epsilon、(u0,u1,u2)=(1+t,3+2t,7+5t) を、上で独立に得た59項へ代入した。
homogeneous quadratic resultant は正確に epsilon で、quintic 判別式の epsilon 次数最小値は2、
epsilon²で割った定数は1466325!=0である。従って R³ の普遍的整除は不可能であり、付値は正確に2。
この fixture は一般 DVR の下界証明に代わる外挿ではなく、付値上界を決定する検査である。

従って D=Disc(F5)/R² は全12係数上の多項式で、二重次数は(8,8)、合計16となる。
Rを可逆にした環では U->F5 が可逆な線形変換で、普遍 quintic 判別式は既約である。
D は R 因子を持たないため、全係数環でも既約である。

## 4. normalizer の偶性に必要な幾何を確認した

P1×P1 の点と、U=V=0および二勾配が従属する条件を記した incidence を考える。
O(2,1)、O(1,2) の各1-jet評価は全射である。点を固定した fiber は二つの独立な値の条件と、
自由な2×2 gradient determinant の既約 hypersurface である。したがって incidence は既約である。
滑らかな V 上で F5 が一般の一つの二重根を持つ点を含むので、その像の閉包は上の D divisor である。
この記述は二つの P1 と U,V の交換に対称であり、normalizer g はこの divisor を保つ。
既約性から g*D=cD である。

g の二つの保存済み fixture では、元と変換後の R は50、37で、D はどちらも53268630696である。
したがって c=1。この一点は既に分かっている比例定数を決めるためにだけ用いる。
連結成分 Spin(4)=SL2×SL2 に対しても divisor が不変であり、半単純連結群には非自明な
character がないため D 自体が不変である。
Spin(5) の nonisotropic line stabilizer N は二成分を持ち、恒等成分が Spin(4)、他方が g の coset
である。従って D は実際の N 全体で不変となる。

## 5. p0-open の同型と regular 降下

先行 Spin5 の Q は det Q=1/16、det J^-1=1である。
Q の全25個の4角 cofactor 恒等式を Fraction 演算で独立に検査した。
Cauchy–Binet とこれらの compound identity を用いれば、Psi の signed maximal-minor 列 chi に対して

`chi^T*Q^-1*chi = det(Psi^T*Q*Psi)/det Q =16*p0`

が任意の Psi で成り立つ。20変数の8次式を展開する必要はない。
n=Q^-1*chi とすると、n は image(Psi) の直交補に入る。
slice では chi=delta*e0、ゆえに p0=delta²/16 となる。

p0が単元なら Gram matrix が可逆なので、Psi は split rank4、signed maximal minors は line bundle
を生成し、normal line は nonisotropic である。したがって normal line は Wのこの開上で
**regular かつ一意に** 復元できる。nonisotropic line と、その直交補へ像を持つ Psi を記した
incidence からの忘却写像は W[p0^-1] への同型である。

G=Spin5=Sp4 は nonisotropic line の空間 G/N に推移的で、e0-line 上の fiber はちょうど
L[delta^-1] である。G->G/N の faithfully flat な principal N-bundle による descent から

`G x^N L[delta^-1] ~= W[p0^-1]`

を得る。これは generic orbit の同定に留まらず、この open 全体の scheme 同型である。
Zariski 局所切断や normal vector の square root の大域選択は仮定しない。

D は先に全 L 上の多項式へ延長済みなので、R=0 の slice 点を除かずに N-invariant regular
function として associated bundle へ降下する。従って W[p0^-1] 上で regular な Sp4 不変式となる。
scalar倍は normal line を変えず、D の次数16と可換であるため、降下した有理不変式の斉次次数は16。

W は affine space、p0は8次 polynomial invariantなので、その関数は P/p0^k と書け、
Pは不変かつ斉次次数16+8kに取れる。これは principal-open coordinate ring の性質である。
**kの上界は証明されていない。** k=1の24次 lift、k=0の可能性の排除、NEW25上の消滅、
既存の24次曲線判別式との一致は全て独立の未解決課題である。

## 再現

`code/q_firstjet_spin4_mixed_discriminant_review_20261103_v63.py`。
60秒・2GB上限で、author moduleを呼ばずに元24係数、small polynomial dictionaries、
Fraction、Sylvester subset DPだけを使う。受理範囲と入力 hash は同名 report に固定する。
