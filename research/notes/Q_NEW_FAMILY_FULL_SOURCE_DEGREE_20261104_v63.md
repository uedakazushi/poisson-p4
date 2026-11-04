# Full normalized sourceの正確な次数：8、tau下降後4

2026-11-04。対象はactual reduced NEW S61像Y16と、そのprimitive β像B15である。
既受理の有理性の議論で未決だった、sourceの一変数有理被覆の**正確な次数**を決める。
新しいMöbius展開・点・CAS・不変式jobは使わず、実際のnormalization dataの一意性と
全normal form同定、generic full stabilizerを数える。

結論はk=QQ,Cで

\[
[k(Sp4)(u,v,c,X,W,B):k(Y16)]=8,
\]
\[
\boxed{[k(B15)(z):k(Y16)]=4},\qquad z=u(t-2)/t.
\]

これは第一jet invariant-fieldの既知指数1と異なるmapの次数である。
またNEW方程式の次数、spectral次数10/20/40、universal曲線の次数8とも区別する。

## 1. 同じ実際のfull source map

元の三parameter族をn1=1とし、modular Gaで評価点y=0へ移した六sourceを
F6=k(u,v,c,X,W,B)とする。点は

\[
p=(X,0,Z,W),\quad
Z=-\frac3{4u}(B+W/1536-W^2+X^2/64+2/9).
\]

同じsourceに全Sp4を掛け、既存のDarboux frameとS61 perspective normalizerを適用する
写像Xold⇢Y16を使う。Xold=Sp4×Spec F6は16次元、既存のfull61 mapと有限gaugeにより
Y16へdenseである。よってgenericにfiniteで、標数0なので分離的である。
Sp4のbig cellを使う実装と全groupを使う函数体は同じ。以下ではgeometric generic
fiberの点を数え、次数を得る。

## 2. sourceは実際のnormalized equivalence dataと一対一

固定された一般Π∈Y16に対し、projective transformation Tとscalar cが

\[
\Pi=c\,T_*P(u,1,v)
\]

を満たすとする。Πはe0でrank4、value Jであるから、p=T⁻¹e0はsourceのdegeneracy divisor
3H+Qの外にある。特にpはHのaffine complement内にあり、finite pointである。

T(p)=e0を満たすprojective mapは、分母のpでの値を1にして、一意に

\[
T=\operatorname{Persp}(\alpha)\,
  \operatorname{block}(A)\,\operatorname{Translate}(-p)
\]

と分解できる。Aはpでのaffine derivativeである。もしq=x−pについて
T(q+p)=Aq/(1+bᵀq)と書けばalpha=A⁻ᵀbであり、隠れたscalar選択は無い。
PGLの作用はprojective quadratic bivector上で線形で、GL5のcentral scalarは
degree2係数と二つのvector indexのscaleを相殺してtrivialに作用する。
Π≠0なので、このTに対応するPoisson scalar cも一意である。

C(p)を元Pの値とする。target value Jの条件は

\[
A(cC(p))A^{\mathsf T}=J.
\]

既存の**有理**Darboux frame RもR(cC(p))Rᵀ=Jを満たすから

\[
g=A R^{-1}\in Sp4
\]

は一意に決まる。cの平方根は使わない。perspective alphaも、S61の四つのgauge条件と
invertibleな4×4 response matrixにより一意である。これはprimitive βだけでなくfull tensorの
normalization theoremである。

逆に(u,v,c,p_y=0,g)から、このR、A=gR、唯一のperspectiveを使いTを復元できる。
z点座標とBはu≠0でbirationalに対応する。valueのPfaffianは8B/9なので、rank4 pointでは
B≠0でありDarboux chartを失わない。従ってsource pointと、y=0を課した上記equivalence dataは
genericに一対一である。

## 3. 正確に四parameter代表×二Aut component

whole normal form orbit theoremから、一般のΠに対応するparameterは正確に

\[
(u,v),\quad(-u,-v),\quad(u,-32u-v),\quad(-u,32u+v)
\]

の四つ。一般J≠0では相異なる。各代表sについて、Tの全選択は
AutPGL([P_s])のtorsorであり、genericな群はGa×mu2である。

modular Gaのsource座標作用はy→y+4u q/3。Tを右からexp(qN)と合成すると、
p=T⁻¹e0のy座標はp_y−4u q/3へ移る。u≠0なので、各Ga component内でp_y=0へ移すqは
正確に一つで、その交点はreducedである。finite Eのcomponentでも、固定Eを先に適用した後の
同じGa平行移動となるため、一度だけ交わる。

従って各parameter代表から二つ、合計**八つのgeometric generic source point**がある。
dominant同次元mapのchar0 separabilityにより

\[
\boxed{[k(Xold):k(Y16)]=4\cdot2=8}.
\]

次の二重計数をしない。

* sigmaの二つの表示M±は、互いの差がGaEに属する。これは上の二Aut componentそのもので、
  さらにfactor2を掛けるものではない。
* GはSp4でありPSp4ではない。Aを固定すればg=AR⁻¹は一意。
  −gは一般に違うprojective derivativeで、primitive βを反転する。
  generic β stabilizerはidentityなので、中心から新たな同じtargetへの重複は発生しない。
* scalar cはTに一意。Darboux frameはcについて有理で、平方根二択を持たない。
* n1=1のtorus gaugeは既に一度固定済み。残るfinite parameter同定は上記四つに含まれる。

点がnormalizer chartを外れる問題も、generic degreeを変えない。
rank4条件はpをaffineかつB≠0へ置き、その他のsource座標のbirational chartの補集合は
16次元sourceの真の閉集合なので、generically finite mapで一般targetを支配できない。

## 4. full tau下降は新たな二次拡大を加えない

受理済みfull-tensor fiber packetは、F6内で有理なraw-pair matching S_tを作り、
g_u=diag(u,u⁻¹,1,1)を使って全tensorをtau固定にした。ここでuは**元のF6のcoordinate**であり、
この段階で新しい平方根を付加するのではない。

perspectiveはS61 quotientで除かれ、effectiveなphysical Sp4 frameの差はg_u⁻¹S_tである。
従ってold group gとnew group hは

\[
g=h g_u^{-1}S_t
\]

というF6上の有理な可逆座標変換で結ばれる。source函数体k(Sp4)F6そのものを変えない。
physical pushforwardをβのpullback anti-actionに合わせるgroup inversionもbirationalであり、
次数を追加しない。

このnew group hを固定して、tauはu→−u、t→t/(t−1)と作用し、全full tensorがtau固定である。
F5=E5(u)、u²=1600r、E5=k(xi,r,s,q,B)、F6=F5(t)に対して

\[
z=u(t-2)/t,\quad t=2u/(u-z),\qquad F6^\tau=E5(z).
\]

tauはuを反転する非自明involutionなので、[F6:E5(z)]=2。従って

\[
S_{new}=k(Sp4)E5(z),\qquad [k(Xold):S_{new}]=2.
\]

既受理の全次数tau同値が、このSnewにk(Y16)を含める。βだけの偶性から
full tensorの下降を推論しているのではない。

## 5. 同じB15上のdegree4 Lüroth map

βbar sectionとgeneric Sp4 freenessにより、同じfield embeddingで

\[
k(B15)=k(Sp4)E5.
\]

physical group variableには必要に応じh→h⁻¹を用いて、既存のpullback規約のbirational productと
合わせる。別のbaseの同型を暗黙に同一視しない。
従ってSnew=k(B15)(z)であり、field towerは

\[
8=[k(Xold):k(Y16)]
 =2\,[k(B15)(z):k(Y16)].
\]

よって

\[
\boxed{[k(B15)(z):k(Y16)]=4}.
\]

K=k(B15)とすれば、既知Lüroth theoremはk(Y16)=K(w)を与え、このsource mapは
degree1ではなく正確に**degree4**のrational curve mapである。
具体的なw(z)の係数、Galois性、全special fiberや全域finite性は主張しない。
QQ上のmapのgeometric generic fiberをchar0で数えているので、QQとCの両方に同じ次数が成立する。

## 保存と未使用の次手

source `code/q_new_family_full_source_degree_20261104_v63.py` は独立gatesとhashの確認、
4×2=8と二次towerの算術、証明schemaの保存だけを行う。
CAS、新しい点・prime、Möbius係数展開、source26、不変式jobは無い。

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_full_source_degree_20261104_v63.py \
  --output work/full-source-degree-replay
```

v=0からの別の五parameter sectionは、δ定数mapのboundary extensionを別途厳密に確認すべき
新しい方法であり、今回のdegree8/4証明には使っていない。
