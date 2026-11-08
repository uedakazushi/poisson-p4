# Concurrent m三・nonzero binary quartic 全源の独立確認

2026-11-08、PRIMARY G、別作者 `/root/quartic_review`。
producer の strict Ore bridge、binary Q 全正常形、central 全補集合、
K211/K31/K22/K4a/K4b の全証明を元式から読解・手で再導出した。
**指定 whole Ore 源の K4 非零部分は PGL5 像次元高々24**と受理する。
加えて linear normal Casimir を持つ unimodular K0 の四指定部分も全空と独立確認した。
新 CAS、operator/module、minor、source probe、jet、GB 又は数値実行は無い。
gate は保存入力の hash を検査するだけで、数学は以下の独立証明確認による。

## 1. 対象、全 higher、strict Ore への矢印

仮定は rank 四 projective Poisson、actualD が被約既約 quintic、
rank 二 invariant plane、normal 位数三、restriction zero cubic が三被約 concurrent lines。
その全 first-cover の Ore 部分 P0=∂c∧W0 を対象にする。
全 raw36/Euler6 P1 と全九 P2 を戻し、初めから higher を c-independent としない。
別 c-dependent first profiles、K4零、他 normal 位数は本 receipt 単独の結論に含めない。

canonical divergence-zero strict lift は normalweight を保つ。
canonical W の normal 部分が零なら、raw M=(tr M/5)I2 から M=0。
この scalar conormal source の parallel G3 は dg=(2/3)g dlog f を満たすが、
各 simple line の valuation が2/3となり polynomial G3!=0 と矛盾する。
従って W の少なくとも一つの normal component は非零。

canonical P1 の全 c 依存は Q1+c∂c∧L+c R、R=ell(u,v)∂a∧∂b。
strict weight 一 Jacobi の非∂c部分 W∧R=0 と domain 性から R=0。
divergence は L=−δ_yQ1、Euler E∧δ_yQ1 を戻すと全 c 項が消える。
P2 は元から c 非依存で、この操作は全九の affine particular を捨てない。
ゆえに全 projective class は

    P=Q(y)+∂c∧W(y), y=(a,b,u,v),
    [Q,Q]=0, E_y∧L_WQ=0

へ入る。逆も c 係数を比較して正しい。W は plane ideal(u,v) を保つ全 quadratic vector field。
この ideal 条件を無視して arbitrary W と呼んだ後に weight 条件を使うことはしない。
Q の全源は Qab=aA1+bB1+C2、Qan=U2、Qbn=V2、Quv=0、
[U2,V2]=A1U2+B1V2。C2 全三係数も残す。
actualD=c Pf(Q)−[E_y∧Q∧W] であり、orientation による全体定数だけを許す。

## 2. nonzero K4 の全根分割と actual frames

K=det(U2,V2)!=0 に対し U(K),V(K) は K の倍数なので、両 field は各 prime normal line に接する。
projective component uQ2−vP2 は cubic。四根を持つと両 field が linear multiple of Euler となり det=0。
三根では全 quadratic log fields は ell E_n+t Jac(g)。その determinant は g times linear なので211。
二根は22/31、一根は4であり、これ以外を除いていない。

二根の linear pencil L=uA+vB には base GL2 の左作用だけを使う。
22 は A/B の rank 一 image lines を正常化して全 ordered P1 row pair、det!=0。
31 は左 A^-1 で uI+vC、tr C=det C=0 とし、C零と全 projective nilcone P1を残す。
根の方向を変える任意 right GL2 を仮定しない。
一根は normal u-column rank 二/一の二場合から K4a/K4b の全二parameter族になる。
全操作は constant block GL と tensor scale、逆は通常の inverse matrix。
複素非零 scale が必要な箇所はその零補集合を別根分割へ残す。

central shift の四列は −V(u)−A1u, −V(v)−A1v, U(u)−B1u, U(v)−B1v。
この式を全 K4a/b 正常形へ代入して rank 三/二/一の境界を確認した。
K4a の β!=1 だけで C2 を全て除け、β1 では元全 C2 を pole 証明へ戻す。
K4b は β!=0: C2零、β0/α!=1: C2=0又はv²、β0/α1:0,uv,v²。
special residual xuv+yv² は actual v→v+ku（必要な paired base shear 込み）で x−2ky と動く。
full W へも同じ可逆作用を施すので、magnetic normalform は fullsource の slice ではない。

## 3. universal full module と simple-root 二枝

chart v1 で piij=Qij+yi Zj−yj Zi、Zi=Qi,v と直接 dehomogenize すると

    P_chart=pi+∂c∧(Y+cZ), Yi=Wi−yi Wv,
    L_Ypi+Y∧Z=0.

最後の符号は plus。strict symplectic Q の inverse は homogeneous weight零かつ closed。
β=±i_E omega の radial quotient は closed、β(pi)=0、β(Z)=1。
local r primitive と Z-flow で pi=r-independent、Z=∂r とし、元 module の r/leaf 成分を比較すると

    Y=H Z+pi#dH, H=β(Y), D=K4(c+Hhom)

を得る。残 leaf 成分は恒等式、逆も正しい。
ここで analytic/log coordinates は恒等式の導出用であって actual PGL ではない。
global polynomial W→Y の kernel は W=ell E_y、ell linear、従って同じ projective tensor class。
これは polynomial gcd から従い、局所に選んだ Wv=0 lift が global regular だとは言わない。

K211 の全 t,C2 で simple root u=v の H=h/(u−1)+regular を使う。
pi_au=0、Za=u、pi_ab|u1=C2(1,1) より Ya の residue は h+C2(1,1)h_b。
有限 polynomials 上 1+d∂b は可逆なので h0、actualD は u−v で割れる。
K31 の全 nilcone/零で chart u1 の pi_av=−v+O(v²)、Za=1+O(v)、
pi_ab|v0=C2(1,0)。Ya の residue は2h+C2(1,0)h_b。
同じ理由で h0、D は v で割れる。
全 C2 と全 W を残した whole factor であり、multiple pole には外挿しない。

## 4. K22：全 ordered rows と全 Laurent 源

p=A−B,q=C−D、全 f=f0a³+f1a²b+f2ab²+f3b³、
全 normal-linear matrices を戻した最高 base-degree module の四式(I)–(IV)を直接微分して確認した。
full higher W/C2 はこの次数に寄与しないので W0 slice の同値主張を使わない。
p q!=0 では I constant/II u と actual swap から f0=f3=0。
IV u³/II u³ は R2=0,A S2=0、続いて A,C!=0 は f1=0、swap で B,D!=0 は f2=0。
reduced f は f1f2!=0 なので diagonal 又は anti-diagonal のみ。
p0 では S=0,f3=0、R=−3q f2/B constant。III constant から D0が必要、radial/axis 一点になる。
q0 は actual swap、p=q0 は determinant零で source外。
diagonal は W0u が u 倍、唯一の higher 非可除項 λv² は元(b,u) module の −2Dλ から消える。
radial/axis の対称項 λu² も元 module の −λ から消える。
全九 P2 と全 higher を保って actual proper factor となる。

残 anti-diagonal を actual row scale/C2 shift で Qab=−av+bu+γuv,Qav=Qbu=uv、γ0/1へ。
β=d(a/u+b)+γ du/u、全 Hhom は numerator normal≥3/degree五、denominator u²v²。
chart H=Σ(k=−2..3)h_k u^k、h_-2..h1 は base-degree≤2、h2≤1,h3 constantの全28項。
u-pole から h_-2=Aa²、h_-1=3Aab+g(a)、a g'−g=2γAa。
従って γA=0,g=Ba。actual swap と homogeneous degree 一 twist を戻すと全必要 H は

    A(a²v/u²+3ab/u)+A'(b²u/v²+3ab/v)
    +B av/u+B'bu/v+ell, γA=γA'=0。

残 v^-2/v^-1 coefficients がそれぞれ A'b²,3A'ab+B'b に限定されることも独立確認した。
γ1 は G3零。γ0 では全 polynomial W0 が producer の四式で実現し f=2ab(Ab−A'a)。
reduced 条件は正確に AA'!=0。ell は c-shift。
actual a→a+αu,b→b+δv は B/B' を行列 determinant −5AA' で全て除ける。
paired torus で A,A'→1、全 inverse は minus shear/非零 scale。
よって全源は一 actual PGL5 orbit。canonical D は c一次で u²v² と残項が互いに素なので絶対既約。
この一軌道の正確な次元や componenthood を主張せず、必要な像上界24だけを使う。

## 5. K4a：全 maximal-pole 再帰と一 actual 軌道

irreducible D は Hhom に正確な u^-4 pole を要求する。
全 H=Σ(k=−4..1)h_k u^k、h_-4..h_-1≤2,h0≤1,h1 constant の28項。
β1 は Yb の residue h−C2(0,1)h_a により全 C2/W で maximal pole 不可。
β!=1 では C2零、v→v+ku,b→b+ka による α→α+(1−β)k で α零へ。
d=1−β!=0 と置く。Yu の三 negative coefficients と Yb top pole を比較し

    h_-4=Aa^n, n=0,1,2, β=(4−3n)/(5−3n), A!=0,
    h_-3=−nAa^(n−1)b,
    h_-2=n(n−1)Aa^(n−2)b²/2−(5−n)Aa^n/(2d)

を得る。途中の linear Euler operators の homogeneous kernel degrees は n−1/3,n−2/3 なので無い。
Ya u^-1 は (h_-1)_b=−(3n²−15n+10)Aa^(n−1)/(6d)。
n0 は a^-1 の非零右辺で polynomial 不可、n1 は全 H の base-degree≤1 となり G3零。
n2/β2 の残 homogeneous kernel は Ba 一方向のみ。全必要 H は

    AΨ+B av/u+ell,
    Ψ=a²v³/u⁴−2abv²/u³+(b²+3a²/2)v/u²−4ab/(3u)。

producer の polynomial W0 を chart Y へ戻してこの Ψ と一致することを手で確認した。
base-gradient f=A(ab²/9−a³/2)、A!=0 で被約。
actual paired shear a→a+tu,b→b+tv は Q を保ち、
(u∂a+v∂b)Ψ=5av/(3u)−4b/3、二次 Taylor 項は linear。
従って B を全て除け、ell は c-shift、A は非零 c-scale で一へ。
global W uniqueness modulo Euler と explicit canonical W があるので full tensor の一 actual orbit が従う。
canonical D は c一次、coefficient u⁴、u0 の残項 a²v³ より gcd一。
一 orbit の像次元高々24。spectral parameter を勝手に有限 samples に置いた証明ではない。

## 6. K4b：全 central 補集合と nonreduced cubic

β!=0 では top H=h/u⁴ に対し Yu が h_a0、Ya が β(5h−3b h_b)=0 を与える。
base-degree≤2 の nonzero h は無い。
β0,C2=δuv では d=1−α!=0 のとき h_a0、
(h_-3)_a=h_b/d、Yb の次係数は5h−(2b+δ/d)h_b=0。highest b-degree で h0。
d0 は Yu から h_b0、Ya の係数が h でやはり零。
残 C2=v² は全 αを保持する。
Ya/Yb top pole から h_-4=A constant!=0、次の全 coefficients は

    h_-3=(4−5α)Ab+B,
    h_-2=5Aa+UAb²+VBb+C,
    U=−(2α−1)(4−5α)/2, V=3−4α,
    h_-1=(6−10α)Aab−(2−α)UAb³/3−αVBb²/2−(3α−2)Cb+4Ba+D。

最後の4Ba は Yb u^-1、全 negative Yu も含む。
degree≤2 を戻すと α=2,1/2,4/5 が強制される。
全 H の base-quadratic 部分は UAb²/u²+((6−10α)Aab+Bα(4α−3)b²/2)/u。
chart Y の base divergence で唯一の Euler gradient gauge へ戻すと

    Wa=K ab+L b², Wb=−K b²/2,
    K=A[((6−4α)(6−10α)−2U)/3],
    L=B α(4α−3)(3−2α)/2。

従って **全 remaining f=b²(Ka/2+Lb/3)**。
A/B/C/D/全 linear H を戻しても nonreduced 又は零なので stated reduced R に矛盾する。
元 higher parameter を零に選んだ因子証明ではない。全 K4b 源が空。

## 7. K0 の四指定部分：linear Casimir/unimodular 補題

`PLANE_CONCURRENT_M3_UNIMODULAR_LINEAR_CASIMIR_EXCLUSION.md` の仮定は
nonzero normal linear L、Q#dL=0、δQ=0 である。
constant normal direction の B_axis/B_other は全 C2 でこれを満たし、
C_tan/C_trans は δQ=t(u∂a+v∂b) なので t0 の全 C2 だけに適用する。
SS/Jordan に linear Casimir があるとは仮定しない。

chart L1 では Z0、全 module は L_Ypi=0。
Q の L-column が零なので δQ=0 はこの chart の div pi0 をそのまま与える。
closed polynomial dual one-form を積分すると pi=Jac(F)、
F=a Q2(r)−b P2(r)+C3(r)、degree d≤3。
全 C2 は C3' として保持され、Q2/P2 の片方が非零なら F nonconstant。
元 module から d(YF)=(divY)dF。
Q2!=0 のとき C(a,b,r)=C(F,b,r) という正確な rational coordinate field を使うと
YF は C(F) に入る。Q2零/P2!=0 は a/b 交換。
coprime A(t),B(t) の Bezout identity を F に代入すれば A(F),B(F) も coprime。
従って polynomial YF は真に P(F)、P polynomial となる。
generic fiber の既約性・Stein factorization・未確認の有理性定理を追加していない。

actualD|L1=±YF、degree五 homogeneous に戻すと D=L^5P(F_d/L^d)。
Y は chart degree≤3、gradF≤2 なので deg P≤floor(5/d)。
d3 なら全 D は L²、d2 なら L で割れる。
d1 は二 linear forms の binary quintic なので C 上可約。
F constant は pi=0、Q_L0 より Q=0、また P=0 なら D=0。
この四指定部分の全 central/W/zero complements が既約 quintic の対象外と確認した。

## 8. 受理した結論と保持する frontier

K211/K31 は全 actual proper factor、K4b は factor 又は nonreduced f で空。
K22/K4a はそれぞれ一 actual PGL5 orbit、従って finite union の像・閉包次元高々24。
このため指定 whole Ore source の **nonzero K4 全枝は G 非優越**。
K0 のうち上の四指定部分以外の非constant-factor whole modules、non-Ore first profiles の別義務、
各 orbit の E所属・ambient componenthood・正確な orbitdimension は本 receipt の結論ではない。
新成分も追加しない。全分類の完成とは呼ばない。
