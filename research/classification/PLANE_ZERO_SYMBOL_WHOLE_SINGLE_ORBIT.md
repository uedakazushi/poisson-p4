# rank零plane・零symbol全39源：一つの実軌道への純還元

2026-11-08、PRIMARY E、root producer。独立読解前の新しい純証明。
新しいoperator、rank、minor、kernel、GB、jet又は点の計算は行わない。
以下は旧algebroid上界を強め、全actual sourceを一つの明示tensorへ戻す。
その軌道のambient成分性・既知親はまだ結論しない。

## 1. 全sourceと既受理のalgebroid入口

`PLANE_RANK_ZERO_ZERO_SYMBOL_LIE_ALGEBROID_BOUND.md` と
`PLANE_RANK_ZERO_ZERO_SYMBOL_SEPARATE_REVIEW.md` の全39係数・二chart・全multiplicityの定理を使う。
rank四・被約既約actual quintic、rank零plane S=(z=w=0)、normal-weight零symbol全零を仮定する。
P4\S=Tot(O(1)^3) on P1_(z:w)で、全Π1はA=O(-1)^3上のtransitive Lie algebroid、
全Π2九係数はglobal CE二cocycleとなる。anchor共通零又は恒等零はactualDのH因子又はrank二以下。

K=ker(anchor)は全fiber aff1であり、canonical trace character χ:K→OのkernelはL=O(-5)。
全Kは非零extension e∈H1(O(-5))をprojective scaleで割ったP3から来る。
split classはH0K≠0を与えてK⊂O(-1)^3に反する。
aff1の全derivationはinner、center零なので、元全algebroidはcanonicalにA=Der_Lie(K)。
この旧全source等式と全excluded補集合を保持する。

## 2. 非reduced binary cubic extensionは実sourceに存在しない

Serre dualityによりP H1(O(-5))はbinary cubicのprojective三次元表現である。
PGL2の全軌道は三相異根、二重根＋単根、三重根の三型である。

ここでH0(Der(K))はprojective eのPGL2 stabilizerのLie algebraと同型である。
その理由を明示する。base上identityのLie bundle automorphismは、χを保存し、
Lの定数scaleとH0(L)のtranslationだけからなる。
H0(L)=0、かつe≠0なのでeを固定するconstant scaleも一で、kernelのLie algebraは零。
逆にbaseのprojective infinitesimal symmetryでeのlineを保つものはliftする。
O(-5)をSL2-linearizeし、eの固有scaleの逆でLをscaleすればextension classそのものを固定する。
同じextension classのbundle同型はχを保ち、aff1 bracket
χ(k1)k2−χ(k2)k1を保つのでLie bundle同型となる。
従ってこのbase infinitesimal actionはDer(K)のsectionへliftする。
odd-degree線束のPGL2-linearizationを無断で仮定せず、SL2 liftとLのscaleを使った。

二重根又は三重根のprojective cubic stabilizerは正次元である。
しかし元sourceのA=O(-1)^3にはH0(A)=0。
よってこの二型は存在せず、eは必ず三相異根型である。
このopenは一つのPGL2軌道であり、全normal frameと全underlying Aの同型はactual block GL5で実現する。
実際、固定base上のAut(O(-1)^3)はconstant GL3だけで、normal PGL2は(z,w)のGL2 liftで実現する。
従って全Π1は一つのactual block軌道に入る。

## 3. CE二cohomologyは実際には零

旧独立reviewのlocal CE計算は、cohomology sheafがq=0,1でconstant C、q≥2で零であることを示す。
q=1のgeneratorはcanonical χで、Aut(aff1)が固定する。
Hochschild–Serre/cohomology-sheaf spectral sequenceの唯一の関係は

    d2: H0(P1,Cχ) → H2(P1,C)。

このtransgressionはLのfirst Chern classであり非零である。
これは次の具体的なpushoutでも分かる。
任意のLie derivationはcommutator line Lを保存するので、Der(K)→At(L)がある。
kernel K上のmapはad(k)|L=χ(k) Id。
従ってalgebroid extension 0→K→A→TP1→0をχでpushoutしたものは、
line bundle LのAtiyah extension 0→O→At(L)→TP1→0そのものである。
local horizontal liftのχ-curvatureのČech classはdlog of L-transitionであり、c1(L)に等しい。
deg L=-5なのでこのclassは零でない。
holomorphic global connectionを仮定したのではなく、その存在を妨げるclassを計算した。

従ってd2は両一次元空間の同型で、全CE hypercohomology H1とH2は零となる。
旧reviewと同じく、各term O,O(1)^3,O(2)^3,O(3)の高次cohomologyは零なので、
analytic/algebraic双方のhypercohomologyは同じglobal finite CE complexで計算される。
従って全九変数のCE二cocycleは全てglobal一cochainのd像である。
H1=0とglobal d0=0から、この一cochain空間六次元上のdは単射でもある。

全global一cochainはH0(A*)=H0(O(1)^3)、つまり六actual base-from-normal shearである。
Π0=0なのでそのfinite作用はΠ1を変えず、Π2をdθだけ動かす。
weight二の二次補正を落とす余地は無く、元Π1はbase-linearなのでfinite代入自体がθに一次である。
逆shearの符号反転で全Π2を除ける。
従って全39源のrank四・既約D部分は、前節の一つのΠ1軌道へ入る。

## 4. 三相異根型の明示canonical model

t=z/w、cusp/node baseではなくsmooth P1上で計算する。
L=O(-5)のtransition g=t^-5、aff1 bundle Kのframesを

    e∞=g e0、h∞=h0+b e0、b=t^-1+t^-4

と取る。各chartで[h,e]=e。
このextensionのbinary cubicは両端二項だけの非零cubic、従って三相異根を持つ。

局所Der(K)のframeを(d,H,E)=(∂t,ad_h,ad_e)とする。
反対chart r=1/tでhorizontal frameは

    d∞=-t²d-5tH-(t²b'+5tb)E、H∞=H+bE、E∞=gE。

Der(K)(1)の次の三sectionは両chartでregularである。

    s_a=(1+t³,5t²,4t)、s_b=(t,1,0)、s_c=(t²,4t,3)。

三columnのdetは定数3で、∞chartでも同じ無零のtrivializationとなる。
これでDer(K)=O(-1)^3を実際に得る。
局所微分と[H,E]=Eだけから全bracketは

    [s_a,s_b]=s_a+t²s_b-4t s_c、
    [s_a,s_c]=-2t s_a+4s_b+t²s_c、
    [s_b,s_c]=s_c。

anchorは(1+t³,t,t²)である。
旧全二chart formulaへ戻すと、座標(a,b,c,z,w)のraw homogeneous tensorは

    Pab=aw-4cz、Pac=-az+4bw、Pbc=bz+cw、
    Paz=w²、Paw=-z²、Pbz=zw、Pbw=Pcz=Pzw=0、Pcw=-zw。

これはTot(A*)上の全linear Poisson bracketなのでJacobi零であり、
表示したglobal quadratic tensorのprojective Jacobiもdense openから全P4で零である。
新しいoperator又は選択点を解いて作ったtensorではない。

## 5. actual rank四・既約性と射程

全signed Euler-Pfを手展開すると、scalar conventionを除くactual quinticは

    F=bz⁴-4cz³w+3az²w²-4bzw³+cw⁴。

三base係数のquarticsは3z²w²、z⁴-4zw³、w⁴-4z³wで、gcd一。
従ってFはC[z,w]上のprimitive linear polynomial in(a,b,c)、ゆえに絶対既約。
非零Pfからactual rank四である。
定数方向微分がFを殺すと、a係数3z²w²のz/w微分からnormal二方向は零、
残る三quarticの独立性からbase三方向も零。
従ってFはessential五変数である。
また三base偏微分の共通零はz=w=0だけなので、Sing(D)の支持は唯一のplane Sである。
plane上のtensorは零で、generic leading normal quarticは四相異根を持つ
（例えばa=c=0,b=1ならz(z³-4w³)）。

以上から、対象の零symbol全sourceは**一つの明示actual PGL軌道**の閉包に入る。
全anchor共通零・可換K・split e・nonreduced cubic eの補集合を上で処理した。
全full39源でのactual正逆を保ち、旧Gの≤24をmembershipと誤読していない。
この軌道のambient成分性、既知成分名、相互相違、全rank零plane又は全分類は未判定。
新成分数を加算せず、独立読解後の次の有限義務とする。
