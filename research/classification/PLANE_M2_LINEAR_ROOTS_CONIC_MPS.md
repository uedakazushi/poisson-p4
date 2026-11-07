# m2 flag(1,1)・conic determinantの全sourceはMPS

2026-11-07。主義務G/E。新しい全source幾何証明、計算jobなし。
実Dが既約被約quintic、S=P2がrestriction rank二、zero cubic Rはreduced、位数二、
normal leading quadraticがgeneric squarefreeであるとする。
`PLANE_M2_SPLIT_ROOT_MATRIX_ENTRANCE.md`のflag(1,1)、
G2=h ℓm、ℓ,mはprimitive linear-base/linear-normal forms、hはbase line、
detM=qはirreducible conic、f=hqという全枝を扱う。

**結論。** hが両root-row basepointsを通らなければ全sourceは空。
両方を通る残全sourceはactual weighted SNC curveによりMPS40へ入る。
特にtangent conic+lineのflag(1,1)は全空、transverse conic+lineの同flagは全MPS。
split determinant、他flags、double normal root、m3/nonredR等は別残余である。

## 1. marked rowsを混ぜない全正常形

Mの各rowは二つのindependent base linear formsからなる。
したがってQ1=(ℓ=0)、Q2=(m=0)はrank四quadric threefoldで、
各rowのbasepoint p1,p2がそれぞれ唯一のvertexである。
detMがsmooth conicなのでp1≠p2、両方ともq=0上にある。

linear image W⊂Mat2は三次元hyperplane。
annihilator matrix Kはrank二である。rank一ならdetのrestrictionはreducible quadraticになるからである。
normal右GL2だけでKをidentityへ戻し、base linear coordinatesを選べば

\[
 M=\begin{pmatrix}x&y\\z&-x\end{pmatrix},\quad
 \ell=xu+yv,\quad m=zu-xv,\quad q=-x^2-yz.
\]

これはmarked rootsをmixするleft GL2を使わない正当なnormal formである。
残るhの位置をfullAut(conic)で無断に消さない。
p1=[0:0:1]、p2=[0:1:0]であり、両方を通るhはxのnonzero倍だけ。
元の全nonzero amplitudes、hの全position、normal/base inverseを保持する。

## 2. normal-order-zero Poisson tensorと閉residue

global normal filtrationの全tensorをP=P0+P1+P2とする。
weight零P0もprojective Poissonであり、その実divisorは正確にG2=κhQ1Q2、κ≠0。
従ってP0はgeneric rank四、inverse rational2form ωはclosedで、
三つのreduced components h,Q1,Q2に沿い高々単極を持つ。
generic pointでclosednessはtangential単極項を消すので、ωはgeneric logarithmic formである。

Q1∩Q2は二つのirreducible surfaces SとTのreduced和。
detM≠0ではu=v=0でS、detM=0ではnormal kernel lineのclosureでTとなる。
qがirreducibleなのでTもirreducible。
上のtrace-zero座標ではTは
`(x,y,z,u,v)=(-a st,a s²,-a t²,b s,b t)`というintegral cubic-scroll像であり、
normal kernelのclosureとvertex側の境界も含む。
generic S/T/h-sectionで各intersectionはtransverseで、bad intersectionsはcodimension三以上。
従ってclosed rational residue αi=Res_Qi ωは、S,Tとh-sectionの各componentに沿うsimple polesだけを持つ。
vertexはQi内codimension三なので新しいWeil divisorを作らない。

rank四quadric threefoldの二ruling plane classesを使うと

\[
 H=(1,1),\quad S=(1,0),\quad T=2H-S=(1,2).
\]

hがpiを避ける場合、h∩Qiはirreducible smooth quadric surface、class H。
hがpiを通る場合、h∩Qiは二ruling planes E+F、classes(1,0),(0,1)。

## 3. ruling linesの留数和でgeneric positionsを全排除

まずh∩Q1がirreducibleとする。
二familyの一般ruling planes内で、vertexとbad intersection curvesを避ける一般直線C1,C2≅P1を取る。
α1を制限したrational1formの留数和は零である。
intersection numbers S=(1,0),T=(1,2),H=(1,1)から

\[
 r_S+r_T+r_H=0,\qquad 2r_T+r_H=0.
\]

ゆえにrS=rT=λ、rH=−2λ。
従ってα1−λ dlog(Q2/h²)|Q1は全divisorial residuesが零。
closed simple-pole formのzero residueはgeneric poleを消すので、この差はQ1のsmooth locus上regular。
全一般ruling lines P1への制限はzero。
二ruling-plane familiesの接方向はsmooth Q1の全tangent spaceを張るため、差はrational1formとしてzeroとなる。
resolutionやreflexive-form延長定理を仮定せず、

\[
 \alpha_1=\lambda\,d\log(Q_2/h^2)|_{Q_1}
\]

を得る。

double residuesのantisymmetryによりα2のS,T residuesは両方−λ。
h∩Q2がirreducibleなら同じ証明でα2=−λdlog(Q1/h²)。
h∩Q2がE+Fへsplitする場合も、ruling linesの二関係は
−λ(S+T)+rE E+rF F=0を与え、rE=rF=2λを強制する。
従ってこの場合も同じα2となる。
**少なくとも片方のh-sectionがirreducibleなら十分**である。

global2form

\[
 \omega_{\log}=\lambda\,d\log(Q_1/h^2)\wedge d\log(Q_2/h^2)
\]

は両Qiでωと同じresidueを持つ。
ω−ωlogはQiのdivisorial poleを全て失い、hだけに高々単極を持つclosed2form。
そのresidueはh≅P3上regular1formなのでzero、残るglobal regular2formもP4上zero。
よってω=ωlog、rank≤2となり、rank四P0のinverseであることに反する。

したがってhがp1,p2の少なくとも一つを避ける全positionは空。
qにtangentなhは相異p1,p2を両方通れないので、tangent conic+lineの本flag全源も空。
残るのはh=xだけであり、finite generic sampleからの推論ではない。

## 4. h=xの全higher sourceを一様SNC曲線へ

残る全sourceのP0実divisorは

\[
 D(P_0)=\kappa x(xu+yv)(zu-xv),\qquad\kappa\ne0.
\]

xはreduced degeneracy componentなので、projective P0のPoisson-invariant hyperplaneである。
raw x-rowはx倍のrowとEuler倍の和だが、base Jacobian part f=xqのx-rowは既にxで割れる。
mod xでEuler係数ρはρy=ρz=0を満たすためρ=0。
したがってnormal mixed entriesも含めraw x-row全体がxで割れる。
任意Euler代表のdivergenceを使うshortcutではない。

整数N≥2を選びweights(x,y,z,u,v)=(1,0,0,N,N)とする。
P0のbase Jacobi tensorはweight非負、mixed x-rowはx因子により非負、他mixed及びnormal-normal quadraticも非負。
完全P1三十方向のweightは少なくともN−1、完全P2九方向は少なくとも2N−1で、全て正。
従ってdiagonal GL5から得る

\[
 P_t=g_t^{-1}P(g_tX)g_t^{-T}
\]

は全coefficientsがtの多項式、t≠0で元PとPGL同値、t=0でもPoissonである。

sum(weights)=2N+1に対し、

\[
 x(xu+yv)(zu-xv)=x^2zu^2+xyzuv-x^3uv-x^2yv^2
\]

の唯一の最小relative weight項はxyzuv。
normal order三以上のactualD項はrelative weight≥N−1>0。
ゆえにlimitの実退化因子はκxyzuvで、真のSNC log-symplectic pointである。

このactual irreducible curve全体を含むambient componentはSNC点を含み、
既受理のMPS Theorem1.4の全component scopeによりMPS40の一つである。
元の任意full Pもその閉包へ入る。
SNC点の所属だけからcurve全体を推論せず、T01の一意性も仮定しない。
全root位置・全period/resonance・全P1/P2を保持したpointwise認識である。

## 5. 残余と証明状態

detM irreducible conicのflag(1,1)全源は、empty又は全MPSへ尽くされた。
split detMではQ1∩Q2のresidual Tが分裂し得るため、上の二surface留数関係を流用しない。
そのtriangle/concurrent枝、flags(0,1)/(0,2)、nontriangle(0,0)、double normal quadratic、m3は未完了。
新算術は一切行っていない。producer証明は別作者quartic_reviewと手照合済み、
正式独立receiptは別ノートへ保存する。
