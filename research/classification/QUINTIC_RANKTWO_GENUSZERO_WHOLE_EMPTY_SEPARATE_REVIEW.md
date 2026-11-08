# Degree五・種数零・rank二の全空性：113／023／014と全marksの独立監査

2026-11-08、PRIMARY E、別作者 `exceptional_coverage`。**純幾何証明 PASS**。
`QUINTIC_RANKTWO_GENUSZERO_113_EXCLUSION.md`（d6a20a2e…）と
`QUINTIC_RANKTWO_GENUSZERO_CONE_EXCLUSION.md`（ea0a357c…）を、
既受理の全source還元と共に読解し、以下の全場合を別に再導出した。
新しいCAS、数値点、operator、rank、minor、kernel、GB又はprofileは無い。

受理する命題は、被約既約actual quintic・一般rank四Poisson構造において、
degree五 integral非正規支持のrestriction rank二、正規化sectional genus零である全源が空、というものである。
全射影、plane上rank零／二の旧場合分け、全marksと重複、small exceptional curveを保持する。
他degree、genus一、rank零、全E又は新component数は本稿の結論ではない。

## 1. 全sourceから本当に得られる入力

旧 `RANK_TWO_QUINTIC_SURFACE_GENUS_ZERO_SEPARATE_REVIEW.md` は次を全源で証明している。
length式はdelta_S=1と唯一のdegree一残余plane P、delta_P=1を強制する。
Sのgeneric cusp、non-split node、005、122、全C_P=2F、014の2Qは既に空。
残る三型は

```
113: CP=Q=L-3F, Q=P(O(1)^2)=P1xP1;
023: CP=Q=L-3F, Q=P(O+O(2))=F2;
014: CP=Q+Fb, Q=L-4F, Q=P(O+O(1))=F1.
```

S側は二つの異なるprime `Gamma1,Gamma2 in |L|` にsplitし、各々Sへbirational。
旧留数計算はP側のresidueを零とし、全closed一形式は

`alpha=lambda*dlog(g1/g2), lambda!=0, Gamma_i=(g_i)`

である。これはlog ansatzの部分源ではない。
この式と同じ全residuesを持つ二形式の差は全divisorで正則、smooth scroll resolution上のHartogsと
H0(Omega1)=0から零となる。

generic nodeでalphaのresidueは、同じdualizing generatorのunit倍のP上tangential bracketに等しい。
ゆえに旧residue零は実際に `pi|P=0` を意味する。
113/023のQは一primeでdegree二、014のQ/Fbは二primeでdegree一ずつPへ写る。
CPは全場合でreducedなので、delta一のunibranch cuspではなくordinary nodeである。
一primeの場合のnodeはfunction field extension後にsplitしてよい。
係数体をPoisson中心に置かず、全tangential bracketを保つ。

## 2. Actual nodeの非零wedge条件

generic Pの有限分離拡大・完成後、D=(uv=0)とする。
actual Pfaffianはuvのunit倍なので、inverse symplectic form sigmaはDに高々simple poleを持つ。
closednessからsigmaはSNC logarithmic二形式である。
これは二位極係数を直接見ても、normal differentialを含まないnumeratorが対応するu又はvで割れることから従う。

log coframe `(du/u,dv/v,da,db)` においてsigmaを

```
A du/u wedge dv/v
+du/u wedge(B da+C db)+dv/v wedge(D da+E db)+F da wedge db
```

と書く。逆行列のcofactorからP上tangential bivectorは非零unitを除いてA|Pなので、
pi|P=0はA|P=0。
一方log matrixのPfaffian `AF-BE+CD` は、uv/Pf(pi)のunit倍だからunitである。
従って二branchの正則restriction
`theta+=(B da+C db)|P`、`theta-=(D da+E db)|P` は

`theta+ wedge theta- != 0`。

branch residueに含まれるA*dlog(other normal)も、A|P=0によりgeneric Pで正則になる。
同じactual anticanonical/dualizing normalizationから来るalphaの二restrictionは、
これらthetaと同じ非零constant規約で一致する。
normal unit、branch交換、finite separable field extensionは非零wedgeを変えない。
以後の全矛盾はこのgeneric必要条件への矛盾であり、全域clean crossingを仮定しない。

## 3. 同じdualizing generatorの反deck residue

Omega_Yをactual quinticの同じdualizing generatorのrational pullbackとする。
full conductorは有限双対性から正確にその負divisorである。
113/023では `div(Omega_Y)=-Q-Gamma1-Gamma2`。
Qのresidue betaは

`div(beta)=-(Gamma1+Gamma2)|Q`。

局所的にOmega_Yはunitをq*g1*g2で割ったvolumeであり、Qはsmoothなdivisorなので、
residueはこの全scheme multiplicityを正確に保つ。
GammaとQのtangency、反復restricted curve、small exceptional curveを消していない。
generic nodeの同じdualizing formの二branch residuesは反対符号なので、
Q->Pのfunction-field deck involution iotaに対し `iota*beta=-beta`。
ambient unit自身がdeck固定であるという誤った仮定は不要。
従って `div((g1|Q)*(g2|Q))` はiota不変。
各restrictionは非零であり、零なら異なるdivisor classのQをprime Gamma_iが含んでしまう。

## 4. 113：全bilinear-factor補題

QのLはO(1,1)で、Q->Pはsmooth quadricのoff-quadric点からの全double projection。
basepoint-free性が中心をquadric外に置く。
全外点はAut(Q)で移り、中心をinvertible skew matrixにすればprojectionはSegre matrixのsymmetrization。
deckは二つのP1 factorsの交換になる。
一般中心だけを選んだ議論ではない。

a=g1|Q、b=g2|Qは全非零bilinear sections。
div(ab)のfactor-swap不変性から、`iota*dlog(a/b)=+dlog(a/b)` 又はその負となる。
全factor場合を次のように確認した。

* 両方irreducibleならiotaは二つの(1,1) primeを固定又は交換する。
* 一方だけirreducibleなら、それを固定し、残る二ruling factorsの和も固定する。
* 両方reducibleならa=l(s)m(t)、b=n(s)o(t)。不変性はbinary quadratics lnとmoの比例。
  l,nが独立なら(m,o)は(l,n)又は(n,l)の定数倍で、ratioはR(s)R(t)又はR(s)/R(t)。
  l,nが比例ならm,oもその同じfactorに比例し、ratioはconstant。

共通ruling、反復factor、比例restricted sectionsもこの場合分けに入る。
proper smooth Q上でdivisor零のrational functionは定数なので、未知のnonconstant unitは残らない。
したがってalpha_Qとdeck pullbackは比例し、wedge零。
§2のactual-node条件に反し、全113は空。

## 5. 023：全quadric-cone sectionとruling場合

Q=F2でL_Q=C0+2f、K_Q=-2L_Q。
|L_Q|はnegative sectionをvertexへ縮めるquadric cone mapで、Q->Pはconeの外点からのdouble projection。
有限性を使うのはX->Dであり、negative sectionを有限に写すとは主張しない。
deckはconeのvertexを保ち、rulingを保つためminimal resolution Qへ延長し、base P1のautomorphismを誘導する。

全L_Q sectionは `a_i*x+B_i(s,t)*y`（a_i定数、B_i全binary quadratic）。
次の三場合が全てである。

* a1*a2非零：各divisorはvertexを避けるirreducible conic。deckは二つを固定又は交換し、ratioのdlogは±。
* 一方だけa_i零：そのdivisorはC0と二ruling fibers、他方は唯一のhorizontal irreducible conic。
  rulingを保つdeckはhorizontal conicを固定し、残りのdivisorも固定するのでratioのdlogは不変。
* 両方a_i零：共通yを消したratioはB1/B2で、元一形式もdeck pullbackもruling base一変数から来る。

二binary rootsの反復、共通fiber、比例marks、二conicの一致を全て含む。
全場合でalpha_Q wedge iota*alpha_Qは零、§2に反する。
smooth quadricのfactor-swapをconeに無条件流用したのではない。
よって全023は空。

## 6. 014：二sheetが持つ同じ第三lineと全二line比

Q=F1、L_Q=C0+fはplane blowdownのpullback。
Q->Pはこのblowdown、Fb->Pはlinear同型で、二つのnode sheetsを与える。
Y上で同じcurve Q intersect Fbは、両方のmapで同じline ellへ写る。
Qではruling fiber、Fbではdegree一lineであり、二つの独立な第三lineを選ぶ自由は無い。

Gamma_iのrestrictionは両sheetで非零。
Q上ではplane line A_iのtotal transform、Fb上ではline B_iである。
A_iがblowdown pointを通ればtotal transformにC0が含まれるが、planeのdivisorial imageは同じA_iだけ。

`div(Omega_Y)=-Q-Fb-Gamma1-Gamma2` から二residuesのdivisorsは
Q上で `-(Fb+Gamma1+Gamma2)|Q`、Fb上で `-(Q+Gamma1+Gamma2)|Fb`。
同じP上のrational二形式としては、generic nodeのtrace-zeroにより互いに反対。
従ってplane divisorsが一致し、

`A1+A2+ell=B1+B2+ell`、従って `A1+A2=B1+B2`。

blowdown exceptional curveはP上codimension二なので新しいplane divisorを作らない。
ellがblowdown pointを通ることとも整合し、pullback canonical divisorのexceptional補正はその点に吸収される。
A_i又はB_iがellに一致しても同じmultiplicityのcycleとしてcancelできる。
planeのUFDからunordered line pairsが一致し、
`lambda*dlog(A1/A2)` と `lambda*dlog(B1/B2)` は±の関係。
比例lineなら両方零であり、やはりwedge零。
§2に反して全014は空となる。
arbitrary plane matchingやgeneric projection閉包の仮定を使わない。

## 7. Schubert補集合と全aggregate

全finite実projectionでCPのimageはplane P全体なので、restriction V->H0(CP,L)の像は正確に三次元。
kernelはV intersect K、K=H0(L-CP)であり、planeの二linear equationsを含む。
したがってdim(V intersect K)=二である。
旧supersource上界に含めたintersection三／四は、実finite imageでは空。
small mapが曲線を縮めてもCPのimage dimension二は変わらず、この結論は全三型で有効。

旧全源が除いた005、122、全2F、014/2Q、generic cusp S、non-split node Sを戻した後、
残る113/023/014の全exact-intersection二sourceは§4–6で全て空。
全marks、全R/交差のmultiplicities、全射影中心・frame、small/non-normal境界を保持した
**degree五・genus零・rank二全分岐の空性**として受理する。
global、main、Gitと旧producer/review bytesは変更しない。
