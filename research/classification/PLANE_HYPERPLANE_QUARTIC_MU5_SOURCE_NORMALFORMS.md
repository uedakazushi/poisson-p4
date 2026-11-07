# μ5全五有限源とA3境界：全sourceと既知103への正則lift

2026-11-07、PRIMARY E。producer exact source batch及び手橋、独立確認前。
全marked geometryは `PLANE_HYPERPLANE_QUARTIC_MU5_MARKED_LIST.md`。
全五finite candidateとA3+A2 k−1controlを一回のsource jobで検証した。

## 全source batch

六case全てで、全28×18 gradient-syzygy行列はrank十六、全kernel二次元。
actual H→h GL3と逆、全quadraticSaito basis、全六磁気sourceを保存した。
strict十・projective三百五十Jacobi及び実Pfaffianは全magnetic coefficient上で恒等的に確認。
全六でcofactorはrank二平方零、hはdualimage外、affineA rank二・constant∈imA、正確にtype6nil22-image。
全6×6磁気translation行列はrank五である。
report `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU5_SOURCES_20261107.json`、
実測1.799秒、最大RSS71,516KiB。
旧全FunctionDef ASTを保存するisolated modeであり、D5旧ambientを再計算していない。

## 全五有限sourceの23軌道閉包

rank五translationの全polynomial normalizerはmu6全source noteと同じ定数5minor公式で構成できる。
未知magnetic coefficientを反転せず、全sourceをR=(I/Ij)*mon_jへ移す。
I!=0をfiber homothety＋tensor scaleで一にし、I0はR=t*mon_jの正則端点で戻す。
各全marked pairのbase automorphismは有限である。
A4+A1と2A2+A1ではnormalization上のcusp/node inverse三点以上を保つから有限。
D5では異なる二branchと二接線を保つmatrixの相対invariance比較が全scalarを強制する。
A5では二branchを固定するr-scaleがy sectionのs³t/st³比を保つ条件で有限となる。
従ってintrinsic fiber kernelとanchor equalityを使うmu6と同じargumentで、
nonzero-I pointのPGL5 stabilizerはtranslation kernelの一次元だけ、軌道は正確に23次元。

canonical RはA5 smoothhyperflex、A4+A1 smoothhyperflex、D5 smoothhyperflex、2A2+A1 smoothhyperflexでh²、
A4+A1 node-inflexionだけhyである。これは保存cokernel vectorの最初の非零係数による。
全五source閉包はそれぞれこの一23軌道閉包に等しい。
全周囲成分性又は既知親への帰属はまだ確認していない。

## A3+A2 k=−1の全sourceを既知103へliftする

whole marked familyのpolynomial quarticを明記する。

```
V_k=hw-y²-4hy+(k-3)h²,
F_k=V_k²-(k-1)y(y+2h)V_k+(k-1)²h²y²,
k!=0,1, H=h.
```

normalizationは元の全familyであり、k−1では正確に

```
F_(-1)=h²(w-2y-4h)²-y³(y+4h).
```

generic kはa=(1-k)²/[4(1+k)²]により既知cuspidal24のF3(a) marked familyへ実Mobius/base変更で移る。
fullSaito basisのconstantGL2 frameと全六magnetic pullbackも移るので、generic全sourceは既知103に属す。
k−1でaを割る式を使わず、次の正則な全source構成を用いる。

g_k=hF_kに対する全28×18 gradient matrix A(k)はkのpolynomial matrix。
generic rankは既知sourceの十六なので全17minorは恒等零である。
k−1では今回もrank十六で、保存したrows
`[0,1,2,3,7,8,9,13,14,15,18,19,20,22,23,25]`、
columns `[0,...,15]` の16minor値は **−3913788948480!=0**。
従ってこのminorが単元となる近傍でrankは正確に十六、全quadratic kernelは正則rank二bundleとなる。
free columnsを二本取り、pivot block inverseでkernel basisを作れば、その係数はk−1で正則。
任意のspecialSaito frameに定数GL2で合わせてから、このkernel basisでliftする。
Saito determinant ratioはspecialで非零なので近傍で単元。
h-Casimir anchorとlinear commutator coefficientsも正則にliftし、全六magnetic Rをそのまま自由に保つ。
よって任意のspecial fullsource pointを端点に持つ実Poisson familyが構成される。

generic fiberが既知103閉成分に属することと閉性から、k−1の**全六magnetic source**も同成分へ入る。
これはdivisorのboundary membershipから別sourceを逆推論したものではなく、
全sourceをrank二kernel bundleから実際に正則liftしたargumentである。
従って全A3+A2 common-tangent branchはgeneric/exceptionalを含め既知103へ被覆される。

## 残る範囲

全mu五marked geometryは、旧A0=D5 cusp-tangent、既知103=A3+A2 common-tangent、
及び上記五finite sourceへ還元された。
五つのfullambient成分性/既知帰属と、この全marked/source証明の独立reviewは残る。
mu六の全source・残余、stable二源、nilshort等を本resultで閉じたとは主張しない。
