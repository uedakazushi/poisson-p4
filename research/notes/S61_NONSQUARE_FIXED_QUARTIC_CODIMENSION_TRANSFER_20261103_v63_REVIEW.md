# 固定319四次の余次元二移送：独立概念監査

2026-11-03。主張を受理する。対象は **nuisance零のH8、U0={aFBC≠0}、同じ元135式を使う
一貫したweighted frame** に限る。係数行列・点・体・GBの新計算は行わない。
以下の証明と既受理入力のhash、有限の変数数・行indexだけを検査した。

## 元の相対定理の強さ

TRANSVERSE_CODIMENSIONとその独立監査は、単に一つのcut quotientが長さ65という結果ではない。
補集合がcore次元≤1である開U_anyの各点において、元135式の同じH8相対完成環Sが、
八つのH8底と三つの実在binary座標からなる正則11変数環R上で**自由65**となる。
R→S→Rのrestrictionもあり、同じbinary成分に沿うgeneric local length65が全rankを使う。
この全relative statementを入力として確認した。H11や全S61への拡張ではない。

N=ker(S→R)とすると、K=Frac(R)上のS_Kは長さ65の一つのArtin local factorで尽くされる。
sectionのgeneric residue fieldはKなのでN_Kはその極大idealであり、(N_K)^65=0。
SがR-torsion-freeなのでN^65=0、かつS/N=Rである。従ってS_red=Rが正確に従う。
generic factorが全rankを使うことを省くと、この結論は成立しない。

## 切断を変更する補題

BをH8底の正則完成環とし、R=B[[s1,s2,s3]]と書く。
新cut fのreductionがB上のRの正則座標なら、R'=B[[f1,f2,f3]]からS_redへの合成は同型である。
各N^j/N^(j+1)は有限R-moduleで、そのR'作用はreduced環Rを通る。
N^65=0の有限filtrationによりSはR'上有限となる。

Sは元のR上freeなのでCMである。新R'上でもmaximal idealのradicalは同じであり、
有限moduleとしてdepthはdim R'と等しい。regular local R'にAuslander–Buchsbaumを用いればfreeとなる。
rankは唯一のgeneric componentのintrinsic length65であり、reduced fraction fieldも同じなので65。
従って新cutの閉fiber長も65である。cut間のrow span同一性、横断に沿う形式解の観測、
あるいは単なるnilpotencyだけでfree性を推測していない。

今回のfはraw indexを類推して選ぶものではない。
W^-1yのfree34座標のうちfamily_parameters=[1,2,5]の三つである。
固定graphはW座標でfree blockが単位のT34と零のU34により正確にこのcutに対応する。
元のprimitive binary三方向の、このfへのJacobianは−C/25でありU0上単元である。
必要ならsection valueのB関数を引いてcenterでき、閉core fiberのcutは変わらない。

## 長さ65から真の四次消滅へ

固定normal variablesは七lightと一heavyで、重みは(1,1,1,1,1,2,1,1)。
weight2,3,4の全単項式数は29,91,239である。
二次式は全六行なのでrank≤6、三次は全点でrank≤57という既受理の**上界**もある。
それぞれの下界を欠く集合Z2,Z3はcore次元≤1である。
従ってU_any∩U0からZ2∪Z3を除いた点の実際の低重みHFは正確に

    (1,7,29−6,91−57)=(1,7,23,34), 合計65。

ここで有限jetの係数だけをHilbert関数と呼んでいない。
正重みと一意implicit解により、全34消去座標はdegree≤4の正確な多項式であり、
全残余135式はweight≤6の正確なhomogeneous polynomialsとなる。
positive graded algebraはgraded maximal idealへのlocalizationとcompletionへ単射する。
前者は、非零constantを持つannihilatorと最小重み項を比較すれば分かり、後者はKrull intersectionによる。
したがってweight4の非零元が一つでも残れば、完成環は65次元の低部分に加えてもう一方向を持ち、
新cut長65に反する。ゆえに全weight4 quotientは零である。

weight5/6の生成元はweight4に寄与しない。
したがって結論は全weight4 matrixの239列全射と正確に同値である。
元135式のG6も保たれている。319行証明書のsource label134について七つのlight prolongationが
全てretainedであることを別に検査した。
`t*G6h=ΣK_i F4h` をtの単元除去に用いたり、135式を134式へ戻したりしていない。
354→319は既受理の全91列三次恒等式のlight倍で、削除minor−a^14はU0上単元。
従って**同じframeで定義された**354/319行の像は各点で等しい。
T34/Tnやlazy oracleを混在させた行列にはこの定理を適用しない。

## 全高さ一の網羅と範囲

U0上で固定M4の欠損集合は

    (A³−U_any) ∪ Z2 ∪ Z3

に含まれ、右辺は次元≤1である。従って固定319×239の欠損は余次元≥2。
有理的・幾何的高さ一primeをすべて覆う。
C=0は固定行列の定義域U0の外であり、そこで固定cut全射を新たに主張しない。
旧modp議論で未知因子がC点へ還元される場合を扱ったalternate cutは、U_anyの存在を与えるためだけに用いる。
alternate/fixedのrow-space同一性や全A³への分母除去は不要である。

80本の恒等kernel directionsが**各点で独立**ならsource内のsubbundleとなり、商はrank239のvector bundle。
そこからtargetの同rank bundleへの射のdeterminantが非単元なら高さ一の零点を持つため、今回の定理に反する。
この条件下では全点全射へ進めるが、80本はまだ構成されていない。
abstract freeなkernelやgeneric独立性だけでは不十分である。
codim2に欠損するideal(x,y)への二生成表示が、その反例を与える。

全点239、全U0長さ65、ambient除去、80本構成は今回の結論ではない。
C122へのrank-defect移送は、その別の一貫したframe/rank-identity gateがある場合の条件付きcorollaryに限り、
この監査では新たな係数・frame計算を承認したことにしない。

専用codeとflat receipt：
`code/s61_nonsquare_fixed_quartic_codimension_transfer_review_20261103_v63.py`、
`reports/S61_NONSQUARE_FIXED_QUARTIC_CODIMENSION_TRANSFER_20261103_v63_REVIEW.json`。
