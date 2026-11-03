# より小さいexact target：delta e0とdelta²

2026-11-03。前のcubic-power条件は正しいまま保持し、さらに小さい二つのmembership targetへ改善する。**一般五点尽くしは、vector membership delta e0∈im Phi、またはscalar membership delta²∈J_gと同値**である。どちらの実membershipも未計算である。

## 1. Canonical kernel

sharp moment constructionのK(a)、formal base R'=K(a)[n42,S1,S2,H]、mを使う。Phi:R'^670→R'^5のtarget basisは1,h,h²,h³,h⁴である。Bが空ならA=E5で、これら五元はK(a)-basisである。したがって

    im Phi=ker(R'^5→E5)=m R'^5

が正確に成立する。E5が五つのK-rational pointへ分裂すると仮定していない。field extensionを含むE5でもvector-space dimensionとprimitive h basisだけで十分である。

## 2. Vector target

e0は実unit1_Aへ写る。Bが空ならdelta∈mよりdelta e0∈im Phiとなる。逆にこのmembershipが成立すればdelta*1_A=0である。deltaはB(a)上の非零因子なのでB(a)=0、faithful field extensionからB=0となる。よって

    B=0 iff delta e0∈im Phi.

deltaのformal parameter次数≤8なので、これが現在の小さいvector targetである。cofactor次数には上界を与えていない。

## 3. Scalar square target

四pivot列をV=(u;D)と書き、

    lambda=(delta,−u*adj(D)), delta=det D

とする。670個のSchur式はlambda*Phiのentriesである。Bが空なら第1節から

    J_g=lambda(im Phi)=m*(lambda_0,...,lambda_4).

delta=lambda_0∈mなのでdelta²∈J_gとなる。逆はJ_g⊂Fitt0(A)⊂Ann(A)とdeltaのB上でのregularityから従う。従って

    B=0 iff delta²∈J_g.

targetのformal次数は≤16となる。旧delta³条件も依然同値だが、小さいtargetが得られた。

## 4. Exact五行toyと範囲

producerは20個の独立変数の5×4 matrix Vと、100本のstandard source columns x_j e_iを用いる。このtoyのimageはmR^5であり、全四つのlambda*V=0、Eulerによるvector target、scalar square targetをinteger polynomial coefficientsで検査する。

これはconditional module proofの小さい検算であり、元S61でdeltaが実際にannihilateすることの標本ではない。一般unit、QQ cofactor、target degree8/16以内でのcofactor boundを証明したとはしない。inhomogeneous generatorsのmembershipにはより高いcofactorが必要な可能性がある。

`code/s61_affine_normal_schur_square_vector_20261103_v63.py` と `reports/S61_AFFINE_NORMAL_SCHUR_SQUARE_VECTOR_20261103_v63/`。先のscalar Fitting、sharp532、独立reviewをpinする。全未知associated/nilpotent structureはcontrolled localizationで保持される。
