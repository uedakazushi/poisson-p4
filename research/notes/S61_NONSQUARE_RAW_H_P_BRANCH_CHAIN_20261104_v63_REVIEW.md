# Raw Hのp分枝五packetの独立QQ監査

2026-11-04。authorコードをimport/rerunせず、変数順をsource(d,e,a)、target(d,z)へ逆転した独自Fraction Laurent辞書で検算する。
対象はP_BRANCH_PROFILES、G1_UNIT_PAIR、THETA_UNIT_PROFILE、THETA_ONE_DEFECT、THETA_RATIONAL_COORDINATESの凍結五packetのみ。

## Sourceとscheme identity

全primary input hashとcertificate hashをstreamで確認する。
raw checkpointの14x5全70係数を、既独立受理GENUINE14 bridgeと直接比較し、row/column scalingが既知FULLU単元だけを使うことを確認する。
pの元係数、750000p=g1*Theta、B²=Theta+C*g1*g2を全QQ係数で検算する。
したがってFULLU上の二分枝はcomaximalな正確なscheme分解であり、交点を落としていない。
全56 source quotient/remainder identityとmonic degree boundを照合した。

## 二つのknown-unit block

g1分枝のraw rows5,8の全四係数と三middle列のliteral零、B−50T/dのquotient、
det=T^6/(750d^4)と既知B単元による代表、左右adjugate identityを独立に検算する。
従ってこの分枝ではheavy rank=2+rank(B7x3)。B全点rank3とは主張しない。

Theta分枝のraw row8について、f=g1*g2²/37500、inverse37500C²g1/B⁴を元係数のquotientまで含めて検算する。
三middle列が元から零なので、実7x4 Schur blockの全28係数を元raw sourceから再構築し、保存quotient/remainderへ掛け戻した。
第一三列は元B7x3のTheta制限そのもの、最後はB⁴*last-Schurである。
そのためTheta上ではheavy rank=1+rank(7x4)であるが、五packetだけからその7x4のrankは決めない。

固定genericB minorと全21Cramer identityも元raw係数から独立に照合した。
U,V,Deltaはglobal unitとして扱わない。
raw11と7の固定一関係、およびその一つのextra-column defectがThetaで零となるidentityを確認した。
これはその一関係だけのnegative結果であり、残りのdefectを検査したとは言わない。

## Rational Theta schemeの全localizer

N²−20g1g2²=d²Theta、z=N/g2と既知B分母のregular代表を直接検算する。
target a,e,B,C,g1,g2,L,K,Fnum,Jnumを公式から独立に再構築し、全forward imageと六reverse compositionを全係数で確認する。
source g1,g2の単元性、z,z+20,L,K,Fnum,Jnumの全逆元witnessを、元Thetaの多項式identityとして確認する。
したがってlocalized QQ[z,d]とのring isomorphismはscheme全体で有効で、reducednessやirreducibilityを仮定しない。

この監査はheavy/pure source transportの範囲を維持する。
full cubic57、quartic239、Theta全heavy rank、未知minorのglobal unit性や新Poisson成分について結論しない。
新しい点・prime・ranksearch・GBは行わない。

Code `code/s61_nonsquare_p_branch_chain_review_20261104_v63.py`。
Report/folder `S61_NONSQUARE_RAW_H_P_BRANCH_CHAIN_20261104_v63_REVIEW`。
