# 構成的逆operator recipeの独立監査

2026-11-04。原著 `Q_NEW_FAMILY_CONSTRUCTIVE_INVERSE_RECIPE_20261104_v63` を監査した。
巨大なgeneric adjugate、J14源展開、normを新たに計算しない。

全埋込みデータを凍結元と照合する：G891項、C1のconstant/linear10,314/10,698項、
N79/N41/Fn/J、C/D全32×136係数、Ωとその逆、実rational orbit representative全16β、
従来の全source-five circuitである。Ωの逆は別に乗算確認した。

既に独立受理した同じ100成分Mだけを用い、保存a,bについて20個のdot productを別に確認する。
`M a=e_r`, `M b=e_1` なら、balanced列 `U^j V^(9−j)` の定義から
`r=V9 Σa_j eta^j`, `1=V9 Σb_j eta^j`。
従ってrの表示したquotientは正しく、generic分母はV^-9で非零。
adjugateは転置cofactorの通常の有限行列式式であり、未知の係数を探索する命令ではない。
sの分母はgeneric C1_linear、q/Bの分母は既受理のdense chartで正当化される。
special algebraでVの逆を仮定する議論ではない。

## 固定したJ14の非零性

選ぶ不変式は任意の未検査wordではなく、**決定的に**
`J14={tr(C4),tr(C3 D)}` とする。
既受理の普遍Inv14 basisのsemantic column0がこのbracketであり、
column `[0,1,4]` の非零minor7584に含まれるため、普遍非零である。
全Inv14 dimension3だけから、この特定wordが非零だと仮定してはいない。
NEW定義degree下界26により、このdegree14不変式は実像上でも非零である。

cyclic traceを形式的に微分すると、C4のC-cofactorはC3が4個、C3DのC-cofactorは
`C²D,CDC,DC²` が一つずつ、D-cofactorはC3となる。この別のword計算と
producer内有限関数のASTを照合した。二次C/Dの微分は**全16 ambient β方向**で先に行い、
その後に源を代入する。shapeへのchain ruleだけで代用しない。

実rational Sp4代表βbarをxi=1にした値へこの有限演算を適用すればh14となる。
代替のu²=1600rというrank2 étale algebraとsource circuitも全て埋め込まれている。
既知τが実Sp4同値なので不変式のodd成分は零で、trace/2は同じh14を返す。
全16源の斉次性から `J14=xi7 h14`, `F4=xi2 f4` であるため、最後のxiの比も正しい。

結論はdense open上の**完全に指定された有限有理operator recipe**である。
generic inverseの巨大な展開係数を保存した、境界へregularに延長した、NEW定義式を得た、
全25次元Poisson成分の有理性を決めた、という主張ではない。
新点・新prime・Krylov再生成・J14源評価・source26展開は使わない。
