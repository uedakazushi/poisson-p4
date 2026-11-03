# 新一次像hypersurfaceの方程式次数は少なくとも8

2026-11-03。新しい15次元primitive beta像Bの方程式はまだ同定していない。
ここでは完全な低次数Sp4不変式の計算から**defining equationの次数≥8**を証明し、
曲線判別式と二つのmoment判別式という自然な候補を除外する。
標本で多項式の恒等零を推測していない。

## defining equationは斉次Sp4不変式

Bは実有理sourceから得た既約hypersurfaceで、全Sp4とS61 weighted scalar作用で安定である。
後者はprimitive betaに一様scalarとして作用するので、Bの素idealの生成元fは斉次に取れる。
Sp4は連結半単純群でcharacterを持たず、安定な主idealの生成元はSp4不変となる。
Sp4中心−Iはlinear二形式に三つのcovector indexを持つbetaへ−1で作用するため、odd次数の不変式はない。

## 全weight-zero空間と全Lie制約

閉かつprimitiveなlinear二形式の24係数を8線形制約で切り、16次元有理基底を作る。
Cartan二元のweightで零になる全単項式を列とし、全十Lie生成元による微分を行に並べる。
char0の連結群なので、これら全Lie制約のkernelが正確に同次数の全不変式空間である。

degree2ではweight-zero12列、degree4では168列、degree6では1448列を全て含め、
QQ rrefから不変式次元をそれぞれ**0,1,0**と得る。
degree4の一意な不変式I4の全係数を保存し、全Lie微分が零であることも再検査した。
実dense mapから得た記録QQ beta上で

    I4=−6287207/1638400≠0。

従ってBはI4の零集合ではない。
odd次数がなく、degree2,6もなく、degree4候補も消えないので、fの次数は少なくとも8である。
複数点で低rankを見ただけの上界証明ではなく、全有限QQ線形系の完全計算である。

## 曲線とmomentの候補除外

同じ実betaはintrinsic elliptic-quintic curveのDelta≠0を満たす。
従ってBはdegree24曲線判別式hypersurfaceではなく、その一般点の曲線は滑らかである。

また三つのtensor indexをsymplectic Jで縮約してW16の非退化な交代形式omegaを作り、
全Lie生成元がomegaを保つことを厳密に確認する。
moment mapを

    mu_A(b)=(1/2)*b^T*omega*rho(A)*b

と定め、trace pairingでsp4行列へ戻す。
これはSp4-equivariantなので、det(mu)とTr(mu²)²−16det(mu)はdegree8不変式である。
記録betaで両方非零となった。
従ってBはこれら二つの零hypersurfaceの成分にも含まれない。
一般betaのmomentはregular semisimpleであり、この認識を特異moment型と混同しない。

## 保存と残件

producer `code/q_new_family_hypersurface_recognition_20261103_v63.py`、
quartic `code/data/q_new_family_quartic_invariant_20261103_v63.json`、
receipt `reports/Q_NEW_FAMILY_HYPERSURFACE_RECOGNITION_20261103_v63/`。
全表現・weight・有限系rank・I4・moment値を保存した。
degree8以上の別の不変式や、既知成分の一次像との同一性は未解決である。
この結果はhypersurface equationそのものの構成や一般fiberの全scheme記述ではない。
