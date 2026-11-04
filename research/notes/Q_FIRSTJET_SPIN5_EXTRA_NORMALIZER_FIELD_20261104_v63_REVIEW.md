# 追加Spin5作用・固定体13回路・stable level10の独立検証

2026-11-04。旧U0-line normalizerのmu4と、元のnormal-line stabilizerの別成分
tauを区別する。元の二biformの全12係数からswap行列とcross4×4行列を再構成し、
tau²=−I、全20項deltaがtauで反転することを確認する。

五つの独立係数p,q,r,s,t上で、quadratic algebra kappa²=ResVのpair算術により
critical-root frameを別に計算する。変換後の両quadraticのmixed係数は零、
detA=1、pure-square係数行列detM=−kappaを確認する。c'²=−kappaなので
B=c'M^-1はSL2、交換後Vの係数はc*c'²*w_i。この事実と旧六不変式の定義から
保存された六imageのscalar・符号・冪が従う。

六変数上では小さいResVとdeltaだけを、別のquadratic pair算術で厳密に下降
する。Norm(J0)は `(I0*I1*I3−I4)²−I0*(I1*I2−I5)²`。I0、I1、ResV、
Norm(J0)、deltaは一般係数体で非零であり、全inverseの一般域が正しい。
それらが全てのspecializationで単元だとは言わない。

1169-node回路のrank4 etale算術・critical frame・六W係数・六image・13固定式
の全recipeを照合する。z²=I0、kappa²=ResVの拡大が体あるいは独立であるとは
仮定しない。basis(1,z,kappa,zkappa)の定数係数は常にtrace/4である。

zの符号を変えた旧frameは、受理済みのsigmaの後にv1を1へ戻すtorus操作を
したものと正確に一致する。この小さい有理恒等式も確認する。kappaの符号は
新frameのcritical根の順序を変えるだけで、二つの正規化は受理済みの**全**
line normalizerで結ばれる。従って六不変式の像は両符号に不変であり、rank4
algebraのbase Kに属する。大きい三成分を展開して零を確かめたとは主張しない。
ここでroot-choiceの完全性は既受理のexact normalizer分類に依存する。

tau²が連結群の中心なのでK上では二乗1。delta非零でtau(delta)=−deltaなので
非自明なinvolutionである。d=delta²、S_i=(I_i+tauI_i)/2、
A_i=(I_i−tauI_i)/(2delta)の13式をF0で生成するとK=F0(delta)、delta²∈F0。
非自明性より[K:F0]=2、ゆえにF0=K^tau。前のassociated-bundle同定により、
これは元のW16/Sp4不変関数体そのものを与える。13生成元の独立性や純有理性
は主張しない。

generic freenessも成立する。最大minorから作るnormalはlineではなく共変な
**vector**で、spinor側のdeterminantは1。従ってgeneric stabilizerはconnected
Spin4に入り、Uのexact stabilizer Tに絞られ、Vのweight1係数でt=1となる。
中心−Iも一般betaを固定しない。必要なdeltaとresultantの非零条件は自由係数
上の非零polynomialで同時に満たせる。新しい点の試行は不要である。

元の群はSO5ではなくSp4である。任意の標数0の体でSp4-torsorは非退化交代
形式を記述し、symplectic Gram-Schmidtで標準形になるためH1(F,Sp4)=1。
generic torsorは自明。Sp4のrational Bruhat cellの次元4+2+4=10より、固定体に
10変数を付加した体は純超越となる。stable rationality level上界10のみであり、
純有理性・最小level・多項式不変式環・NEW方程式を結論しない。

旧6変数image展開の60秒timeoutとsourceは保持し、再実行しない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_firstjet_extra_normalizer_review_20261104_v63.py --output work/extra-normalizer-independent
```

60秒・2GB。小さいframeと有限回路recipeだけを検証し、巨大な自己合成や展開は無い。
