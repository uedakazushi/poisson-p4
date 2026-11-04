# Firstjet B_NEW 上のgeneric Sp4作用はfree

2026-11-04。対象はactual NEW sourceのfirstjet像B_NEWであり、quadratic Poisson成分25全体やB_NEW全点の主張ではない。同じ既存OCTIC_BASE_DOMINANCE index0のQQ beta16と、同じprime32003だけを使う。新しい点・prime・search・GBはない。

## 実sourceと32×15 commutator certificate

元 `Q_FIRSTJET_SPIN5_STRUCTURE_20261103_v63` に保存された二つのquadratic coefficient matrixからC(beta),D(beta)を再評価する。両16×136係数と全beta16分母が32003で単元であることを確認する。QQ Fractionでの全評価と、係数・betaを先にmod32003へ還元した評価を別に作り、両者および既存KrylovのC,Dが一致することを確認する。

sl4の基底はrow-majorの十二offdiagonal Eijと、三つのEii−E33(i=0,1,2)。

    L_beta : sl4→Mat4⊕Mat4, X↦([X,C(beta)],[X,D(beta)])

は32×15 matrixである。全matrixと選択十五行のminorを保存する。minorが非零ならrank15である。rank不足なら別点を試さず停止し、generic failureとは解釈しない。

このminorはZ_(32003)[beta16]上の次数30 polynomial Deltaである。元C,Dの全係数と評価betaがp-integral、Deltaの同じ点での値がp-unitなので、DeltaはQQ上でも非零であり、B_NEWに恒等的には消えない。これはmodp rankからchar0 fiberを無根拠に推定するものではなく、実QQ polynomialの非零minorをp-unit評価で証明する操作である。actual sourceとの結合はOCTIC_BASE_DOMINANCEの既存QQ証明を使う。

## Lie中心化より強い全group中心化

rank15 open上で、C,Dと可換な任意の4×4行列Tを取る。T−tr(T)I/4はsl4のkernelに入るので零。したがって共通GL4 commutantは**全体がscalar**であり、単にLie stabilizerが零という主張ではない。この議論は選択minorが単元な係数環上でも、4を逆にすれば同じ線型代数で成立する。

C,Dは既受理の同時共役Sp4-covariantsである。従ってbetaを固定するg∈Sp4は双方と可換であり、scalarとなる。symplectic条件からg=±Iである。一方Sp4中心−Iは元の三つのcovariant indexによりW16上で−1として作用する。rank15ならbeta≠0なので−Iはbetaを固定しない。従ってgeneric stabilizerはidentityである。char0ではgroup-scheme stabilizerにも同じ結論となる。微分計算だけから有限stabilizerを落としたものではない。

選択Delta chart自体をSp4不変と仮定しない。rank15というintrinsic条件は同時共役で保たれるので、B_NEW内に非空なSp4不変openを得る。B_NEWはactual source像の既約閉包なので、この性質がgeneric点に成立する。

このpacketはfirstjet作用のgeneric freenessだけを証明する。rational orbit section、商函数体index1、B_NEWの有理性、fullPoisson成分25の有理性はここで合成せず、それぞれ別gateに残す。
