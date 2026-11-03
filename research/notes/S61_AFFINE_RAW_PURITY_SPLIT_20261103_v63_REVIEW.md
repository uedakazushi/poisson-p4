# raw compact Yの非純性と孤立無限遠のsheaf分離：独立検証

凍結済み `S61_AFFINE_RAW_PURITY_SPLIT_20261103_v63` を、一般体 (K=\mathbf C(W_{16}))、または必要な非零分母を避ける十分小さいbeta開集合上の結論として受理する。raw Yは純3次元ではない。選択した単純孤立無限遠は正確なsheaf saturationと補成分colonで分離できるが、残余の空性や全w-torsionの除去を意味しない。

## 元134式への独立な照合

元の各式の45個のq変数に関する次数を直接数え、全てを次数2へhomogenizer wで斉次化した。beta係数を保存されたQQ点へ代入し、別の係数・occurrence微分で全134×46 Jacobianを作った。最高q項の全134値は零であり、保存された45行は追加の135番目の式を一切使わない。

固定した非零q座標を除く44列とw列の45×45 determinantは、独立なQQ計算で

\[
-124102513209883356037120000000000000000/59049
\]

に一致した。全46列のJacobianがprojective radial vectorを殺すEuler関係も確認した。従ってglobal homogeneous degree45 minor Δが非零なら、削除したq座標も非零である。実際、その座標が零なら、残る45列に非零radial kernelが残りdeterminantが零となる。よってこの開は正しいprojective chartの全rank45相対Jacobian開である。

古いaxis族も、単に保存された16×16行列式を再計算するだけではなく、物理的なquadratic tensorからS61正規化して再構成した。明示的な10変数rational sourceで全135最高q式が恒等的に零となること、QQ sampleが保存61座標と一致することを確認した。primitive betaの三つのtensor slotに作用するSp4微分を別実装し、source10方向と合わせた全16×20微分を再構成した。保存された16 minorは厳密に569である。

従って既存のSp4-equivarianceと合わせ、このsourceの像はbetaを優越する。45 minor非零の開と交わるsourceは稠密であり、一般betaでも被約孤立無限遠点を与える。beta全射を使わずに、相対Jacobianだけからtotal spaceのflatnessを推論してはいない。

producerを新しい出力先で一度実行し、証明書のbyte一致も確認した。旧sourceと関連証拠のhashは維持した。

## normal-direction compact incidenceと同じsourceである理由

この議論の射影空間は元134式を斉次化したactual raw modelである。normal=0の既知三変数planeについて、保存された五つの最高二次式に三変数を掛けた15×10 cubic matrixを独立に再構成し、選択した10×10 blockと保存QQ inverseの積が単位行列であることを確認した。従って全10 cubic monomialがそのidealに入り、一般betaのこのplaneにはprojective infinity点がない。

よって無限遠の点ではnormalは非零であり、そのprojective方向は一意である。normal-direction compact incidenceとの比較はここで同型になり、孤立axis成分と既存の135式の3次元raw boundary族の両方を失わない。135式の零点は134式の零点に含まれるので、134側にも3次元以上の部分がある。従って一般betaで0次元の成分と3次元以上の部分が共存し、raw Yの純3次元性は否定される。134側の最大次元を正確に3と主張する必要はない。

全YでΔを引き戻した開を使う場合にも、単純点上に余計なnormal-direction fiberを作らない。normal=0の点で134式をF5とG129へ可逆な行操作で分けると、Gのplane接微分は零である。全相対rank45とplaneの三変数からGのnormal微分はrank42となり、そのkernel内の非零方向は存在しない。これは元の単純点をexceptional方向へ重複させる可能性も排除する。

## clopen選択と二つのsheaf ideal

ここからはgeneric体K上で議論する。Δ非零の点はprojective相対tangentが零であるから、actual零点schemeでは被約な孤立点である。有限型の0次元scheme (D(\Delta)\cap Y) は有限かつ被約で、特性零ではfinite étaleである。Yはproperなので、その有限部分はYの中でも閉じている。従って既に開であるこの部分はclopenとなる。

このfinite étale algebra上でw=0を課すと、各field factorでwは零または単位なので、さらにclopenな部分 (Y_{\infty}\) を選ぶ。Yのideal sheafをIとし、

\[
I_{\infty}=(I+(w)):\Delta^{\infty},\qquad
I_{\mathrm{rest}}=I:I_{\infty}
\]

とする。**最初は単一colonではなくΔ-saturationである。** 単一colonではΔがnilpotentな別の点を残す可能性があり、同じ結論は出ない。Δがline bundleのsectionであることは、局所trivializationの単位倍がsaturationを変えないので問題にならない。irrelevant saturationもsheaf化に含める。

最初の式は (D(\Delta)\) の選択部分のscheme-theoretic closureである。選択部分は既にclopenなので、その正確なidealを与え、余計なembedded closureを加えない。\(\mathcal O_Y=\mathcal O_{Y_\infty}\times\mathcal O_{Y_{\mathrm{rest}}}\) と書けば、\(I_\infty/I\) とそのannihilatorは相補なidempotent idealsである。従って補成分には単一colonでよく、sheafとして

\[
I=I_\infty\cap I_{\mathrm{rest}},\qquad
I_\infty+I_{\mathrm{rest}}=\mathcal O_{\mathrm{ambient}}
\]

を得る。(D(w)) では最初のidealは単位idealとなるため、(I_{\mathrm{rest}}=I) が厳密に成り立つ。未知affine点、正次元成分、非被約構造を消していない。

この記述はfinite étale factorの次数、展開したcolon基底、全w飽和を計算したものではない。残余が純次元であること、他の孤立無限遠を全て選んだこと、一般affine fiberの空性も主張しない。generic体での分解を、全betaや任意の特殊fiberへ無条件に拡張しない。

独立codeは `code/s61_affine_raw_purity_split_review_20261103_v63.py`、receiptは `reports/S61_AFFINE_RAW_PURITY_SPLIT_20261103_v63_REVIEW.json` である。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/s61_affine_raw_purity_split_review_20261103_v63.py \
  --replay work/raw-purity-producer-replay \
  --output work/raw-purity-independent-review.json
```
