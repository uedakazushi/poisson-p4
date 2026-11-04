# Full S61 fiberの独立概念監査：有理16次元

2026-11-04。`Q_NEW_FAMILY_FULL_TENSOR_FIBER_PREFLIGHT_20261104_v63` を受理する。
結論は、actual reduced S61像C16がQQ/C上で有理16次元であること、同じB15上の
generic reduced curveがK=QQ(B15)上で有理であること。25次元projective成分へ戻す
global normalizationは別gateであり、本監査から自動には採用しない。

## 1. βだけからαを推測しない独立のpair argument

元raw value6＋linear24の30×26 matrixはgeneric rank25、唯一kernelのsource部分は
非零、perspective4部分は恒等零である。従ってsource6＋GL4の22列には正確に
一次元kernelがあり、rankは21。これは新たなrank計算でなく、受理済みmatrix/kernel
の部分空間についての厳密な線型代数である。

GL4-saturated rawjet像をYrawとすればdim21。valueが非退化な6次元open上で、
既存の有理Darboux frameはvalueをJへ移す。同じGL4-stable像について、このframeの
逆によるproduct構成があるので、valueJのnormalized RAW pair像Yは **dim15**。
そのデータはprimitive β16とtrace covectorα4である。

t=2のpair sectionは、既受理のβbarと

`alphabar=(0,0,(s+1)/(240B),−32xiBr/3)`

によりE5上で有理である。このpairのSp4-saturationはβへの射影がB15へdenseなので
dim≥15。Yに含まれYは既約15次元だから、このpair sectionのsaturationはYへdense。

すると `Sp4×S5 -> Y -> B15` の合成は既受理のbirational product mapである。
同じ関数体の両端に挟まれるため **Y→B15自体がbirational**。
したがって一般部ではprimitive βのSp4 matchingがraw αも合わせる。
primitive βのinvariantsだけからαの整合性を仮定していない。

このdimension/field-sandwich proofは、著者のEpair=E5/Rosenlicht argumentの独立な補強。
generic orbit閉包の境界にt=2が落ちていないかという懸念も、pair section自身のdensityで
排除する。ODEに接するから有理的に積分できる、という推論を使わない。

## 2. t方向のorbit constancyと有理frame

全source fieldは既受理のF6=F5(t)、kerδ=F5。primitive invariant fieldE5はF5内にあり、
tに依存しない。t=2 sectionは同じ五定数を保つ。B15の既受理birational productと
generic freenessにより、generic βtとβ2は同じorbitで一意のframeを持つ。
§1のbirational pair liftは、そのframeがαtとα2も合わせることを示す。

このことはrawkernelから直接も整合する。値RCRt=Jと
`δC+A C+C At=0` を微分すると

`B=δR R^-1−R A R^-1` は `BJ+JBt=0`、従ってsp4である。
perspective補償は零だから、normalized raw pair全体がそのSp4方向で動く。

physical pushforwardSに対し、primitive βは `(S^-1)*β`、αはS^-tα、
C,DはS^-t C St、S^-t D Stへ移る。従って
`C2P=PCt,D2P=PDt` の一意なintertwiner lineはS^-tのscalar倍。
α2の第4成分は一般部で非零なので、lambda=(Pαt)4/(α2)4により
`Q=P/lambda`, **S=Q^-t** が一意に決まる。幾何的存在と一意性からGalois descentで
F5(t)上有理となり、Cramer chartによる有限recipeとも一致する。
選ぶ15minor chartのtau不変性は不要で、一意なSそのものにequivarianceが成り立つ。

## 3. Full tensorの像を失わない理由

α2=S^-tαtから、5×5 perspective matricesについて

`Persp(α2)*block(S)=block(S)*Persp(αt)`

がexactに成立する。任意16-entry Sにαt=Stα2を入れた全25成分を、独立の小さい
symbolic checkでも照合した。従って著者の新pipelineは、旧full source正規化へ
physical Sp4 matrix `gu^-1*S` を掛けただけで、parameterごとに同じSp4 orbitである。

歴史的S61.normalizeの4 selected-raw gaugeとprimitive Lee=0 gaugeを同一視しない。
perspective U4のactionはJ上で自由、normalizerは一意cross-sectionで
`n(U*T)=n(T)`。Sp4はU4を正規化するので、physical Sp4作用はこのquotientへ降りる。
この既存normalization theoremと上のmatrix identityから、全61像のSp4-saturationが
変わらず、受理済みdense C16が保存される。

六sourceへの二つの削減（n1=1とmodular Gaで評価点y=0）は、既受理のfull tensor有限
同値である。dense18-source mapのrank16をsix-sourceへ単に外挿したわけではない。
scalar cと全評価点の変換も受理済みsix-source gateに含まれる。

## 4. τ下降は全tensorについて成立する

旧全次数parameter同値のeps=a=−1で、cも反転すると
`scaledP_tau=D0_*scaledP`。点もp_tau=D0pなのでtranslationと整合する。
さらに `Rtau D0=H R`, `alphatau=H alpha`, `gu_tau=H gu`。
一意なraw-pair matching frameから `Stau=H S H` が従う。
perspectiveの自然性を使うと、全projective tensorのpipeline中のHは全て相殺する。
これはβだけの偶性ではなく、全次数のfull tensor同値からの結論である。
その後のEuler projectionとQQ定義のS61.normalizeもfield descentを保つ。

旧conic parameterの定義
`t=(rf/sf−1)/(rf/256−1)` とrf↔sfから、τt=t/(t−1)が直接従う。
`z=u(t−2)/t`、`t=2u/(u−z)` は両向き有理でτz=z。
F5=E5(u)、u²=1600r、τu=−uだから **F6^tau=E5(z)**。
これらの有理恒等式は独立に記号検査した。平方根を有理関数と置かない。

operatorの入力も有限に指定できる。E5(z)[u]/(u²−1600r)内でt=2のseedを

`(u,v,c,X,W,B)=(u,−32us/(s²−25rq²),xiBu/40,s/(48u),(q+1)/3072,B)`

とし、保存first integralsを評価して五定数を得る。その定数を保ったまま
t=2u/(u−z)を旧有理逆へ入れる。seedがt=2 relationを満たすことも独立に確認した。
従って「何らかの回収」を仮定したrecipeではない。

## 5. SAME baseの包含とLüroth

下降full tensorをT(a,z)とする。Sp4×SpecE5(z)のfull S61像はC16へdenseである。
physical group変数gを使うとβはrho(g^-1)βbarとなる。既存B15のpullback規約へ揃える
ためgroup変数h=g^-1を使えば、base mapはそのままrho(h)βbarである。
このrational group inversionを明示し、別のbase embeddingを同一と仮定しない。

既受理Sp4×SpecE5とB15のbirational同定の下で、dominanceから

`K=QQ(B15) ⊂ QQ(C16) ⊂ K(z)`

という **同じK上の** field inclusionsを得る。16>15なのでQQ(C16)/Kの超越次数は1。
Lürothの定理によりQQ(C16)=K(w)。KはQQ有理15変数だから、
**actual reduced C16はQQ上有理16次元、従ってC上でも有理16次元**。

z-sourceからgeneric curveへの次数1は証明不要であり、今回その次数は未決。
full source→C16のindex1やw(z)のexpanded formulaを主張しない。
generic reduced fiberの有理性と、全fiberの境界/特異性、full Poisson schemeの
nilpotents・全体の被約性を区別する。projective25へのtransferも別gateである。

## 保存

本checkerは著者を再実行せず、既受理gatesのhash、保存kernel perspective零、
上記の小さい普遍matrix/τ/seed identitiesを確認する。新しい点・prime・source26計算・
full150展開・GBは無い。決定的な新しい監査は§1のnormalized-pair dimensionと
field-sandwich argument、および§3–5のfull dominanceとsame-base fieldの確認である。
