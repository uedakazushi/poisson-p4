# 八HF全源の独立確認：六空、二全族はT36・T38

2026-11-08。PRIMARY E。別作者による純証明の読解と、独立な標準ライブラリ
`Fraction`／疎多項式計算。producer又はCASをimportしていない。

## 結論と範囲

既受理のR=2L+M／m二constant-double全34正常形・全51壁の入口を入力とし、
そこに残る八HF全源を確認した。六源は被約既約quinticを持たず、残る二つの
**全A10族のPGL飽和閉包**は、それぞれ既知の閉成分T36、T38に含まれる。
一般higher parameterだけの帰属ではない。新成分又は新しい分類上の欠落は追加しない。

全first-symbol分類と元の全normalizer逆を今回再計算したとは主張しない。
それらと全像次元上界23は、先行独立
`PLANE_DOUBLE_LINE_M2_CONSTANT_DOUBLE_SEPARATE_REVIEW.md` を明示的な入力とする。
今回、保存済み全30／全9 basisと完全K1から八源の係数を独立に戻し、
FULL9 **inhomogeneous** P2の解・全homogeneous核・元350式を確認した。

## 八源と全parameterの保持

全K1の各列とraw150係数を照合した。保存したP2の同じ非零最大minor、
全free coordinate、全weight二compatibilityを保持したまま、元350式と
signed Euler-Pfaffianを全parameter多項式として再展開した。

ε一の四源では、元式の一つがそれぞれ
`256m²/3, 25m², 16m²/3, 64m²` である。
したがってC上の被約支持ではm零であり、全Dはu因子を持つ。
**全raw P2九には最初からnormal行がない**ことも確認したので、
homogeneous P2だけへ制限した議論ではない。

ε零のk零／k四源では、全Dがそれぞれ
`−2b+p3u+p6v`／`−6b+p4u+p7v` で厳密に割れる。
b係数は非零定数なので、higher parameterの零境界を含めproperな因子である。

ε零のk=−1/2,p=2/5（A）、k=−2,p=−1/2（B）は、
いずれも全十parameterで元350式が恒等零のactual Poisson A10族である。
Aの二homogeneous P2方向を落としていない。Bのc依存項も残している。
P0の非零plane restrictionは固定なので、射影像への写像は全A10で定義される。

## 同じ二anchorの真のsliceと、全次数の検証

whole族から指定した同じ二点を取り、整数scale20／2だけを用いた。
全150 basis・元350式から微分Dを再構成し、Euler25、実GL25＋scalarを加えた
実群51列も定義から再構成した。全行列・完全核・free identity・既存最大minorを照合した。
両点でrank(D,Euler,M)=(101,24,48)、全DM零である。
ここだけからsmooth24又は成分性を推論していない。

|族|normalized torusの座標重み|tensor scalar|入力character範囲|真normal|全graph次数|元350式の全次数|
|---|---|---|---|---|---|---|
|A|(2,−1,3,0,1)|s⁻¹|[−8,6]|−2、一dim|四|八|
|B|(2,2,3,0,1)|s⁻¹|[−6,4]|−2、一dim|三|六|

全150入力・全350出力・全51群列のcharacterを再構成した。
同じD完全核から取ったNはcharacter−2である。全51列を殺す保存functionalは
N上でA:75/16、B:2を取り、Nが群方向でないことを確認した。
同じ101minorの座標方向UはDで独立なので、raw150は
`G48 ⊕ U101 ⊕ CN` の直和である。抽象normal quotientだけで代用していない。

選択101式はD上で可逆なため、真slice上の解は一意なformal graph
`P+tN+Uf(t)` である。f(0)=f′(0)=0。一意性とtorus安定性から、
各入力characterχに現れる項はχ=−2jを満たすtʲだけである。
したがって上表の次数は**全completionそのもの**の上限であり、有限jetからの外挿でない。
保存された全graph係数を再構成し、全補正が同じU内、character−2jにあることを確認した。
元式は二次式なので、全350の全次数零～八／零～六を再展開して零を確認した。

48独立群方向に対応する実一parameter群の積と真102sliceからraw150への写像は、
微分が上の直和ゆえétaleである。Euler加算、GL、tensor scaleは射影Poisson方程式を保つ。
よって局所completionは48滑らかな群方向と全sliceの積になり、
rawではsmooth49、Euler24と射影scale一を除いてsmooth unique24成分になる。
閉軌道又はreductive stabilizerを仮定するLuna論は使わない。

## 実T36・T38認識と閉像の橋

Aはt=144e、Bはt=3eの**同じ全graph**を使う。
e一で、完成済みMPS sourceの全correctionを元JSONから戻し、次を検証した。

|族|係数体|親|固定PGL determinant|完成source parameters (p0,t0,t1,t2,t3)|
|---|---|---|---|---|
|A|Q(√6)|T36|−311040|(300√6,−1800,−1800,−300√6,300√6)|
|B|Q|T38|−288|(24,1,−1,12,−12)|

Q(√6)は`θ²=6`で二項へ還元するFraction演算であり、CASの代数体機能を使っていない。
全M／左右逆、全raw150のforwardとinverseの実Euler恒等式、全10のE∧差、
signed Pfaffianの`D_target=det(M) D_source(M⁻¹X)`を再検証した。
divisorだけの一致をtensor帰属として使っていない。

normalized torusはt→s⁻²tを与えるので、e非零の全graph点は同じ実PGL軌道である。
e一が閉T36／T38上にあるため、非零e全体と、その極限e零のanchorも同じ親に入る。
**この帰属の橋にはsmooth親の成分性すら必要なく、literal全曲線・実同定・閉性で十分**である。
前段の全sliceは追加の局所完全性を独立に認証している。

最後に、既受理の全PGL×A10像上界23と今回のanchor実軌道23を合わせる。
PGL×A10は既約なので像閉包Yも既約であり、23軌道閉包を含みdimY≤23である。
従ってYはそのanchor軌道閉包に等しい。閉PGL不変なT36／T38はこの全Yを含む。
この議論はA10の各点のtransitivity、一般fiberのflatness又は境界での軌道次元一定を要求しない。

## 訂正の扱いと保存

Bの元anchor coneへの全10のE∧P引戻しを再計算して零を確認した。
raw Puv非零からrank二を推論した旧文言は不採用であり、正しくはprojective rank零である。
この誤推論は上の全source／slice／実PGL恒等式に使わない。
旧V1 graphのsyntax failureも数学の根拠に使わず、固定したNEW V2だけを照合した。

一回のgate PASS後の全有限再構成は2.477375秒、RSS101624KiB、全PASS。
checkerは `code/classification_plane_m2_eight_hf_review.py`、reportは
`CLASSIFICATION_PLANE_M2_EIGHT_HF_SEPARATE_REVIEW_20261108.json`。
reportに実行checker全文・SHA、全30入力binding、各scopeの結果を保存した。
共有producer/helper、旧証明書、global台帳、主稿、Gitを変更していない。
