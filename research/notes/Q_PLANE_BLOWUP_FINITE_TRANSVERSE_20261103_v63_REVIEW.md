# 独立review：actual blowupによる二つの15平面の一般方向除去

2026年11月3日。対象packet
`Q_PLANE_BLOWUP_FINITE_TRANSVERSE_20261103_v63` を独立に検査し、下記の範囲で **ACCEPTED** とする。

## 受理した結論

二つの指定されたaffine15平面L0,L1は、Qstar=Q(2/9)のactual full35普通接錐の台に含まれない。
従ってそれぞれと台との交叉はproper closed affine coneで、次元は14以下である。

これは全平面の全点の除去でも、指定有理点p0自体の除去でもなく、full35局所環の全体次元上界でもない。
別途の全quartic-initial定理と組み合わせれば「五次以上の真の初期式が一般方向を除く」と言えるが、
今回のblowup除去自体は、そのquartic完全性を使わずに成立する。
ordinary arcの二次補正の不在だけを根拠にせず、Puiseux directionsも含むactual strict transformを扱っている。

## exactな独立再計算

原producerをfreshに再実行した後、producerおよびその`unpack` helperをimportせず、
凍結済みQ350、raw C350、二つの35-coordinate plane graphから独自のFraction sparse計算で次を再構成した。

* 14 base変数、20 normal変数、1 radial変数が全ambient blowup chartの35座標を与えること。
* p0上で全350 Qが0であること。
* 全350×21 linear行列 `[dQ|normal,C(v)]` の全entry。
* 保存19行minorと補19列から独立に解いた全21×2 normalized kernel。
* kernelのradial座標が全て厳密に0であること。
* quadratic polarizationから直接作った全350×3 raw係数と、19 implicit equationsを除去した全350×3 corrected係数。
* 各planeの全非零minorと、保存certificateとの全行一致。

二つの独立kernelと19次非零minorでlinear rankが正確に19、
corrected quadraticの3次非零minorで残る二変数の全quadratic spaceを張ることが証明される。
有限体標本・浮動小数・選択行だけの推測ではない。

保存minorは原packetと一致した。

| plane | linear determinant | corrected quadratic determinant |
|---|---|---|
| L0 | −956393211332282501953125/100352 | −27221241815728515625/311328730570208324171328 |
| L1 | −18169764766174127197265625/205641981272 | 9471284996723/21652080156000 |

独立reviewの測定4.64秒はfresh producer replayと全行再構成の両方を含む。
raw Q/Cがactual full350の二・三次部分であることは、それぞれの先行したexact packetに依存する。
今回それらの全full-model展開を再作成したという主張ではない。

## formal proofの監査

### 1. 全weak transformとclosed fiber

rの二乗を除いたactual full350式は

    G=Q(v+n)+r*C(v+n)+r²*D(v+n)+...

である。a-only denominatorはa=2/9+r*p_uでもunitなので、この式は実際の正則formal weak transformである。
linear rank19により、選択した19 unknown（rを含む）はR=QQ[[p−p0]]上でimplicitに消去できる。
残り二unknownのlinear termsは全て0である。

kernel上ではr=0だから、二次部分はQ(kernel)とlinear pivot eliminationの補正で正確に尽くされる。
r*Cのvariation、r²*D、高次tailを落とすことによる見落としはない。
corrected quadraticsが二変数maximal ideal mの平方全体を張るため、残るideal Hについて

    H ⊂ m²,  m²=H+m³

となる。有限module m²/HへのNakayamaによりH=m²である。
よってclosed transverse fiberは厳密にQQ[[s,t]]/(s,t)²、長さ3となる。

### 2. complete Nakayamaによる相対有限性

weak transform completion SはR[[s,t]]のquotientで、S/m_R Sが長さ3である。
従ってmaximal ideal m_Sの十分な冪はm_R Sに入り、m_S-adic topologyとm_R-adic topologyは同値である。
Sは後者でもcomplete/separatedとなる。
closed fiber basis 1,s,tを持ち上げ、complete Nakayamaの逐次近似を適用すると、SはR上高々三生成の有限moduleとなる。
R-flatness、generic rank3、rのS上のregularityは仮定していない。

### 3. saturationとcompletion

actual chartのweak idealをJとすると、strict transformはJ:r^inftyである。
Noetherianityによりascending principal colonsは有限Nで安定する。
flat completionは各r^Nのmultiplication kernelと可換するため、
安定したcolonをcompletionの前後で取って同じidealを得る。
従ってstrict transform completionは正確に

    T=S/(0:_S r^infty)

である。Tは有限R-moduleのquotientで有限のままである。
Tが0でなければ、その定義からrはnonzerodivisorであり、local maximal ideal内にある。
よってdim T≤dim R=14で、dim T/(r)≤13となる。

### 4. actual tangent coneとの比較

blowupのexceptional divisorはProj(gr_m(actual local ring))である。
このchartでそのcompletionはT/(r)である。
仮に指定15平面全体がactual接錐の台に含まれれば、そのprojective14-chartはp0を通ってexceptional divisorに含まれる。
するとこのcompleted local ringの次元は少なくとも14で、上の13以下に矛盾する。
従って平面全体は含まれず、閉性と既約性から交叉coneの次元は14以下になる。

## 凍結した独立review

* `code/q_plane_blowup_finite_transverse_review_20261103_v63.py`
* `reports/Q_PLANE_BLOWUP_FINITE_TRANSVERSE_20261103_v63_REVIEW/report.json`

raw入力・producer・証明ノート・certificateのSHAを検査し、このreviewのreportにも保存した。
full35の全枝、残るproper closed residual subsets、全initial idealの決定は今回の結論ではない。
