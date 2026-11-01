# 基点の移動に対する primitive 一次部分の普遍的輸送

2026-11-01。ここで証明するのは座標変化の公式である。曲線不変量の保存や、追加式の
全 S61 イデアルへの所属は仮定しない。

## 1. 公式

アフィン双ベクトルを `P(x)`、`P(0)=J`、`L_i=∂_iP(0)` とする。
基点を `ε e_h` へ動かし、

\[
 A_h=\tfrac12 L_hJ,\qquad F_h=I+\varepsilon A_h,
 \qquad z=F_h(x-\varepsilon e_h)
\]

と置く。`A_hJ+JA_hᵀ+L_h=0` なので、変換後の定数部分は ε の一次まで J である。
線形部分の変分は連鎖律から

\[
 \delta_hL_i=\partial_h\partial_iP(0)+A_hL_i+L_iA_h^t
 -\sum_k(A_h)_{ki}L_k. \tag{1}
\]

最後に、既存の線形 primitive projection と `β=−JLJ` を適用する。
Darboux frame の別の選択との差は infinitesimal Sp4 作用であり、Sp4 不変量の微分には影響しない。

S61 再構成では二次部分は q にアフィン線形で、βだけの部分は二次である。
従って (1) は普遍的に

\[
 \delta_h\beta=M_hq+B_h(\beta), \tag{2}
\]

となる。`M_h` は β に依存しない16×45有理行列、`B_h` は β の二次ベクトル。
四つを積んだ64×45行列の階数は45である。

全16一次変数と45二次変数を記号のまま残し、各変分が閉・primitive であること、
S61 の16基底で一意に表示できること、(2)の次数、定数行列の階数を厳密に確認した。
これは残り134式を仮定せず、低次 Jacobi を解いた既存の再構成上で成立する。

## 2. 不変量による消去候補との境界

内在的な曲線 `Cβ` の Fisher 不変量は β 次数8の c4 と次数12の c6。
(2)をそれらの微分へ代入すると、q にアフィン線形な四行ずつを得る。
重み付き次数はそれぞれ9と13である。
`2c4 dc6−3c6 dc4` を用いる j の微分の分子は β 次数19、輸送後の重み付き次数21。

保存済みトーリック一次部分で、c4 の四行と c6 の四行を併せた階数は8、
j の四行の階数は4となった。独立実装の証明書は
`reports/PRIMITIVE_CURVE_INVARIANT_GRADIENTS_20261101.json` を参照。
ただし、これらを未知成分まで含む全 S61 の消去に使うには、全 Poisson イデアルでの
保存則または明示的な所属証明書が必要である。既知成分上で消えることだけを根拠に
追加の定義式として課してはならない。

その後、射影 Poisson 構造で Cβ が滑らか、すなわち Δ≠0 であっても、基点の移動で
j が変化する厳密な反例が得られた。従って全 S61、さらに Δ による飽和における
普遍的保存則という予想は反証された。
[反例と全 S61 の確認](CURVE_INVARIANT_SMOOTH_COUNTEREXAMPLE_20261101.md)を参照。
この反例は輸送公式(1)(2)を反証するものではなく、その値を常に零と置く使い方を排除する。

## 3. 再現

固定一次部分での輸送:

```sh
.venv/bin/python code/s61_jet_basepoint_transport_20261101.py \
  --output reports/new-jet-transport-sample.json
```

全16一次変数を保持する輸送:

```sh
.venv/bin/python code/s61_jet_transport_universal_20261101.py \
  --output reports/new-jet-transport-universal.json
```

保存済み結果は `S61_JET_BASEPOINT_TRANSPORT_20261101.json` と
`S61_JET_TRANSPORT_UNIVERSAL_20261101.json`。q は中心化していない元の重み2座標である。
各行列・二次ドリフト・変数対応・入力ハッシュを保存する。
