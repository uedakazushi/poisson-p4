# 固定T01一次部分はG28上にない

2026-11-06。Gの固定fiber診断で既に使ったT01の一次部分だけを対象とする。
元の点は`reports/S61_KNOWN_FIBER_AUDIT_20261031.json`のT01行である。

一回の有理数体上の正確な評価により

    G28(beta_T01) = 671334843179656461002307382635849 / 9039207968 ≠ 0

を得た。従ってこの固定βは、既知NEW25の閉像V(G28)には入らない。
他の48個のproper既知閉像を除外したという結果ではない。
一般βの尽くし、固定fiberのunit所属、標数零への有限体結果の移行も主張しない。

実装は`code/classification_frozen_beta_g28.py`、
結果は`reports/CLASSIFICATION_FROZEN_T01_BETA_G28_20261106.json`。
audited C16で元の重み1座標をprimitive b16へ移し、canonical B24×16との一致を確認した。
muはC,DとOmega/rhoの二通りで一致し、Hess(F/2)も117項の四次式とC,Dの微分からの
二通りで一致した。G28の78項はA14生成元を含まないため、その生成元は評価していない。
使用した全入力hashと11個の中間不変式値は結果JSONに保存した。

独立SOURCE査読とROOTGOの後、同じ既存点一つだけで実行した。
60秒・1GiBの上限で0.468秒で終了した。別の点・素数・因子・不変式を探索していない。

再現用の小SOURCE manifestは
`reports/CLASSIFICATION_FROZEN_T01_G28_SOURCE_20261106.json`。
元producerを変更せず、次のように指定のwork位置へコピーしてから新しい出力先で評価できる。

```sh
mkdir -p work
cp reports/CLASSIFICATION_FROZEN_T01_G28_SOURCE_20261106.json work/frozen-beta-g28-source.json
.venv/bin/python code/classification_frozen_beta_g28.py \
  --source-sha256 939da0876dcb538ea4417c0aa6daadcb0f4fbc7b6b42a59eb22436cf649cb42d \
  --output work/new-frozen-beta-g28-value.json
```

独立の最終78項検査は、Fractionによる元和とは別に、逆順の整数積と
全分母のLCMを用いて同じ有理数を得た。同じ法32003での値は25057≠0。
`reports/CLASSIFICATION_FROZEN_T01_G28_FINAL_SUM_REVIEW_20261106.json`を参照する。
