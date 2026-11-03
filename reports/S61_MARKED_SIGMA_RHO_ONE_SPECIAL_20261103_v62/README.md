# 有限 h の ρ=1 切断と追加点の探索入力

`rho_one_evidence.tar.gz` は今回実行した固定有限体 h の M2 入力・出力、
原点外の飽和点、線形 Smith 計算、p 進探索のコードと結果を保存する。
各ファイルの SHA256 と実行の射程は `manifest.json` に記録した。
数値出力と標数零の証明は区別する。

```sh
mkdir -p work/rho-one-replay
tar -xzf reports/S61_MARKED_SIGMA_RHO_ONE_SPECIAL_20261103_v62/rho_one_evidence.tar.gz -C work/rho-one-replay
source /workspace/cloud-setup/activate.sh
M2 --script work/rho-one-replay/h-4499-extra.m2
```

同じ固定 h の GB は約60秒で完了し、各回90〜100秒の実行上限を設けた。
四つの追加点の実行終了値は個別に0を回収した。初回の6ケースはバッチ全体の終了値0と
各完全出力を回収したが、個別終了値を別ファイルへは記録していない。

これらの M2 入力は、凍結した全 h Smith 模型の `global-rows` の有理係数を
F32003 に還元し、ρ=1 とした全135式である。第136行は元の t の表示で、
追加方程式には使用しない。元の全係数の SHA は manifest に保存した。

探索で得た二次数体と法線座標を使い、別の正式 producer
`code/s61_marked_sigma_finite_extras_20261103_v62.py` が元135式への厳密な再代入と
全周囲横断を確認する。p 進再構成の低精度の誤候補も探索ログへ残し、
最終的に採用した正確な座標と混同しない。
