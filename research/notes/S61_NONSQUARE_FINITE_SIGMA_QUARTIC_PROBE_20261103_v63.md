# 有限Sigma上のraw四次：一回のGF rank probe

2026-11-03。rootが指定した bounded actual test。既存p=32003、有限algebra Sbarの次元28と、独立検証済みの二つのregular multiplication operatorsを再利用する。新しいGB・field・primeは生成しない。

入力は独立finite-Sigma reviewの `own-operators.txt` と、元804項から作り全24行を照合した `own-raw-modp.txt` である。後者は24×56、q_i*x_jの係数であり、各q_i*x0行から四quadricsを剥がす。その四式を用いて元24行の全係数を再構成し、一致を確認してから21個のquadratic monomial倍を作る。

84×126の各係数を28×28regular operatorに置き換え、実2352×3528行列を作る。aの作用素は既知DB=0の式 `a=−P(20P+3e+20)/180` で定義する。ブロックはSbar^126→Sbar^84の列作用であり、乗算operatorを恣意的に転置しない。

一回だけFLINT `nmod_mat_rank` を呼ぶ。期待rankは78×28=2184。60秒・2GBのcapを設け、timeoutなら未完のfrontierだけを保存する。全実blockはlittle-endian uint32 header/entriesでhashし、公開receiptにはSHAとcompact source bindingを保存する。二回目rank、巨大inverse、別fieldを自動起動しない。

rankが2184なら、六Koszulがsplitであることと既存の21×4rightinverseのprovenanceを別に監査した後、有限proper model/Nakayamaへの移行を提案する。このprobeだけではQQに昇格しない。rankが低くてもQQでの欠損を結論しない。

実行結果はreportのnative_completedとrankを見る。実行前に成功を仮定しない。

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_finite_sigma_quartic_probe_20261103_v63.py \
  --output work/marked-square-v63/finite-Sigma-quartic-probe-new
```

main C122のactual rank64問題はoracle構築後の別タスクとして保持する。
