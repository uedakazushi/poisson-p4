# 四次全射に対する286×238の十分条件

2026-11-03。親の提案した有限 corollary を、凍結済みの319行縮約の独立監査と heavy-square 全開 pivot から導く。実際の rank238 は未証明である。

全 aFBC-open の同じ weighted normal chart を用いる。低重みの prolongation は、6二次式×29単項式と21 genuine cubic×7 light の321行である。二次式は pure light、三次式の単項式は light³ または heavy×light なので、これら321行の heavy² 列は全て零である。

独立監査済みの35恒等式は、全てこの321行だけを使っている。従って同じ unit minor `−a^14` により低重み321行を **286行** に減らして row span を保つ。heavy² を除く238列に射影しても等号は保たれる。

既存の独立監査済み actual 行77・131の線形結合は heavy² 係数 `(50aBC)^2` を持ち、これは全 aFBC-open の単元である。今回の小 verifier でも保存された実二係数と行係数を別に代入し、この恒等式を確認する。従って286低重み行が238列を張る点では、この heavy-square 行を追加すると weight4 の239列全体が張られる。

これは **286×238 の十分条件** であり、全四次 ideal と同値な319×239判定をさらに286×238へ等価変形したものではない。省いた33 genuine quartic 行が不足を補う可能性を残す。d、J、product pointwise41、heavy7、full cubic57 の仮定は、この低重み row-span 等号に必要ない。新しい source 展開・小行列式探索・GB・体・未知因子の局所化は行わない。

再現:

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_quartic_286_sufficient_20261103_v63.py \
  --output work/marked-square-v63/quartic286-new
```

公開 report はこの有限 recipe と恒等式を記録し、rank238・長65・ambient removal の成立は主張しない。
