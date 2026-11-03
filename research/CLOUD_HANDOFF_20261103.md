# Cloud 研究継続の記録 — 2026年11月3日

公開済み Cloud 環境で研究を再開した。開始時の `main` と fetch 後の GitHub main は
`4c3f623fb698f9adc6077015b25e2325b11072af`。この記録は移行作業ではなく研究作業の記録である。
11月2日のローカル移行用停止を、Cloud の停止指示には適用していない。

## 環境と再現入口

既存の `/workspace/cloud-setup/activate.sh` を各 shell で source する。
Python 3.12.14、固定 Python 依存、FLINT 3.1.0/GMP 6.3.0/MPFR 4.2.2、
Macaulay2 1.24.11、日本語 LuaLaTeX を再利用した。セットアップ全体は再実行していない。
[軽量確認と互換性](../reports/CLOUD_RESUME_20261103.md)、
[文献・GitHub 接続確認](../reports/CLOUD_CONNECTIVITY_20261103.md)。

ネットワークを使う `exec_command` には、正規の
`sandbox_permissions=with_additional_permissions` と
`additional_permissions.network.enabled=true` が必要だった。
inherited proxy と TLS 検証を保つ。7文献サイトと複数の一次資料に接続成功した。
Wiley への DOI 転送はホスト未許可で403、export arXiv API は一度20秒 timeout。
これらは現在の局所代数・幾何の研究を妨げない。

## 完了した研究単位

**十分一般の一次ジェット上で、非平方 `K=a=0` の全層を除去した。**

- 中央8例を QQ と F32003 で元135式から再構成し、明示した七 ray の和との根基等式、
  原点支持の二例の長70・65を厳密に確認。
- 七代表の同じ入射族の軌道型、全最高式、p単元切断、一般βへの優越性を別実装で確認。
  指定β自身の微分階数は全て15であり、そこでβ写像が滑らかとは仮定していない。
- 横断長上界と既知の同一入射族の一般長の下界から、混合標数上の全周囲の冪零性を導いた。
  非有限中央環は実際の t 飽和商の有限性と重みで処理した。
- 正しい `c=0` の無限遠端点を含む完全試験曲線を用い、固有性と P² の半不変因子の
  排除から全層へ接続。不良像の次元は15以下。

証明：[全層への持上げ](notes/S61_BINARY_A_ZERO_LIFT_20261103.md)。
独立点検：[計算・幾何・証明の監査](notes/S61_BINARY_A_ZERO_COMPUTATION_REVIEW_20261103.md)。
台帳の新規 ID は `S61_BINARY_A_ZERO_CENTRAL_SUPPORT`、
`S61_BINARY_A_ZERO_RAY_GEOMETRY`、`S61_BINARY_A_ZERO_STRATUM_REMOVAL`。

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/s61_binary_a_zero_replay_20261103.py --output work/a0-fresh-report
.venv/bin/python code/s61_binary_a_zero_ray_geometry_20261103.py --output work/a0-geometry-fresh.json
M2 --script reports/S61_BINARY_A_ZERO_INFINITY_MODP_20261103.m2
```

出力先は毎回新しい名前にする。正式入力・全出力は
`reports/S61_BINARY_A_ZERO_CLOUD_20261103/`、独立 geometry の JSON、
`reports/S61_BINARY_A_ZERO_INFINITY_MODP_20261103.*` にある。
531×386 の行加群は522662演算、法32003で行列式次数17。
有限体の因子を QQ の因子と同一視する代わりに、混合標数の固有性を使った。

## 継続中の課題

次の研究単位として `d=0` の旗端点の混合標数の橋と、平方 `c=2a` の有限値を扱っている。
それぞれの新しい計算・証明の採否は STATE/NEXT_TASKS と最終保存時に更新する。
一般アフィン五点の尽くし、未知成分数の上界、軌道・交叉・非被約構造の全分類は未解決。
**成分数57は引き続き下界である。**

v60 の日本語原稿統合と v61 への更新も進めている。最終 PDF ページ数・検証記録は
`reports/MANUSCRIPT_SYNC_20261103.md` に実測結果を記す。
ChatGPT Library への汎用 ZIP 保存機能は見つからず、Library 保存は未実行。
GitHub main への通常 commit/push を保存手段として使う。
