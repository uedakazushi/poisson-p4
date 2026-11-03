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

## 後続の完了単位：非平方 d=0

最初の push 後も研究を続け、**十分一般β上の指定非平方 d=0 全層を除去した**。
固定した (β,支持面) 上の実際の第二因子 P¹を明示し、有限部を元135式から
QQ[m] Smith 枠で再構築した。非定数対角因子を逆にせず、全561×269行加群が全列を張る。
再構築行列の SHA256 は `7fd3baf38b83fb9af4a43f377f9021959b9cb2b8995f1e1473dd0ab61ce49fa0`、
引継ぎと完全一致。全過程177.33秒。
実際の平方旗端を真coreへ有限群で移し、平坦な平方零商と実kernelの混合標数の議論、
飽和支持の正重み有限性を使う。全P¹の消滅から既約16次元の支持面因子上の proper image を
取って全層へ進めた。[最終証明](notes/S61_BINARY_D_ZERO_STRATUM_REMOVAL_20261103.md)、
[独立レビュー](notes/S61_BINARY_D_ZERO_BRIDGE_REVIEW_20261103.md)。

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/s61_binary_d_zero_prime_bridge_20261103.py --output work/d0-prime-new
.venv/bin/python code/s61_binary_d_zero_geometry_20261103.py --output work/d0-geometry-new.json
.venv/bin/python code/s61_binary_d_zero_pencil_20261103.py --output work/d0-pencil-new
```

## 現在の証明上の障害と次の操作

平方 `c=2a` の中央環 S の閉商 H=S/(ρ,七light) を二つのQQ数体上で厳密に計算した。
`1000h²+1127=0` には重み2の中央軸、`h²+2=0` には重み(2,3)の中央平面がある。
両Hで t=x4,t²=0だが、この消滅は周囲への証明ではない。
したがって全hで S/(ρ) 有限とする単純な判定経路は使えない。
[根基・行加群・数体証明書](notes/S61_MARKED_SIGMA_PENCIL_PROGRESS_20261103.md)を保存した。

次は同ノート§6に従い、新しい軸の x3=1 と平面の重み付き射影直線を、凍結した Smith
陰関数式へ戻す。元の (β,[q])、全135式、同じ滑らかな入射族と一般重複度を同定してから、
新しい横断上界を求める。固定中央 t²=0から周囲の消滅を結論しない。
全S5の別の有限体次元跳びも同ノート§5に記録した。
汎用 QQ(h) GB の既失敗計算を同条件で再実行する前に、この小さい閉層を処理する。

一般アフィン五点の尽くし、未知成分数の上界、軌道・交叉・非被約構造の全分類は未解決。
**成分数57は引き続き下界である。** a=0 と d=0 の除去だけで全 marked 非平方を
尽くしたとも結論していない。[shifted polar の厳密な係数模型と残る有限基準](notes/S61_NONSQUARE_MARKED_REMAINING_20261103.md)
を記録した。この模型は Jacobi 解でも実境界点でもない。新しい中央枝での周囲消滅に必要な同一入射族が、現在の
具体的な未完成の橋である。アクセスや権限の障害とは区別する。

v60 の日本語原稿統合と v61 の a=0,d=0 更新を行い、最終 PDF のページ数・検証記録は
`reports/MANUSCRIPT_SYNC_20261103_FINAL.md` に実測結果を記す。
ChatGPT Library への汎用 ZIP 保存機能は見つからず、Library 保存は未実行。
GitHub main への通常 commit/push を保存手段として使う。

## 保存済み checkpoint

最初の研究 checkpoint は `043955dc8ba47e0186ba8db3de87c2c43c983758`。
`research: remove the nonsquare a-zero stratum after Cloud resumption` を main に
通常 push し、GitHub 側で成功を確認した。この保存後も同じタスクで研究を継続した。

第二の研究 checkpoint は `46fb2456cbb1d0d2fee44657c7f5f4ad2d83b4ca`。d=0 の全層除去、平方中央支持、限定的係数模型、主原稿110ページ・評価14ページを通常 push した。

同コミットの clean な一時 checkout から継続ZIP `poisson_p4_continuation_20261103_v61-s61-cloud.zip` を作成した。2155ファイル、259871673 bytes、SHA256 `f69102a2af6dfd298b0dc4f7a8ab4ea1063220f1d6c9594fedfef40a11a2d38a`。
[保存記録](../reports/BUNDLE_V61_CLOUD_20261103.json)。ZIPは既定の ignored `bundles/` にあり、Library upload は未実行。再現コード・証明・原稿・全原資料と保存記録はGitHubにある。

正式な再現計算と原稿ビルドは完了し、未完の常駐研究計算は残していない。研究を継続する常駐ジョブ・timer・GitHub Actionsは作成していない。開始時の未追跡セットアップ報告はそのまま保持し、配布ZIPには混入させなかった。次の研究は上の数学的な橋から再開する。
