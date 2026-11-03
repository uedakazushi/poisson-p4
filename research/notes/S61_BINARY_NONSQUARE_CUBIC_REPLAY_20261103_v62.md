# 全135式の三次行列を Cloud 上で再構成

2026-11-03 / v62。Macaulay2 のバージョン差とは独立に、既存の QQ 多項式 producer と公開済み Cloud の FLINT 3.1.0 を使い、元135式から weighted cubic 行列を再構成した。変更した実行条件は native compile 呼出しだけで、線形 section、陰関数の三次補正、既知単元による21 pivot は元 producer の全 assertion を通した。

再構成された63×91行列と42×70 Schur 行列の SHA256 は、それぞれ

    dae59904925260605628d925ae0565de084c9af0ae1bccf1774f1af067830ad6
    0946ea9d6b5ce952d12c6ac7095c91b4f7ef60f20c05fd2b4f558f5c00c0a0a2

であり、既存 `S61_BINARY_HYPERPLANE_CUBIC_MODULAR_20261102` の元行列 hash と一致する。

引継ぎ archive `HYPERPLANE_CUBIC_D_DEGREE_20261102.tar.gz` の member `native_data/hyperplane-full-cubic-schur36-input.txt` も比較した。こちらの header は求める階数36を記し、その後に2940個の多項式 entry と一つの minor の選択が続く。現行 Schur file の header は除去した既知単元 pivot 数21を記し、次の二行は pivot 列・残行の metadata である。この format 差を明示して取り除くと、**全2940 entry の文字列が同じ順に完全一致する**。archive 全体や header の hash 不一致を、行列自体の不一致と扱う必要はなくなった。

実行は25.33秒。新 [producer](../../code/s61_binary_nonsquare_cubic_replay_20261103_v62.py) と [report](../../reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/report.json) に source hash と全検査を保存した。同ディレクトリの `matrices.tar.gz` は新しい `full`, `schur`, `weights`, `row_meta` を保持し、native binary や45MBの中間展開は含めない。

    source /workspace/cloud-setup/activate.sh
    .venv/bin/python code/s61_binary_nonsquare_cubic_replay_20261103_v62.py \
      --output /tmp/nonsquare-cubic-replay-new

これは正確な入力同一性の検証である。全 core における階数36、四次全射、横断長65、ambient t 冪零性の新しい結論は単独では含まない。既存 archive を使った探索計算を証明へ昇格する場合も、その計算自体の assertion と適用範囲を別に保存する。
