# Cloud 引継ぎ — 有限候補と固定二座標切断

2026-11-02。移行指示で研究・検算を停止した時点の記録。ここで新しい数学的再実行は行っていない。親が commit を担当する。

## 保存済みの結果

* `S61_FIXED_PAIR_SINGLE_COLON_20261102` が有限長評価の出発点。十分一般の β と任意の c=(u6,u42) に対し、全135式をこの二座標で切った斉次 ideal I_c の Proj(I_c:t²) は零次元以下で、アフィン部分を保存する。一般アフィンファイバーが五点に尽きること、この射影候補の長さ、またはその t 可逆部分の長さは未決。
* `S61_PROJECTIVE_TOP_CYCLE_20261102` の次数2300は未飽和三次元成分の値。一般三切断で t=0 商は2270。これらは有限アフィン候補の次数上界ではない。
* `S61_AFFINE_HOMOGENIZER_FINITE_ACTION_20261102` は二座標を切った元134式の20法線上の有限191生成模型。各射影法線方向における実代数を交換子の最大同時不変核から求め、t の191乗の階数がその方向の t≠0 部分の長さ。これは135式への追加が t を反転した chart で同じ理想であることを使う。一方向の値や幾つかの標本の和を全方向の長さにしてはいけない。
* `S61_AFFINE_FIXED_PAIR_UNIVARIATE_20261102` は166次の一変数表示。巡回行列式を反転した底上の同型。全方向には上の191生成法を使える。
* `S61_FIXED_PAIR_INFINITY_CURVES_20261102` は **親の最終取り込み・独立確認前の原稿**。code/report/note の三ファイルは既に正式ディレクトリにあり、入力/ハッシュは manifest 記載。C 上の15次 section と15本の垂直線、旗上の6本の有理曲線（計22）の候補を記述。三つの Plücker 切断の単純性・対の交点回避、両写像による旗グラフの相異性、重複値の次数和の確認が重点。原稿にはこの論証を記載済みだが、Cloud 移行後に親が確認してから採否を決める。既にある JSON だけを理由に凍結済みとみなさない。

## 途中計算

`finite_candidate_scratch.tar.gz` に29ファイルの入力・出力・ソースを収容。全件を manifest の members に記録。元のレイアウトはアーカイブ内の `work/` と `tmp/`。ソースが古いローカル絶対パスを参照する場合は抽出後にそのパスを直す。コンパイル済み macOS バイナリは除いた。

1. `work/s61_affine_full20_colon_linear_20261102.py/.cpp` は指定 β*,c=0 の全135・全44斉次座標による F32003 の三次検査。補助的な出力 `tmp/s61-affine-full20-colon-linear-20261102/output` の最終4行は 713、92、0、0。該当次数の92個の組合せは全係数が恒等零となり、この有限体標本の未飽和コロンに一次生成元は出なかった。標数零、飽和後の長さ、高次に関する結果ではない。134式版には691と91が記録されている。
2. `work/s61_affine_full_commutator_20261102.*` と `work/s61_affine_homogeneous_t_matrix_20261102.*` は191生成正式版の前段。最終判定は正式版を使用する。有限体の中間入力は既に各ディレクトリに保存され、再生成方法も含む。
3. まだ着手していない案は、Supp(t²³)∩∞ の次元上界を H11 の実 t² での除去範囲、旗と二平面の交叉の一般 t²、T5 の t²³ などから点検し、**一般の別の二線形座標**で無限遠を完全に避けられるかを調べること。未査読の着想であり、既証の固定 (u6,u42) に対する t² コロン零次元定理と区別する。指数23で一般ファイバーの長さが既に求まったという主張はない。

## Cloud で再開する場合

リポジトリのルートから実行する。以下は再開コマンドの記録であり、移行時には実行していない。既存 report を上書きしない。

```sh
python code/s61_fixed_pair_infinity_curves_20261102.py --output /tmp/cloud-fixed-pair-curves.json
python code/s61_fixed_pair_single_colon_20261102.py --output /tmp/cloud-fixed-pair-single-colon.json
python code/s61_affine_homogenizer_finite_action_20261102.py --output /tmp/cloud-affine-finite-t.json
tar -xzf research/handoff/20261102/elliptic_jets/finite_candidate_scratch.tar.gz work
python work/s61_affine_full20_colon_linear_20261102.py
```

旧 Mac の native ビルド関数 `code/s61_binary_hyperplane_cubic_matrix_20261102.py::compile_native` は Homebrew と `.dylib` の位置を使う。Linux Cloud では対応する FLINT、GMP、MPFR、C++17 と動的ライブラリの位置を先に用意する。各 report はローカルの正式 source hash を保持する。移植時の差分と再実行結果は別に記録する。
