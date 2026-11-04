# B18 四 core: 七つの実二次初期式

2026-11-04。元の全次数の core 関数は `NORMALIZER_CORE4_REDUCTION` の有限入力に固定する。
ここで計算するのは、その七関数の二次初期式だけであり、全 ideal や socle ではない。

`CORE4_QUADRATIC_PREFLIGHT` は、七つの修正出力 F'（row72 では row60 を引く）が
Schur 九変数で二次以下零であること、実98 unit 方程式 E のヤコビ行列 J と
七行 lambda に対する lambda*J=dF'/du の全686恒等式を固定した。
従って保存済み source 二次近似 u_[2] だけを使って

    cubic(F' after exact unit solve)
      = cubic(F'(u_[2])) - lambda*cubic(E(u_[2]))

が厳密に成り立つ。未知の三次 source 補正を零と仮定していない。
この adjoint によって全98個の未知三次係数を相殺した。

四変数の順序は x112,x117,tauA,tauB。全35単項式を持つ QQ jet interpreter で
元 scalar DAG の必要17085祖先を評価する。全98方程式が二次以下零であることを再検算し、
七補正式を

    L1=5*x112-20*x117+24*tauA+6*tauB

で正確に割る。各35係数で積が元の式へ戻ることを確認する。
こうして得る七つの二次式は、五つの G の nonlinear IFT 後の七 core 関数の
実二次初期式である。G' は一次以下零なので、その IFT の二次加速度は影響しない。
N-pair の二次・三次の新しい係数も使わない。

全七式、7x10 QQ 係数行列、span rank、非零 minor は添付 report と
`seven-actual-core-quadrics.json` に保存する。ここで単に零だった係数を高次へ
自動的に引き上げない。二次初期式から全 ideal・nilradical・socle・E4 または D0 の深さは
推論しない。source 三次求解、新 source 四次、GB、追加点は一切行わない。

再現:

    .venv/bin/python code/ambient_B18_core4_quadratic_initial_20261104_v63.py --output work/B18-core4-quadratic-replay

これは許可された一回の四変数三次評価であり、既存 source 二次キャッシュ以外の
新しい source 解は生成しない。
