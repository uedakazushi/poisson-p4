# 現行研究課題 — Poisson接錐の三次初期関係と成分交叉

2026-11-10 / v16。再開はRESEARCH_POLICY.mdから。GitHub main正本、生成PDF非追跡、GitHub Actions不使用、Libraryは実push失敗時の緊急退避のみ。

**証明付きの到達点**：完全等方Log点pの接空間30次元、局所次元28、真の接錐に非零平方零元h=a(s2+c)、二次初期イデアルは7生成元、Hilbert–Samuel重複度1。最大次元の28次元Log1112成分はpで滑らか。25次元long rank4成分との共通23次元Poisson族を明示し、交叉の局所次元は23か24。固定Aを保つ5パラメータの滑らかな二次式には普遍Poisson延長があり、Aと基底三形式直線を保つ二次形式はrank(q|W)=0/1の二軌道。g_ρ=diag(1,1,ρ,1,ρ^-1)は(η,λ)→(ρ²η,ρλ)と作用し、rank2完全等方型へ退化する。

**新しい同変分解**：30接方向の重み(-2,-1,0,1,2)別次元は(1,5,12,10,2)。95正規方向の重み(-2,-1,0,1,2,3,4)別次元は(2,7,17,26,28,12,3)。7二次初期関係の重みは(2,1,2,3,1,1,1)、平方零hは重み2。30固有ベクトルと重み表はコードcheck_nilpotent_tangent_weights.pyで独立検算。

**次に着手**：95正規方向を重みブロックごとに二次まで消去し、三次残差と7本の二次関係のsyzygyを組み合わせて全初期イデアルの三次部分を求める。同時に固定A・三形式の対称性を破る共通Poisson変形の有無を調べ、交叉次元23/24を決定する。26/24次元の二次模型成分は実接錐の成分と未同定。ambient局所環の被約性も未解決で、接錐の非被約性から推論しない。

原稿：papers/isotropic_log_tangent_cone_ja.tex と共有 papers/isotropic_nilpotent_pencil.tex。計算：code/subvariety_incidence/check_isotropic_quadratic_cone.py、check_nilpotent_torus_attraction.py、check_nilpotent_tangent_weights.py 等と各certificates。手動入口：bash scripts/verify_github_source.sh tangent-cone。厳密計算、LuaLaTeX組版、GitHub保存は別々に検証する。
