# 四 core 二次初期式のための source 三次補正消去

2026-11-04。新手法の有限入力。ここでは三次係数そのものはまだ評価しない。

前提は実 N-pair hypersurface と `NORMALIZER_CORE4_REDUCTION` の正確な入力である。
その正則九座標環 S と十二関数 F_j=L G_j を用いる。
既存二次キャッシュから G の一次項は五次元で、選択 rows は18,19,60,28,29。
残る七 rows は67,72,63,64,68,73,74。このうち row72 のみ row60 を引く。
これにより七つの修正関数 F'_j の全二次 jet は S 上で恒等的に零となる。
従って G'_j の一次 jet も零であり

    (F'_j)_3 = L_1 * (G'_j)_2

が成り立つ。L_2 や L_3 は不要である。
五つの G 選択関数を解く形式 IFT は、その消去変数の二次以上の加速度を変えるが、
一次項零の G'_j の二次項を変えない。従って五物理 tangent 座標
x108,x109,x113,x118,x119 を単に零としてよい。
残る変数は x112,x117,tauA,tauB の四つであり、三次以下の単項式は35個。

source は元の scalar DAG と実98 unit 方程式 E(u,c)=0 を使う。u は96係数未知数と
rho27,rho61、c は四 core 座標である。保存済み二族57 source 二次 jet を
既存の物理 R0 逆写像とモニック Schur 二次写像へ代入し、u の二次近似 u_[2] を得る。
未知の三次補正を u_3 と書き、J=E_u(0)、B_j=(F'_j)_u(0) とすると

    J*u_3 = -[E(u_[2],c)]_3,
    [F'_j(u(c),c)]_3
      = [F'_j(u_[2],c)]_3 - lambda_j*[E(u_[2],c)]_3,
    lambda_j*J=B_j.

J は実98 unit 行列。今回の preflight は七つの lambda 行を QQ 上で解き、
全686係数について lambda*J=B を直接照合した。
従って未知 source 三次係数を零と仮定せず、それらへの依存を厳密に相殺する。
source 三次解を新しく求める必要はない。元 source の四次係数も使わない。

有限入力には全98方程式・未知数、七つの literal residual 結合、98ヤコビ行列、
七 adjoint 行、二族の元57二次キャッシュから四変数への全代入、必要な元 DAG 祖先を保存する。
祖先17085、乗法9750、35単項式の乗法畳込みは165 ordered pairs なので、
最悪の dense 乗法処理も1608750 scalar pairs。60秒・2GB・32MiB の一回評価に十分小さい。

次の一回計算の完了条件は、全98方程式の二次以下零を検算し、上の adjoint で補正した
七つの三次式を既知

    L_1=5*x112-20*x117+24*tauA+6*tauB

で正確に割り、七つの QQ 二次式・全35係数の除法恒等式・その span rank を保存すること。
除法余りが非零ならそこで停止する。二次式の結果から全 ideal・socle・深さを推論しない。
自動の core 三次・source 四次への引き上げは行わない。

再現入口:

    .venv/bin/python code/ambient_B18_core4_quadratic_preflight_20261104_v63.py --output work/B18-core4-quadratic-preflight

この入口は preflight のみであり、三次係数評価を実行しない。
