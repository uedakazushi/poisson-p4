# actual direct P4＋H5：型付き chart の独立監査

2026-11-04。全151節点の型と外部 frozen sourceを独立に読み、全70 heavy-root
対応、P4/H5の行列順、hybrid単位、元25 sourceへの covector輸送を確認した。
この receipt は指定一成分の新しい算術評価を含まない。別の bounded arithmetic
検査は48秒ガードで停止したため、その7022係数一致は未受理である。

P4の実 current15行は `[0,1,3,6]`、original32標的は `[22,11,13,25]`。
各P4行のheavy5は構造的に0。H5の実行は `[2,4,7,13,10]`、heavy標的順は
`[27,31,28,30,29]`。H5の5行だけrhoで割り、heavy31だけgammaで割るので、
P4は未変更であり、既受理の actual `D4=a^4 C^3 R434` をそのまま使える。

H5の2＋3 block順は既受理の pair2 / 選択B3 minor の順と一致する。
`Hnum=[[delta3 Adj2,0],[-Adj3 C32 Adj2,delta2 Adj3]]` は、一般の可換係数環上の
cofactor恒等式により、H5との両側積が `Delta5 I5` となる。
したがって H5 inverseは `Hnum/Delta5`、P4 inverseは `Adj4(P4)/D4`。
h22を逆にせず、領域は正確に `FULLU∩D(R434 Delta5)` である。

実9 pivot blockを `K=[[P4,0],[C54,H5]]` と書くと、そのinverseは
`[[P4inv,0],[-H5inv C54 P4inv,H5inv]]`。
全行列 `M=[[K,B],[C,D]]` に対して表示された可逆な左右作用素は
`diag(I9,D-C Kinv B)` を与える。従って元 actual15×22 と residual6×13 は同一cokernel。
全8個のQ9追加pure tailが元の実 root のまま残13標的に保持されている。

残6 sourceの元25への輸送は、既存Gaussian stepのforward行作用を時間順に合成し、
actual残17選択、S17と削除、S16と削除を続けたもの。covectorはこのforward作用素を
右へ掛けるのであって、transposeや逆行列に置き換えない。

R434=0やDelta5=0の側は除去したとはしない。残行列の階数6、global cubic57、
quartic239、未知kernelの自由基底、他成分の係数・小行列式は主張しない。

再現（source/type読取のみ、scalar workerを起動しない）：
`.venv/bin/python code/s61_nonsquare_prePOWER_direct_schur_entry_review_20261104_v63.py --chart-only`
