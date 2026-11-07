# 残る12場合の OLD4 二次応答を同時に消す係数核 — SOURCE

これは既存の有理定数係数の rank/kernel 計算を準備する SOURCE であり、rank、上側応答、未知の二次係数はまだ評価しない。

残る場合を `(579,580,593,594,597,598,605,606,611,613,615,616)` とする。各場合はすでに563/609の元の source cycle を引いた正確な617側の source である。4本の OLD pivot の一次補正係数の r,ell 成分を並べた12×8行列を W とする。元の full617 OLD2 記号証明と一致する affine SOURCE ac8497 の corrected plane-a1 だけを用いる。

各未知 ordinary rr 応答の係数は `-a1_ell*Gamma`、各 rell 応答の係数は `-a1_r*Gamma` である。したがって `lambda W=0` は全24記号の両チャンネルへの寄与を同時に消す。これは各関数の二次応答を0とする主張ではない。得られる source は元の617 corrected source の線形結合であり、lower4 は0である。上側2成分の完全な既知部分は

`[Gamma_t*(Xc3-Xp2*a1)]_(r^2 ell) + D_t*(Xc2-Xp1*a1)`

である。full Psi、全27/99の消去、native individual L2*S0、既知 FIRST、full563/609 source correction を含む。既知部分自体はこの gate では計算しない。native source は full S0 kernel の制限後にのみ actual であり、raw auxiliary の jet を actual と呼ばない。

行ごとに `D_j=lcm(den W_j)`、`g_j=content(D_j W_j)`、`P_j=D_j W_j/g_j=s_j W_j` とする。さらに `L=lcm(den s_j)`、`t_j=L*s_j` と置く。すべて非零の既知有理数による source 行変換である。P の行 L1 norm を1以上で切り上げ、大きい8個の積 B を取れば、Hadamard により P のすべての r×r minor (r<=8) の絶対値は B 以下である。これは rank を求めずに計算できる。

準備する一回の worker は P の fraction-free full-pivot Bareiss、選択した original row/column minor、全12×8の Cramer factor identities、全 free row の source kernel を返す。Bareiss の committed entries は original minors、乗積差は `2B^2` 以下で、除算は exact と検査する。基底 minor の determinant d と、行 i を行 j に取り替えた determinant n_i に対して original source kernel を `d*t_j*e_j-sum_i(n_i*t_i*e_i)` と保存する。全8座標を W の共通分母で整数化した元の行列でも検査する。free block の対角 `d*t_j` は非零なので核ベクトルは独立であり、全12行の factorization が完全性を与える。

rank<=8 なので少なくとも4本の complete opaque-free source combination が得られる。これは値や rank の推測ではない。579/613/616を除く9行を使う代案なら少なくとも1本となり、true conductor third も避けられる。ここでは完全な W12 を優先する。必要な conductor third は元の own-chart の真の leading cubic と既知の linear chart map だけで決まり、RAW45 の literal-existing-cubic-source-input に既にある。global m^3 と pointed composition により新しい physical Phi2/Phi3 は不要であるが、third 自体を0とはしない。gA は prefix0[7790]、gB は own-CB conductor であり unrestricted f63 とは区別する。これらの pairing も今回の worker には含めない。

整数の前処理、prior、メタデータ比較だけを今回の SOURCE 実行に含める。rank worker は immutable code/custody、独立 SOURCE/code peer、ROOT の明示 GO がそろうまで実行しない。固定 cap は committed384bit、temporary1153bit、60秒、2GB、numeric12000、RAW32MiB。compute55秒と snapshot5秒を分け、STOP時は exact matrix、pivot/cursor、未確定演算、factor/kernel prefix を保存する。timeout は外部60秒でも制限する。新しい source/physical/gauge/IFT、上側12値、未知二次値は範囲外である。

V2 bound clarification: worker が各 guarded binary result を保存するため、committed prior には minor/kernel だけでなく product/partial-sum prior も含める。V1 の79bitは前者のみ、V2では全結果129bit（temporary129bit）で、hard384/1153は変更しない。V1 bytes は保持する。
