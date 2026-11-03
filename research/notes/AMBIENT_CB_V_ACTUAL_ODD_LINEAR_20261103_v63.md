# CB generic V：modulusを残した全actual odd firstjet rank9

2026-11-03。作者の厳密QQ(V)計算が5.43秒で完了した。独立監査は未了であり、以下の全完成環への適用は監査対象として分ける。
CAのV1 single-orbit argumentをCBへ流用しない。

## 全sourceの有限計算

`AMBIENT_CB_V_FAMILY_FRAME_20261103_v63` の19 GL変数と10 fiber変数を使う。
finite sourceはordered GL exponentialsをfull CB ratio/f/H tensorへ作用させるものである。
特にfiberの二次項 −delta_q*sum(delta_f_i*G_i) をHessianとpure-kernel cubicの両方に含める。
全150係数のkernel-mixed HessianをQQ[V]で作り、別のunivariate ordered Lie expansionのkernel二次係数と一致を確認する。
kernel三次係数も全150個保存する。

28 actual target座標Yと一actual branch座標からquadratic base Rを取る。
kernel-complementをJ_r、そのY-minor逆行列をDとし、rem(M)=M−J_r*D*M_Yと書く。
W=F_kk/2、H=F_kr、Z=F_kkk/6に対し

    Wc=rem(W),   Zc=rem(Z−H*D*W_Y),   E=rem(H*D)

とする。branchのWc係数をw、linear branch correctionをLとする。
normal(M)=M−Wc*M_branch/w で補正すると、actual odd firstjetは

    [normal(E)−normal(Zc)/w*L, normal(Zc)/w]

である。全操作はQQ(V)上で行い、Vを数値代入しない。

この全150×29行列のrankは **9**、whole source critical cotangentのrankも **9** である。
保存した九次minorは

    157837977*V^28*(V²−2)^9 /125440000000000000000000000000000000000

で非零である。全column factorizationも保存し、rank上界9を厳密に確認する。
これは九つのsample値からの零推測ではなく、QQ(V)上の全行列identityである。
frameのV=0、V²=2除外はこの証明用frameの除外であり、その点でactual型が変わると主張しない。

## 全完成環への適用鎖（独立監査用）

P0の全proper CB flag fiberは一台長2である。
proper写像がfiniteでfiber長2以下となるtarget開集合を取り、その逆像をV曲線に制限すればV=0を含む非空開集合となる。
従ってgeneric Vでは全proper fiberが有限で長2以下である。
source Jacobianのcorank1により、その全fiberは一台長2となる。
この議論はV=2などの任意の一点自体がその開集合に属すると仮定しない。

既知の全CB source critical idealはP0近傍で滑らかな19因子と
(WU,WV,UV,U³)であり、generic V枝では正則20次元となる。
今回のwhole critical cotangent rank9もその正則性を与える。
quadratic base Rは正則29次元、finite sourceはB=R[t]/(t²−nu)、actual imageはA=R⊕I*tである。
J=(I,nu)の商はこのcritical quotientなのでJは正則九生成idealである。
actual Iのfirstjet rank9はJ/mJを全て張り、NakayamaによりI=Jとなる。
I≠0からt=(ct)/c∈Frac(A)を得て、ここで初めてfinite normal Bがactual normalizationであると結論する。

この鎖が独立に受理されれば、generic VでのCB full cone completionも

    C[[g20,s8,t²,s1*t,...,s8*t,t³]] ⊂ C[[g20,s8,t]]

となる。CAの一点の型からCBへ移した結論ではない。
本packet単体はfull geometryの監査前に全CB型を確定しない。
B18原点のfull odd ideal、三変数cubic、actual CA/CB intersectionの完成型も未解決である。

保存:

* `code/ambient_CB_V_actual_odd_linear_20261103_v63.py`
* `reports/AMBIENT_CB_V_ACTUAL_ODD_LINEAR_20261103_v63/{report.json,source_jets.json.gz,certificate.json.gz}`

再現:

    .venv/bin/python code/ambient_CB_V_actual_odd_linear_20261103_v63.py --output /tmp/cb-v-odd
