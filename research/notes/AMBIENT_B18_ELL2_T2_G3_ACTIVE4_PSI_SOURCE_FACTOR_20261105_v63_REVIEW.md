# 残る active 四方向・八つの形式 psi contraction の独立 source 監査

固定 projected rank25 証明の active16列は plane条件 rank20 証明の選択列に含まれる。全238 active 列の保存 factor が最初16 basis だけを用いることを確認した。従って active projected kernel は次元222であり、plane kernel218に対する補方向は77,101,117,125の四つとなる。各方向の全102×9一次項を、保存 selected99 定数列で独立再構成した。新しい rank 探索をしていない。

未知 post27 Ystar 一次微分を独立記号として保持し、全64152形式 slot と145個の非零 label を照合した。既知上側一次商の自由座標25,26は top2*r、top2*ell である。その商を使うと各方向の寄与は `−NF(top2*r²) A−NF(top2*r*ell) B` に因数分解される。A は b1_r*psi0 の和、B は b1_ell*psi0+b1_r*psi1 の和である。全八つの重みを照合したが、値は求めていない。

selected99は全て元741個の定数零 source 結合である。元768 X0で各結合が零であること、native S0だけの結合も別に零であることを確認した。ただしこの事実だけで新しい重みにおける old honest second 感度が消えるとはしない。

**post27 の carry。** `Ystar_j=Yj−Yp*a_j`、a_j0=0である。G2Yj=MXjを二次まで展開し、`G2_1Yp0=M0Xp1+M1Xp0` を使うと厳密に

`G2_1Ystar_j1=M0(Xj2−Xp1*a1_j)+M1(Xj1−Xp0*a1_j)−G2_2Yj0`

となる。M0Xp0=0により未知a2が落ちる。未知Yp0を零とせず、この等式でその carry を保持している。Xj0=0によりM2項もない。

**coherent frame。** eps_M(u,b)=n_M G2u+Hfree_M b と書く。M'=M+G2Aに対し、`U_A(u,b)=(u+A eps_M(u,b),b)` は形式単元であり、eps_M' U_A=eps_Mを満たす。さらにU_Aはpure G3列を厳密に固定する。選んだ contracting homotopy による別のliftとの差はG3境界であり、R線形性を仮定しない。

27と99のunit消去後、定数frameの残りは `top += PiY0 A0 Hfree0 E_free3 bottom` という三角変換である。従ってcanonical psiだけをnormalized L48のnormal formへ代入してはならない。normalized M'_1=M1+G2_1A0と対応する全post27 carryを使うか、top/bottom二次関係とcoverをこのtyped三角frameで共に運ぶ必要がある。canonical残り二muとnormalizer A0は未計算のまま保持する。G3 liftの一次曖昧さはPiTop2が消し、再消去で生じる差も既知一次関係の一次倍として二次normal formで消える。

この監査はsource/formulaのみを受理する。新しい八重みでのsource availability・honest second感度・normalizer依存の閉鎖は別gateである。実psi、HF2、全長は未決であり、新しいphysical/gauge/source係数は求めていない。
