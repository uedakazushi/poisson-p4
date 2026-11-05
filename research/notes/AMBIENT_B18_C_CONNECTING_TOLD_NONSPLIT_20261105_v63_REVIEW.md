# 元 Told の UCT は R0 線形非分裂：独立純粋検証

独立 SOURCE peer は PASS した。報告 SHA256 は
`e4abd6f3315e8f9f05d039aa8b8eba86f96f9dae3213503809908ae697c23956`。
新しい係数値や source/jet を求めず、既存の実 HF と source 関係から
元の完全列 0→K/ellK→Told11→真正 C[ell]4→0 の非分裂を確認した。

ここで用いる r は元の rA であり、DVR parameter ell ではない。
実 Told の HF は (5,3,2,1)。gr¹ の基底を
u=rh3、v=z5h4、w=rh4、gr² の基底を
q1=z5rh4、q2=r²h4 と取る。これは単なる既知初期上界ではなく、
独立証明された正確な HF と literal十一基底により実加群の基底である。

ITold の四次元 QQ 基底は v、z5rh4、r²h4、r³h4。
後三つは次数二以上なので、gr¹ での像は厳密に v の直線だけとなる。
M7=Told/ITold の実関係 rh0=rh1=0、rh2=rh3/4 から、元 Told では
五 head の r 倍の初期式はそれぞれ

  alpha0 v、alpha1 v、u/4+alpha2 v、u、w

となる。三つの alpha は未計算のまま保持した。
既存の実初期関係 G38 は r²h3+(11767/3652)q1=0 を gr² に与える。
従って全五 head に対する r² の gr² q2座標を別に計算すると

  (0,0,0,0,1)

となった。

真正 target の cyclic top g を chi(g)=1 と正規化する。
実 connecting の head chi は (0,0,0,0,d)、d=−36/1681 である。
任意の持上げ x について、定数 head4係数は1/dに強制される。
残りの mTold 部分を r²倍しても m³Told に入るので、r²x の gr² q2座標は
常に1/d=−1681/36≠0である。これは明らかな head4 持上げだけでなく、
**全ての実持上げ**についての障害である。

もし R0 線形 section があれば r²s(g)=s(r²g)=0でなければならない。
真正 target では r²=0であるため、上の実 gr² 非零と矛盾する。
よって元 UCT 完全列の R0 extension class は非零である。

障害 q2 は ITold に入り、商 M7 では消える。
従って既に独立確認した M7 内の S による分裂と矛盾しない。
この非分裂から ITold→Kbar7 の別 extension が分裂するか否か、Kbar の全型、
元 H0 の長さ・指数、nil 三乗は決めない。

実行0.0123秒、所有上界11,741で完了した。
再現 checker は `code/ambient_B18_C_connecting_Told_nonsplit_review_20261105_v63.py`。
固定した小さい形式行列と source metadata のみを使い、凍結出力先への再実行は拒否する。
別作業コピーで再検証すること。
