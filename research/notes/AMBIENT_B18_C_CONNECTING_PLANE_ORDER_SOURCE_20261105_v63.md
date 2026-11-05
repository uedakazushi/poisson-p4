# 五 head 接続写像の平面次数：sigma3 は不要、sigma2 の一 adjoint は保持

前段 FIVE_FORMAL の168個の実 source reference は正しいが、その注記の
「第八生成元 f63」という表現は二つの source frame を混同していた。
honest C0 の U9 は実際には zA0,…,zA6 と零断面への制限 gA505938 を使う。
一方、元の640関係の Q0 multiplier は物理係数 f63=7790 を使い、
厳密な telescope f63=gA+sum(e_i zAi) で前者へ移す。この区別を source ID と
既存 telescope の全七係数で固定する。凍結済みの参照、符号、sigma、RegDivell
の式は変更せず、前段の注記だけを本 companion で訂正する。

Fj=delta+bA*tauA−bB*tauB と二つのモニック二次式から、四 basis 倍の
三つの奇成分を直接展開する。既存 C2 と pure-ell cache、および実 first-frame
証明により、全12物理係数の定数・平面一次は零、delta と bA の平面二次も零。
nuA,nuB の定数・平面一次も零である。U9 の奇成分にも同じ次数条件を確認する。
したがって全56列で V0=V1_plane=0、さらに V2_A,plane=0 が成り立つ。

ただし V2_B,plane は零ではない。物理18,60,67の Fj そのもの、列0,16,20に
既知の bB 平面二次が残る。高い B 座標の分子は
N_(r²ell²)−(96/275)N_(r³ell) なので、
(V2_B,ell²−(96/275)V2_B,r ell)*sigma2_rr を保持する。
これを単に次数だけで消すことはできない。

q7 が要求する各分子 monomial の全積経路を列挙すると、ell に依存しない
sigma のうち必要なのは sigma0、sigma1_r、sigma1_z0/z4/z5、sigma2_rr のみ。
sigma3 は実一次平面係数の零によって厳密に落ちる。混合方向 sigma1_z は
B の二次平面項から生じるため省かない。

V4 の二つの平面方向への contraction、V3*sigma1、および残った sigma2 adjoint
はまだ未評価である。新しい高次 source/IFT/gauge/jet、元 connecting、H0 や
nil 三乗の値を計算したとはしない。
