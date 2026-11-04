# Actual core/(d) の DVR 模型と一有限部分の独立監査

凍結 NPAIR_RESTRICTION_DVR を別実装の sparse QQ 演算だけで検査する。まず元の
actual12-flag-divided-polynomials.json の十二 Hi を直接 d=0 に制限し、Conn の八変数
入力と全係数比較する。そこで初めて五つの定数 unit graph relations と全七残余式を
再生する。作者の recorded H12 だけを前提にしない。新 GB・source jet は実行しない。

u=a+v2,x=v2,t=v6 の七式が F=u²-6ux,G=(23x-6t)u+2x²+4xt,H=ux² に属する全七
cofactor identities と、最初の三式から F,G,H を復元する恒等式を検証する。
q=ux に対する u*q,x*q,t*q の三つの source cofactor identities も確認する。

QQ[t]^3⊕QQ[t]/(t)q 上の literal table を読み、別演算で全64 associativity checks、
全source関係の零性、q≠0を確認する。逆方向の写像も、table の各関係が F,G,H と
三 annihilator identities から成り立つことで well-defined である。両写像は u,x と
q=ux を固定し、従って同型。これは単なる source からの全射模型ではない。

有限長部分は q の一方向に等しい。有限長元は t の冪で消えるので free QQ[t]^3 部分を
持たず、q は u,x,t 全てに消される。従って actual core/(d) の H0長さは1、深さ0。
graded basis degrees0,1,1,2から Hilbert series は (1+2z)/(1-z)+z²、重複度3となる。

t saturation は q を除いて (u²,ux,x²-3tu+2tx)。この free rank3 商から、double line
(3u-2x,x²) と single line (u,x+2t) への写像を検査する。自由基底上の gluing matrix
determinant は4t/3≠0なので単射であり、両 ideal の intersection と正確に一致する。

対象は actual divided core を d で割った環である。E4/(d)、D0、full D の深さには
この H0結果を移さない。結果・時間・入力 hash は report と independent-proof に保存する。
