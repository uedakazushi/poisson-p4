# 五つの実 core 三次式へ向けた adjoint/dual preflight

2026-11-04。これは source/theorem と有限入力の監査のみ。
三次・四次係数の新評価、source 三次・四次解、GB はまだ実行しない。

## 1. 対象と零勾配

実98 unit 方程式 E(u,c)=0 と、その元 scalar DAG を保持する。
u は96係数未知数と二モニック消去座標。S の九パラメータを
五つの消去 tangent 座標 U と四 core 座標 C に分ける。

    U=(x108,x109,x113,x118,x119),
    C=(x112,x117,tauA,tauB).

残る五関数の元出力差 Fi は rows63,64,68,73,74。
修正 row72-row60 はこの五関数には含まれない。
独自 reverse-linear AD で、これら Fi の全98 unknown に関する原点勾配が零であることを確認する。
また既受理の実二次 core packet で、五つ全ての Gi,2|U=0 が零であることも確認する。

## 2. source4 が消え、source3 は linear adjoint で消去できる

保存二次近似を u_[2](C)、未知の真の三次・四次補正を u3,u4 と書く。
J0=E_u(0) は同じ実98 unit 行列。E3 は u_[2] を元方程式へ代入した三次残差である。

    J0*u3=-E3.

Fi_u(0)=0 なので未知 u4 は Fi の四次項に寄与しない。
u3 の寄与は Fi_u(u1,C) の一次部分だけである。
この5x98の linear-in-C 行列を Hi,1(C) とし、

    lambda_i,1(C)*J0=Hi,1(C)

を解く。この恒等式を全5*4*98係数で直接確認する。
従って真の四次出力は厳密に

    Fi4 = [Fi(u_[2],C)]4 - lambda_i,1(C)*E3.

未知 u3 を零と置いた計算ではない。u4 の情報も不要。
Hi,1 は元 source の symbolic Hessian 全体を展開せず、
保存 u1 で評価した scalar DAG に reverse AD を行い、四変数一次 ring の係数として得る。

## 3. 五 IFT の二次加速度と dual 方向

五つの選択関数は rows18,19,60,28,29 で、G の一次行列 D5 は対角
(1/12,1/12,1/24,1/12,1/12)。
これらの raw Fi の原点勾配には零を仮定しない。定数 adjoint lambda_j,0 を別に作る。

    U2_j=-D5_jj^(-1) * ([Fj(u_[2])]3-lambda_j,0*E3)/L1.

ここでは U=0 とし、L1=5*x112-20*x117+24*tauA+6*tauB。
G選択関数の一次項が U のみなので、L2 の補正は消える。

五つの残余 Fi3 の U 微分は、保存済み**全九変数 source 二次**を
四 core 変数と五つの独立な first-order dual U 方向へ代入すれば得られる。
未知 source u3 の dual 係数は Fi_u(0)=0 を掛けるので消える。
従って source3 の新求解を行わない。
dual ring で E が core 次数二まで零だとは仮定しない。
保存二次近似が保証するのは E の base 次数二以下、および dual 係数の core 次数一以下の零である。

五 IFT 後の真の core 三次関数は

    Ki3 = (Fi4|U0 + sum_j [(partial_Uj Fi3)|U0]*U2_j)/L1.

Gi1=0 と Gi2|U0=0 により、この式にも L2 は現れない。
Npair の新しい高次係数も不要である。

## 4. 有限入力と budget

二族の保存 source57二次を、既存物理 R0 逆写像と二つのモニック Schur 二次写像へ
実際に合成し、同じ全九パラメータで98未知数を bind する。
そこから四変数二次と五 dual 方向の入力を literal に保存する。
選択五つの constant adjoint と残余五つの linear adjoint を合わせ、
全25*98=2450個の QQ 恒等式を確認する。

提案する一回係数評価の ring は、四変数四次以下70単項式の base と、
五つの dual 方向それぞれ core 二次以下15単項式、計145係数/node。
base 乗法495、各 dual 項の乗法45 ordered pairs。
必要 DAG 祖先・乗法・unit inverse の数と conservative scalar-pair budget を report に固定する。
2GB/60秒/32MiB上限で止める条件を保ち、次数や求解を自動的に増やさない。

次の有限完了条件は、五つの Ki3、全70係数での L1 除法恒等式、
source 98残差の base<=2/dual<=1 零、および adjoint correction を保存すること。
さらに既受理 nodal quartic の四つの binary quartics Phi へ Ki3 を代入すれば、
binary degree12 の厳密な section になる。
その一つが非零なら、(f,g) の初期 ideal から Ki は実二関数の整域曲面 S2 内で非零。
従って full core の次元は高々1へ下がる。
全て零だった場合は membership や socle を肯定せず、そこでこの方法を止める。

今回の preflight はこの係数評価をまだ行わない。
