# generic B18でのCB正規化射：全ratio×fを含むcritical scheme

2026-11-03。同じB18代表P0で、CB normalizer sourceの全ramification Fitting schemeを決定する。
target CBやactual CA/CB交叉Dの全完成環を求めたという意味ではない。

## 全sourceとassociated ratio

旗座標(h_N,k_N,y,z,w)=(x,h,y,z,w)で、fixed sourceを

    P(f,q,H)=canonical(raw_CB(f,−q*f,H))

とする。fはbinary3係数、Hは16係数、qはassociated P1のendpoint近傍のratioである。
中心はf=−h_N²、q=0、H=−h_N*k_N*w+k_N²*z+z³/3である。
source fiber20微分には、∂fの列F_i−qG_iと∂qの列−Σf_iG_iをすべて含める。
q*fのmixed項を落とした接線模型ではない。

旗の同じ16次元parabolicについて、enlarged22空間への全作用と、
ratio character χ=M33+M44−M00−M22を含む全19fiber列の安定性を検証する。
具体的にはembedding E(q)へ

    rho(M)E(q) −χ*q*E'(q)

を作用させたものが再びE(q)の像に入る。
これによりassociated ratioごとの実際のsource actionを使っていることが確認される。

P0でparabolic fiber orbitはrank10であり、CAと同じ十係数への射影にdet−2/3を持つ。
その十係数を固定したaffine10断面とP10局所切断の作用写像はetaleである。
旗9と合わせて19個のsmooth source因子を取り出す。

## exact全critical ideal

この断面上の可変150×20 fiber Jacobianの選択blockは恒等的にdet−864/25である。
それを除いた全150×9旗行列の8次blockは恒等的にdet8/3である。
従って残る列の全Schur係数は、全source Jacobianの最大minor idealと正確に等しい。
全70非零Schur多項式のidealをQQ上で計算する。

断面変数s0,...,s9について七つのmonic消去式は

    e5=10*(s3*s7+s7²)/27,
    s0=2*s3/9−16*s7/9,
    s1=−3*e5/5,  s4=27*e5/10−s7²,
    s5=e5, s6=0, s8=s7², s9=s7³/3.

残る変数をw=s2、u=s3+s7、v=s3−s7とすると、全idealは正確に

    (wu,wv,uv,u³)

となる。全Schur式とこの七消去式＋四生成元との両ideal包含をQQで検証し、
別のMacaulay2実行でもideal等号を確認する。

したがってcone normalizer ramification完成環は

    C[[g1,...,g19,w,u,v]]/(wu,wv,uv,u³)

であり、projective版はg変数18個である。
CAのsource critical完成環と抽象同型だが、有限sourceの変換・消去係数は異なる。
両者のtarget image idealが等しいと推論することはできない。

primary分解(u,v)∩(u,w)∩(w,v,u³)から、二つの被約枝とjunction embedded素を得る。
nilradical(u)は三乗零・二乗非零で、cone dim20/edim22/depth19、projective dim19/edim21/depth18、重複度2である。
全odd idealやactual Dのnilpotence・associated primesを決めた結果ではない。

* code: `code/ambient_B18_CB_critical_scheme_20261103_v63.py`
* evidence: `reports/AMBIENT_B18_CB_CRITICAL_SCHEME_20261103_v63/`
