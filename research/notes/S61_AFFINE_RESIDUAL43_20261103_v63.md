# 全未知残余の五determinant方程式と、一つのdegree43方程式

2026-11-03。[大域polarization splitter](S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63.md)
を変えずに、その残余idealを**一つのdegree43多項式**で正確に記述する。
これは全一般fiberの尽くしではない。未知の孤立点・正次元・非被約支持は全て残す。
未知q、cyclic determinant、quintic Hによる局所化は行わない。
一般係数多項式は有限五次代数の逆行列・行列式recipeで指定し、巨大な展開表を保存したものではない。

## 1. 五つの実determinant係数ideal

A=K[q]/I、Eはその実商である既知有限étale5代数、pはE係数の普遍点とする。
既存packetの45元行selectorで

    D_p(q)=det J_selected((q+p)/2), d(p)=D_p(p),
    e=Tr_E/K(D_p/d(p))

とする。既存証明からA=E×A_res、A_res=A/(e)である。
分離元h=Σ(j+1)q_jの五元基底で

    D_p(q)=D0(q)+hD1(q)+…+h^4D4(q)

と一意に書く。分裂後eは五つの直交idempotentの和なので、
**e*D_p=D_p**がA⊗Eで成り立つ。係数ごとにD_j=eD_jとなる。
逆にe=Σ_j Tr(h^j/d(p))*D_jだから、正確に

    I+(e)=I+(D0,D1,D2,D3,D4)

となる。determinantの根だけを使うradical等号ではなくideal等号である。
この五式は全未知支持と非被約構造を保つ。

## 2. 普遍三平面による次数低下

E1=span_K(1,p0,…,p44)⊂Eとする。
全16betaを保った普遍三平面の証明からdim E1≤4。
今回の中心mod32003の評価でrank4なので一般にはdim E1=4である。
参照：TORIC_UNIVERSAL_THREE_PLANE_20261101、S61_FIVE_COVER_PRODUCT_20261101。
有限体rankだけから一般のrank上界を推論していない。

一方Jはqに関して一次であり、

    M_p(q)=J_selected((q+p)/2)=L(q)+B(p)

と書ける。Lはq一次、Bはpに関してaffine一次である。
determinantのq次数45部分はE係数として1のspanに入り、
q次数44部分はE1に入る。
従ってE1を消す非零linear functional ell:E→Kを選べば

    R43(q)=ell(D_p(q))

のq次数は43以下になる。

ellを一般betaで指定するには、Eをhの五元基底で表し、
lambda4=1として四つの独立な関数(1と選択p座標)を消す4×4系を解く。
その全selectorと非零minorを保存した。
残る全p座標の消失は、既存普遍三平面によるrank4上界と四独立条件から従う。
従ってこれは一般beta上で定義された有限逆行列recipeである。

中心mod32003では

    ell=(0,0,22810,21387,1)  [h^0,…,h^4係数に作用],
    splitweights=(12072,12072,12072,7859,19931)。

固定したq方向でq次数45,44,43の係数を行列式のtrace公式から直接求めると

    (0,0,15015)

となった。よってR43は零多項式ではなく、一般betaで**正確に次数43**である。
raw trace idempotent eの同じ方向の次数45係数は10048≠0であり、
traceを取るだけの平均化では次数は下がらない。

## 3. R43一式が残余ideal全体を与える

Eを分裂し、各点p_iのLagrange idempotentをL_iとする。
D_p(p_i)=d_i L_iなので、r=phi(R43)∈Eの各値はd_i*ell(L_i)である。
中心での五値は

    30124,15759,23800,4157,26406、 norm=17591≠0 mod32003。

従ってbetaだけの非空開でrはEの単元である。
R43=eR43は既知因子だけに支持し、その既知因子上では単元だから、
A内で(R43)=(e)となる。明示的には既存写像sigma:E→eAで
R43*sigma(r^(-1))=eである。
従って全未知残余の実idealは

    **I_res=I+(R43)=I+(D0,D1,D2,D3,D4)=I+(e)**。

この等号はスキームとして成立し、未知残余の非被約部分も落とさない。
betaだけの追加normと4×4minorは全て中心でp単元である。
五式を残す場合は、D4をR43に置き換えた生成系の次数上界は
(45,44,44,44,43)になる。

## 4. 有限検査と次の実際の目標

producer `code/s61_affine_residual43_20261103_v63.py`、
receipt `reports/S61_AFFINE_RESIDUAL43_20261103_v63/`。
保存した五次代数の座標、四次小行列式、全ell条件、既知点でのD係数行列、
五つのr値とnorm、trace公式による最高三次数を独立に再計算できる。
大域splitterの全134QQ polarizationと25 midpoint検査は凍結先行packetを用い、hashを再確認した。

一般fiber尽くしの新しい具体的目標は、**全元134式とR43の生成idealがunitであること**。
そのunitは未取得である。次元零・長さ5・成分上界をこの変形から結論しない。
