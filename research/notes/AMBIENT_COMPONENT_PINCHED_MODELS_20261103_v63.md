# 指定 Q における被約 CA/CB/N の全完成局所環と conductor

2026年11月3日。これは **個々の被約成分** の完成局所環を決める定理である。
全 S61、全 CA∩CB の scheme-theoretic intersection ideal、他の枝や埋込み構造は決めない。
指定 Q は v24 から継続している a=1 の点で、a=2/9 の別の共鳴点とは区別する。

## 定理

    Ac=C[[u1,…,u10,s1,…,sc,p²,p³,s1*p,…,sc*p]]

と置く。S61 の被約成分の指定 Q における完成局所環は

    Ohat_(CA,Q)≅A8,    Ohat_(CB,Q)≅A8,    Ohat_(N,Q)≅A7

である。同じ記述は下記の7次元臨界 rank2 一次ジェット軌道上の全ての非正規点に適用される。

|対象|次元|embedding dimension|depth|multiplicity|
|---|---:|---:|---:|---:|
|CA,CB (S61)|19|28|11|2|
|N (S61)|18|26|11|2|

Ac の正規化は

    Bc=C[[u1,…,u10,s1,…,sc,p]],

正規化商は `Bc/Ac≅C[[u1,…,u10]]*p` で、conductor は Bc 側で

    (s1,…,sc,p²)

である。三つの成分の被約非正規層・被約特異層は、この近傍では共通の滑らかな10次元部分である。
全射影模型の値・perspective の滑らかな九方向を戻すと、次元は28,28,27、
embedding dimension は37,37,35、depth はいずれも20となる。

この環は通常 cusp product ではない。既存の multiplicity2、analytical irreducibility、R1、非 CM 性を
全局所環の表示で具体化している。

## 1. 有限正規化と全 source slice

入力は閉 proper birational な滑らかな解消 EA,EB,EN と、
基底像 Z の7次元臨界軌道に沿う全 slice 定理
`AMBIENT_BASE_PINCHED_SURFACE_20261103_v63.md` である。
基底点を

    b0=(τ,A,B,C,D,F)=(0,0,1,1,0,0)

とする。この点の全基底旗ファイバーは一つの台を持つ長さ2なので、各成分の旗ファイバーは
その部分スキームとなり、いずれも長さ1または2の一点である。
従ってこの基底軌道上の全点の近傍で、各成分の解消は有限となり、正規化である。

基底 slice の証明で使った七つの実際の β 座標を固定する。
Sp4 同変性と étale な H7×T9→V16 により、各成分とその正規化は形式的に七次元 smooth factor と
この slice 上の成分・正規化の積となる。
正規化側は基底正規化の二変数(s,p)と、九または十の affine ファイバー変数の全形式的 product である。
基底の任意の linear section と scheme-image を交換していない。
ここでも七次元群作用の厳密な product 分解を使う。

slice の全 source は

    τ=p, A=0, B=1, C=1, D=−p/2, F=s,
    旗=(0:1:0:p)

である。N では固定線 CA の r=0 と九つの三次係数k0,…,k8を動かす。
CA,CB ではさらにrを動かす。
三次式の順は

    G=k0*x²z+k1*x²w+k2*xz²+k3*xzw+k4*xw²
                      +k5*z³+k6*z²w+k7*zw²+k8*w³.

この全11または12変数の source を、実際の有限 symplectic 作用と perspective 正規化で61座標へ送る。
原134式の全て、全150係数の再構成、projective Jacobi の全係数は恒等的に零である。
特定の整数点での照合や有限体の標本確認ではない。

## 2. 実際の target 座標を正規化の regular 座標に選ぶ

この source の β の実座標は、既に

    s, p², p³, s*p

を含む。残る45個の q はファイバー係数に affine linear である。
p=s=0 におけるファイバー方向の係数行列 M0 はNで階数9、CA/CBで階数10である。
従ってその階数個の **実際の target q 座標** を K と選ぶ。
選択行は45個の q 内の0始まり indexで

    N:     (19,21,24,25,27,28,32,33,34),
    CA,CB: (0,19,21,24,25,27,28,32,33,34)

である。選択部分 L は定数可逆である。
K,s,p は正規化の全 regular system of parameters を与えるので、正規化環を

    B=R[[p]],    R=C[[s,K−K(Q)]]

と書ける。target 完成局所環 A の像は R と p²,p³ を含む。
従って **p²B⊂A** である。これが以下の一次奇数係数の計算を全局所環の計算へ上げる鍵である。
p 次数2以上の未知項は全て既に A の中にあり、切り捨て誤差にはならない。

## 3. 欠ける一次係数の ideal を正確に計算する

p=s=0 の q を q0、p 微分を q1と書く。
選択したKを固定して p を動かすと、元のファイバー変数は一次で

    −L^−1*(q1)_selected

だけ補正される。従って45個全ての q の一次 p 係数は厳密に

    q1−M0*L^−1*(q1)_selected

となる。β が既に sp を含むので、この係数 ideal は s を法として計算すれば十分である。
K による regular 座標の定数 affine 変換で、p=s=0 の残差ファイバー座標を
bar-r,bar-k0,…,bar-k8 と書く。この変換は target 座標だけを使い、p を使わない。

保存した全係数計算から、上の45個の一次係数が生成する ideal は

    N:     (bar-k1,bar-k3,bar-k4,bar-k6,bar-k7,bar-k8),
    CA,CB: (bar-r,bar-k1,bar-k3,bar-k4,bar-k6,bar-k7,bar-k8)

と正確に一致する。係数行列の階数は6,7,7で、定数項は零である。
よって s を加えた欠損 ideal は

    IN=(s,bar-k1,bar-k3,bar-k4,bar-k6,bar-k7,bar-k8),
    IA=IB=(s,bar-r,bar-k1,bar-k3,bar-k4,bar-k6,bar-k7,bar-k8).

各 ideal は regular 座標の一部であり、codimension は7,8,8である。
残る三つのファイバー座標は bar-k0,bar-k2,bar-k5 である。

一般に A⊂R[[p]] が R と p²B を含むなら、A は一次係数の ideal I によって

    A=R+I*p+p²B

と **全体が** 決まる。実際、右辺への包含は生成函数の一次係数から従い、
逆向きは I の生成元を一次係数に持つ target 函数の p²以上を取り除けば得られる。
像は有限部分加群として閉じているため、形式的無限和を取っても包含は保たれる。
従って高次未検査項は残らず、上の I が全 target 完成局所環を決める。

非正規点では I が極大 ideal に含まれ、七次元群因子と三つの自由ファイバー係数を u1,…,u10 と取ると、
定理の Ac が得られる。I が単位を含む点では A=B となり、その成分は滑らかである。

## 4. 指定 Q への適用

指定 Q の基底は (τ,A,B,C,D,F)=(0,0,4/3,3/2,0,0)である。
これは上の b0 と同じ Sp4 軌道に入る。
実際、x'=a x,y'=a^−1 y,z'=b z,w'=b^−1 w の torus 作用で

    C'=C/(a*b²),    B'=B*b/a²

となり、b^5=C²/B、a=C/b²を選べば C'=B'=1である。
必要な非零五乗根は C 上で存在する。
既に独立に検証した指定 Q の長さ2旗ファイバーは、この変換後も長さ2である。
したがってその点では上の欠損 ideal の全生成元が零であり、非正規 local model Ac が適用される。
具体的な五乗根を有限体へ落としたり、そのような計算に依存したものではない。

## 5. conductor と数値不変量

R=C[[u10,s1,…,sc]]、I=(s1,…,sc)、B=R[[p]] とすると

    B/A ≅ (R/I)*p ≅ C[[u10]]*p.

conductor はB内で(I,p²)、A内では(I,p²,p³,I*p)である。
これの被約零集合は滑らかな10次元部分となる。その外でA=Bは正則であり、
その上ではB/Aが非零なのでAは正規でなく、従って正則でもない。
ゆえにこれが被約非正規層・被約特異層の双方である。

完全列0→A→B→B/A→0で depth B=d、depth(B/A)=10、d=18または19なので
depth A=11となる。
regular normalization で m_A B=(u1,…,u10,s1,…,sc,p²) はcolength2 parameter ideal であり、
有限 birational rank1 の multiplicity 加法性から e(A)=2となる。
極大 ideal の独立生成元はd−1個の regular base座標、p²,p³とc個のsi*pなので
embedding dimension はd+1+cとなる。

## 証明書と範囲

* producer: `code/ambient_component_pinched_models_20261103_v63.py`
* 全61出力、target q 選択、可逆行列、全45一次係数と ideal の階数・零化条件:
  `reports/AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63/report.json`
* 三族全てで元134式、150係数再構成、全 projective Jacobi を QQ 上で検査。測定25.61秒。

個々の被約 CA/CB/N の指定点での局所環はこの定理で決まる。
全 CA∩CB の scheme-theoretic intersection ideal、全 S61 の完成局所環・他の枝・primary 構造は未解決である。
