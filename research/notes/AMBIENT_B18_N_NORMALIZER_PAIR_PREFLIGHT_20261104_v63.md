# B18: 同じ N 正規化源上の二つの実トレース根

2026-11-04。新しい一次 packet。凍結済み source・原稿・台帳は変更しない。

## 1. 範囲と有限入力

係数体は QQ、必要なら C へ拡大する。全て B18 の指定点での形式局所環である。
`N_SAME_R0_BINDING` の実正規化源を B_N と書く。これは正則八次元環であり、
物理的な同じ九座標の順序は

    (rho27,rho108,rho109,rho112,rho113,rho117,rho118,rho119,rho61)

である。七つの一次座標

    X=(rho108,rho109,rho112,rho113,rho117,rho118,rho119)

とその packet の元 source 変数 p_N が B_N の形式座標となる。
source28・実20切断・56係数 unit 系は、元の有限 DAG と既存の独立レビューをそのまま用いる。

両成分の零 weight 正規化は、同じ R0 九座標上で

    B_i=R0[tau_i]/(tau_i^2-nu_i),
    physical_q_j=a_ij(R0)+b_ij(R0)*tau_i

と書かれている。a,b,nu は `ACTUAL_TOP_COEFFICIENT_PREFLIGHT` の実48係数
unit 系の根であり、新しい任意の正規化座標ではない。
`ZERO_NORMALIZER_QUADRATIC` は、その関数の既存二次キャッシュを同じ物理 R0 へ変換したもの。
今回の検算は N の保存済み17単項式の mixed 二次係数を再評価するだけで、新しい source 求解を行わない。

N の w=1 tail の同じ物理零 weight 出力は零 weight N へ写る。
この点で用いるのは、既に source-bound された有限 w 多項式模型の正次数座標を
lambda の各重みで掛ける代数的族である。N は閉で Gm 不変なので lambda=0 も N に入る。
従ってこの零 weight 出力は実 CA,CB の双方へ写る。w=0 完成環で w=1 を代入する操作ではない。
この零 weight 射は、以下の非零係数が示す通り、両成分の conductor 全体には入らない。

## 2. 固定旗での補助的な単位変換

CA の生表示を

    P=f(h,x) dx^dy + g*h^2 dz^dw + dy^X_H,
    f=h*(A*h+B*x)

とする。ここで wedge 記号は上付きの ^ で略記し、H は許された全16個の三次係数を持つ。
A が単位の開では h'=A*h+B*x, x'=x に変える。h=(h'-B*x')/A とすると、
pushforward に Euler vertical

    B*h * (dy wedge Euler)

を加えた生テンソルは CB の表示に正確に等しい。その係数は

    f_CB=A*h^2,   mu_CB=g/A,   H_CB=H((h'-B*x')/A,x',z,w),
    -(partial f_CB/partial x')/2=B*h.

全25生行列成分を全16係数のまま検算した。反転は A だけであり、B,g,conductor 係数は反転しない。
B18 の交換済み固定旗中心では A=-1 が単位である。
この式は具体的な proper-source lift を与えるが、以下の同じ R0 における根の構成は
この旗の名前を現在の core p_i と同一視せず、実物理出力から独立に行う。

## 3. 消える分母を単位扱いしない正則根

両側とも物理 row18 の実 b の B_N への制限は非零である。一次項は

    b_A18|N = 2*X108 + O(m^2),
    b_B18|N = -X108/2 + O(m^2).

この非零関数を f_i とする。実成分の像 A_i は R0 と全 a_ij+b_ij*tau_i が生成する環で、
f_i*tau_i は A_i にある。f_i を反転すると tau_i=(f_i*tau_i)/f_i が回復するので
A_i[1/f_i]=B_i[1/f_i]。実 N 零 weight 射をこの開へ制限すれば、同じ分数体 K(B_N) に

    t_i=(physical_q18|N-a_i18(R0|N))/b_i18(R0|N)

が定義され、t_i^2=nu_i(R0|N) を満たす。B_N は正則、従って整閉なので t_i は B_N に属する。
これは非優越な N→CA/CB に対する正規化の普遍性の誤用ではない。
まず非零の実 f_i 上で分数体射を作り、その後にモニック関係と整閉性を使った。
中心で f_i を単位として反転しない。この式は一意な正則商を指定する有限代数的関数入力である。

nu_i(0)=0 より t_i(0)=0。元の48係数 unit 系の正則根とこの一意な正則商は、
具体的な full source のトレース根を同じ B_N 上に結ぶ。

## 4. 実一次写像と hypersurface

mixed 二次キャッシュで、q18-a18 の X108*p_N 係数を上の b の一次係数で割り、根の符号を決める。
nu_i の p_N^2 係数と X_j*p_N 係数から全一次項を得る。結果は

    t_A = -X112/8 + X117/6 - p_N/2 + O(m^2),
    t_B = -X112/3 + 8*X117/3 + 2*p_N + O(m^2).

12出力全てで mixed 係数の整合を検算する。保存 source28 の全20切断も同じ17係数で零。
E4 の二つのモニック方程式は既存の unit minor -20/9 により rho27,rho61 を消す。
残る正則九次元環 S の座標は (X,t_A,t_B) である。
写像 S→B_N のヤコビ行列は rank8、(X,t_A) の minor は -1/2。
形式逆関数定理で S→B_N は全射で、その核は一つの正則元 L_formal で生成される。
一意な一次項の規格化を

    L_formal = 5*X112 -20*X117 +24*t_A +6*t_B + O(m^2)

とする。この hypersurface は実 N-normalizer pair の閉像であり、単なる二次共通因子ではない。

## 5. 全次数の因子と未計算の部分

二つのモニック式を消した後の12個の実出力差 F_j は、同じ N の物理出力なので B_N 上で零。
従って各 F_j は厳密に L_formal*G_j と一意に書ける。
この正則除法は S の形式座標 (X,t_A,L_formal) における一変数 Hadamard 除法であり、
L_formal を局所単位として反転する意味ではない。
全12個の G_j は実入力から一意に指定されるが、この一次 packet では高次係数を展開しない。

E4 の全環、四 core 方程式、socle、D0 の深さ、T0 の有限部分は決めていない。
正規化 tensor 比較 kernel の議論とも区別する。新 source jet、GB、beta 点は一切追加しない。

## 6. 再現と依存

    .venv/bin/python code/ambient_B18_N_normalizer_pair_preflight_20261104_v63.py --output work/B18-N-pair-replay

主要入力は出力 report の SHA 表に固定する。加えて零 weight 射の幾何は
`AMBIENT_B18_W1_TAIL_BINDING_20261104_v63_REVIEW` と
`AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63` の有限 w の範囲、
および全 proper N 像を扱う `AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63_REVIEW`
を用いる。後者の rank<=2 境界の特定だけを新しい正規化普遍性へ置き換えてはいない。
今回の非零 pullback と minor は別添 evidence の全有理係数で確認できる。
