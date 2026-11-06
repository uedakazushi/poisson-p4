# 階数4接二次曲面のanisotropic代表：MPS閉包への独立確認

2026-11-07（Asia/Tokyo）。対象は次の一つの射影Poisson軌道である。
この軌道が既知MPS40成分の和に含まれることを、二つの実PGL極限と全係数恒等式で確認した。
階数4接二次曲面層の全点がこの軌道または別のisotropic族に帰着することは、別の正常形証明への依存である。
新しい成分、一般ファイバーの尽くし、完全分類を主張しない。

## 1. 全射影Jacobiを満たす実際のpencil

座標を(h,u,v,z,w)とし、hを中心とする代表Pの六係数を

    P_uv=hu,                 P_uz=-hw-3u^2-v^2,
    P_uw=-2uv,               P_vz=-2hz-2uv,
    P_vw=-3hw-u^2-3v^2,      P_zw=4uw-4vz

とする。ほかに次の七係数を持つ二次bivector Rを取る。

    R_hu=-hu/2,    R_hv=-3hv/2,   R_hz=hz/2,   R_hw=3v^2/2,
    R_uv=-uv/2,    R_uz=uz/4,     R_vz=vz.

書かれていない係数は零、逆順の添字は反対符号とする。
Euler場をEとすると、全パラメータtについて

    E wedge [P+tR,P+tR]=0

が成り立つ。Rにはh方向の項があるため、h-Casimirの制限付き変形とは異なる。
ここで証明するのは射影Jacobiであり、[P+tR,P+tR]そのものの零を必要としていない。
退化五次式は正確に

    E wedge (P+tR)^2/2
      = -3/8 (tuz+4hw+4u^2+4v^2)
              (thvw+tv^3+2h^2w+2hu^2+2hv^2) partial_(huvzw).

特にt=0では-3h(hw+u^2+v^2)^2となり、指定されたanisotropic代表の退化因子を戻す。
一般tでの二次・三次因子の存在だけから、MPS所属を推論しない。

## 2. P+Rから五つの独立超平面へのPGL極限

整数重みa=(1,0,0,-3,-2)について、二次bivectorへの作用を

    (rho_lambda A)_ij(x)
       =lambda^(-a_i-a_j) A_ij(lambda^(a_0)x_0,...,lambda^(a_4)x_4)

と定める。lambda≠0では実際の線形座標変換、従ってPGL5作用である。
P+Rの全単項式の重みは非負であり、零重み部分はRからR_hwだけを除いたR_logである。
R_logの各非零係数は対応する二座標の積なのでlogcanonicalであり、

    E wedge R_log^2/2 = -3/8 h u v z w partial_(huvzw).

従ってその退化因子は五つの独立座標超平面の被約SNC因子である。
非零lambdaの族はP+RのPGL軌道に属し、lambda=0でR_logへ延びる。

P+Rを含む任意の全Poisson既約成分Cを取る。PGL5は連結なのでCを個別に保ち、
閉性からR_logもCに入る。MPS Theorem1.4（既存40図式との同定を含む）により、Cは既知MPS40成分の一つである。
この論法の向きは「元の成分がSNC極限を含む」ことであり、既知集合への極限の所属だけから元点の所属を推測していない。

MPSの適用範囲と一次資料の依存は `BETA_ZERO_COVERAGE.md` および
`LOWER_DEGREE_SNC_COVERAGE.md` と同じである。外部定理自体を今回再証明したとはしない。

## 3. 元のanisotropic軌道への収縮

別の重みb=(0,1,1,2,2)では、Pの全項の重みは-1、Rの全項の重みは0である。
従って

    lambda*rho_lambda(P+R)=P+lambda R.

射影Poisson空間では非零定数倍は同一点を表すので、lambda≠0ではP+Rの軌道にあり、
lambda=0ではPに延びる。前節で得た既知MPS成分は閉じているから、Pおよびその全PGL軌道も入る。
この二段階の極限により、anisotropic代表の既知成分被覆が成立する。

元のPを含む全ての成分がMPSであるとまでは主張しない。
支持の被覆に必要なのは、Pが既知閉成分の少なくとも一つに入ることである。

## 4. 独立計算の範囲

`reports/CLASSIFICATION_ANISOTROPIC_MPS_DEFORMATION_REVIEW_20261107.json` に、
再実行可能なchecker全文・hash・入力hash・全P/R係数・重みを保存した。
producerとSymPyをimportせず、Fractionによる疎多項式積・偏微分・外積から
全パラメータ射影Jacobi、五次式、二つの重み、SNC極限を独立に再構成した。
有限契約のgateを通過し、CPU30秒・仮想メモリ1GiBの下で0.021秒で完了した。

SNC極限の対数係数行列をw=1上で反転し、biresidue行列

    [[ 0, 4,-1,-2,-1],
     [-4, 0, 2, 6,-4],
     [ 1,-2, 0,-2, 3],
     [ 2,-6, 2, 0, 2],
     [ 1, 4,-3,-2, 0]]

を得た。smoothable辺は04、24、13であり、最大装飾図式は置換
(0,1,2,3,4)→(2,4,1,3,0)でT33に一致する。
これはSNC点の図式照合であって、それだけでPの所属先が一意にT33だとは主張しない。
上のMPS40の和への包含には、この追加の一意同定は不要である。
