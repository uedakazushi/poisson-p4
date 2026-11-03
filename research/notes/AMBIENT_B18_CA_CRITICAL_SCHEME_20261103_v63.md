# generic B18でのCA正規化射：全ramification schemeの完成環

2026-11-03。対象は、既存のproper旗模型から得たCAの正規化射の**source上のramification scheme**である。
実CA/CB交叉Dの全完成環やCAのtarget完成環を、この結果だけで同定しない。

## 点と有限正規化への依存

座標は(h,x,y,z,w)、L=<∂y>とし、generic B18の代表を

    P0 = ∂y ∧ (x²∂h + h*x∂z + (h²+z²)∂w)

とする。CA/CB/Nの全proper旗fiberの双対数長2と、有限近傍でのnormalization解釈は
受理済みの`AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW`を使用する。
CA sourceの非零cone chartは旗9変数と固定旗20係数からなるsmooth29次元空間である。

## 19個の滑らか変数を除くetale断面

固定旗はL=<∂y>, W=<h,x>, H=<x>である。
その連結parabolicのLie algebraは、全25行列成分から以下の9成分を除く16次元空間である。

    (0,2),(0,3),(0,4),(1,0),(1,2),(1,3),(1,4),(3,2),(4,2).

全16生成元について、固定旗CA20次元線形空間の全20基底への作用が同じ空間に入ることを
全150物理係数で検証する。このparabolicが連結なのでLie安定性は群安定性へ積分される。
P0でのfiber係数上のparabolic orbit微分は階数10である。
20係数の番号0,4,5,6,7,8,11,12,13,16への射影で非零10次minorを持つ。
従ってこれら10係数をP0での値に固定したaffine10断面と、parabolic内の10次元局所切断との
作用写像はP0でetaleである。

全normalizer sourceは、さらに旗の9局所変数を持つ。
ramification Fitting idealはG-equivariantなので、その完成環から旗9とparabolic10という
計19個の滑らか因子を取り出せる。
これはtargetへ任意の座標切断を行う操作ではなく、smoothなnormalizer source自身のetale分解である。

## 全critical idealを小Schur列で計算する

固定旗の全150×20係数行列Cは一定で単射である。
その単元20次blockでfiber20列を除いた、旗9列の全150×9行列Qを作る。
全source Jacobianの最大29次minorのidealは、Qの全9次minorのidealに等しい。

上記10断面の座標をs0,...,s9とする。
Qの行5,22,23,27,31,35,50,65、列0,1,2,4,5,6,7,8の8次blockは
この**全affine断面で行列式1**である。
従って残る列3をSchur消去した全150成分が、critical idealを正確に生成する。
多項式の一次部分だけではなく、全成分を使う。

## 全QQ idealと三変数表示

七つの座標は、全idealを法として次の多項式で消去できる。

    s0 = 3*s3/5 - 4*s7/5
    s5 = (12*s7² - 3*s3²)/5
    s1 = 7*s5/12
    s4 = 35*s5/24 - 2*s7²
    s6 = 0
    s8 = s7²
    s9 = s7³/3.

残る三変数を

    w=s2,  u=s3+2*s7,  v=7*s3-6*s7

と取ると、全残余idealは正確に

    (w*u, w*v, u*v, u³)

となる。逆変換はs3=(3u+v)/10、s7=(7u−v)/20である。
全Schur生成元と、上の七消去式＋四生成元との両ideal包含をQQ上で検証する。
さらに別のMacaulay2計算で二つのidealの等号を検査する。

したがって、全normalizer ramification schemeのcone完成環は

    C[[g1,...,g19,w,u,v]]/(wu,wv,uv,u³)

である。非零coneのscalar方向を除いたprojective版ではg変数は18個となる。

## schemeの帰結とscope

三変数idealには正確なprimary分解

    (wu,wv,uv,u³) = (u,v) ∩ (u,w) ∩ (w,v,u³)

がある。従って被約critical locusは二本の軸のnodeで、junctionにembedded primeが一つある。
nilradicalは(u)、その三乗は零、二乗は非零であり、nil部分はjunctionに長さ2で支持される。
三変数環の次元1、埋め込み次元3、depth0、重複度2なので、cone完成環は
次元20・埋め込み次元22・depth19、projective版は次元19・埋め込み次元21・depth18である。
この接空間次元22は、先のintrinsic ramification cotangent rank7=29−22と一致する。
先の接空間だけから実critical次元を推測したものではない。

これらはsourceのramification schemeの結果である。
CAのtarget完成環、rank2底上の全odd ideal、実Dのnilradical・associated primes・depthには
追加の有限像計算が必要であり、今回解決したと主張しない。

* producer: `code/ambient_B18_CA_critical_scheme_20261103_v63.py`
* evidence: `reports/AMBIENT_B18_CA_CRITICAL_SCHEME_20261103_v63/`

全150係数、全parabolic作用、etale minor、全Q、全Schur式、両ideal包含、七消去式を保存する。
