# 実cubicのheavy二係数単位とDelta34チャート：独立監査

2026-11-04。係数環は R=QQ[a,e,d,1/(adFBCJ)]。
authorの係数抽出・solverを再実行せず、Fraction係数の疎辞書で再構成した。
新しい点、素数、minor、pivot、GB、四次係数は使わない。

## 実sourceとの結合

受理済みrow-basisのhashに一致する凍結57×91を直接読む。行42,47はnormalized43,49、
元135式のcubic20,61である。二行の正規化単位はともに1/625000000。
91単項式からHを持つ七列を独立に抽出し、全14係数から差を作る。
最初の41 product行は七H係数が全て零であることも実行列で確認する。

物理変数と旧座標の実際の代入

`(v0,v1,v2,v3,v4,v6,v7)=(x0,x1,x2,W,x3,Z-3x0+x1,x4)`

から係数変換行列を転置で独立生成し、保存された7×7行列と逆行列に照合する。
差の最初の五係数は厳密に零で、残るZ,W係数は

`kZ=-aFB*g1*g2/10^7, kW=-FB*(Cg1-B)/10^7`

となる。全七係数を原sourceから確認しており、多項式Bézout式だけの監査ではない。

## 単位と同じcokernel

g1=(3e-2)d+20、g2=(6e+1)d+40、h=Cg1-B、ell=6e+1、
K=600a+126e²-63e-14-ell*C と置く。実係数から独立に得た
cZ=KC、cW=a(2ell B-Kg2)、delta=-64000C²について

`cZ*delta*kZ+cW*delta*kW=192*a*F*B²*C³/125=eta`

をQQの全係数で検証する。etaとdeltaはRの既知単元である。
g1,g2,hや未知minorを反転していない。
T2と保存逆の両方向の積の全8項、det(T2)=1、係数対が(1,0)へ移ることを
分母払いした正確な多項式演算で確認する。row係数行列へは右からT2の転置を掛ける。

行47をdelta*(行47-行42)へ置き換える変換の行列式はdeltaで単位。
純84列は保持し、Hの七列だけに上記可逆変換を行う。得た係数1のpivotを消去すると
元の実57→91と同じcokernelを持つ自由56→90となる。
product41行はH係数を持たないのでこの操作で変化しない。
56×90の全係数を展開した主張ではなく、全域で定義された有限可逆recipeである。

## 階数57は追加Delta34開集合上

受理済みの全点product41とZW product19から、全R上でvector bundleの全射
B50→V21があり、kernelはrank29。compact source relationが、projected genuine16の
primitive核ベクトルG61-G20を全域で与える。
既に受理済みの実restricted36×40の固定34minorをDelta34とする。
Delta34非零上では、product部分rank19とminorのlower34からprojected genuineのrank≥15。
既存の核でrank≤15なのでちょうど15で、核はこのprimitive lineである。

上のheavy functionalはproduct41を消し、delta*(G61-G20)を1へ送る。
従ってfull genuine16にfiber上の関係があればprojected核lineへ入り、その係数も零となる。
よって**R∩D(Delta34)の全点で、実full cubic rankは57**。
Delta34の元の保存値が非零という既存証明は再利用し、新たに評価しない。
Delta34をRの全域単元とはみなさない。

全FULLUでのrank57、四次rank239、U78、CI、全kernel exhaustionはこの監査に含まれない。
既知uniform kernel38と残り42も変えない。
