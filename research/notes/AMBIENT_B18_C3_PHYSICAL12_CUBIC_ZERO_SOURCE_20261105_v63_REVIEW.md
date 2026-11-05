# actual C の HF366 と physical12・tauAB三次 subset の独立 SOURCE 証明

2026-11-05。PRIMARY report
`AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63/report.json`
SHA256 `f2f0cfbc738cd58167741b9ce6cd2a73cec5b884dfe0f47fbb1acf94df19d72c`
に対する独立監査である。REVIEW report は
`reports/AMBIENT_B18_C3_PHYSICAL12_CUBIC_ZERO_SOURCE_20261105_v63_REVIEW/report.json`
SHA256 `276b80bf7ea5edd5fd0abf23e33021a6dcb68bc294d64d1c47a15d300f432b75`。
author producer の import、三次 scalar 値、source/gauge/jet solve、Gaussianは使っていない。

## 完全な source と同じ低次係数

対象は同じ (R_0=\mathbb Q[[z_0,\ldots,z_6,r,\ell]]) 上の比較加群
(C=\operatorname{coker}(D_0\to E_4)) であり、Cprimeや環の積ではない。
既受理の raw cokernel は三つの target (e_A=\tau_A,e_B=\tau_B,e_H=\tau_A\tau_B) を持つ。
source image の ideal 記述と右完全性から、全表示は

\[
I_Ae_A\ (8),\quad I_Be_B\ (8),\quad I_AI_Be_H\ (64),\quad
F_j(e_A,e_B,e_H)\ (12\cdot3)
\]

の **116列** である。Fj*1はすでにraw imageに属する。実際の作用行列は

\[
\begin{pmatrix}
\delta&0&-\nu_Bb_B\\0&\delta&\nu_Ab_A\\-b_B&b_A&\delta
\end{pmatrix}.
\]

全12行の C2 の delta,bA,bB は、保存された元164 root の値から直接再構成して一致を確認した。
新しい pure-ell 二次値も同じ元root番号で接続した。
全delta,bA,bBの定数、両nuの定数、deltaの全9一次係数が0であることを明示的に確認した。
nu,bA,bBのplane一次係数は0、deltaとbAのplane二次係数は三方向とも0である。

C2の74 normal formが普通の次数を下げないことも確認した。
従って高次tailから定数・一次・plane二次成分への混入はない。
monic z積の二次初項はzzだけで、r²・rℓ座標のfunctionalsは実際の低次plane成分を読んでいる。

## 一次21と完全な二次 source

accepted all-order graph frame (I_A=(z,g_A),I_B=(z+\beta,g_B)) を用いる。
gA,gBはorder3、betaはorder2以上である。

\[
R_{A,i}=z_ie_A,\quad R_{B,i}=(z_i+\beta_i)e_B,\quad
R_{H,i}=F_i e_B\quad(i=0,\ldots,6)
\]

の一次部分は、それぞれ zi eA,zi eB,zi eH で、literal I21 を作る。
RHの二次plane部分は実際のdelta2,bA2の零性から0である。
RBの二次plane beta2 eB は残す。

独立 checker は116列の全一次写像を組み立て、選択21列のI21を確認した。
非選択列に対して「その列−選択21列×同じ一次係数」を取り、free-tail I95を持つ
**完全な95個の constant kernel**を構成し、全一次multiplybackが0となった。
内訳は二curve、64product、12tauA、残り5tauB、12tauABである。

M0=0なので、任意の二次 initial relation は、constant first-kernelの二次像と
線形係数×M1から生じる。後者はzを法として消える。
95列のうちplaneに残るものは12tauA列の (-b_{B,j,2}e_H) のみである。
その一次消去はRHしか使わず、RH2plane=0なのでcarryは追加されない。
二curveはorder3、64productの二次項はzz、残りtauBはdelta2/bA2plane=0、
tauABのoff-diagonal二次項はnu1*b1のzz、H成分はdelta2plane=0である。

bB2の行18,60,67と列r²,rℓ,ℓ²の小行列式を独立な六permutation和で確認し、
値は5/768となった。投影二次像はその負符号なので、その三列minorは−5/768である。
いずれも非零で、二次spanはplane² eHの全体である。
よって

\[
\operatorname{HF}_{0,1,2}(C)=(3,\ 3\cdot9-21,\ 3\cdot3-3)=(3,6,6).
\]

ℓ³Cはm³Cに含まれるため、C3=C/ℓ³Cも同じ次数0〜2を持つ。

## 指定12列の coherent な三次 clearing

Mj=Fj eHはM1=0である。M2のA,B成分はzz、H成分delta2はzを含む。
全98種類の quadratic monomial/targetについて、係数を未知のままにして
unit一次frameへの clearing を独立に構成した。

- Aのzzはz倍のRAで消し、RA2=0なのでplane三次carryはない。
- Bのzzはz倍のRBで消す。RB2のbeta2 eBを含む全二次成分にzが掛かるため、
  plane三次への射影は0である。
- Hのzzとzi*r,zi*ℓは、それぞれz,r,ℓ倍のRHで消す。
  RH2のplane成分は実際に0であり、nonplane部分はzを保つ。

これで未知のzz係数やmixed係数を捨てずに、

\[
\pi_{plane}(M_{3,j}-M_{2,frame}L(M_{2,j}))
\]

を扱える。off-diagonalのnu*bは、両因子の定数とplane一次が0なので、
plane上ではorder4以上である。残るのは **未知のdelta3 eH**だけである。
その四plane単項式r³,r²ℓ,rℓ²,ℓ³はすべて既知のplane² eHの一次倍なので消える。
12×4の未知delta3 sensitivitiesを形式的に残したまま、全ての投影が0となる。

従って、指定した12tauAB列の二次をcoherentに消した後の三次クラスは全て0である。
元の列自体が零、primitive V3自体が零、全C3三次関係が零という主張ではない。

参考の plain relation model は二つのCI(3,3)とplane/m² eHを合わせたもので、
HF(3,6,6,4,2)、長さ21となる。これは実際のC3の長さ・torsion指数を決めない。
元H0やnilcubeにも新しい結論を与えない。
