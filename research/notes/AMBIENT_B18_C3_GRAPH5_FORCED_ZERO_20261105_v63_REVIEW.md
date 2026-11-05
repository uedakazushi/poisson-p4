# 五 graph 三次係数の source と次元による強制消滅：独立監査

2026-11-05。PRIMARY `AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63/report.json`
SHA256 `a35b929755230acc1621df59d25d21a55d2617892925a1b2a027d7b5b1451d8b`
に対する独立 REVIEW は
`reports/AMBIENT_B18_C3_GRAPH5_FORCED_ZERO_20261105_v63_REVIEW/report.json`
SHA256 `04e85dc1ec29a4387ce77958932cb9b42cdef2f676d4c41cadaf2e09b3a2f92c`。
graph rootのreadback、陰関数・gauge・jet solve、rank search、author importはない。

## 同じ実際の C と長さ下界

対象はR0-module (C=\operatorname{coker}(D_0\to E_4))、(C_3=C/\ell^3C) である。
accepted generic DVR peer の実rankは4。accepted C2 の次元14とℓ作用rank6から

\[
\dim(C/\ell C)=14-6=8
\]

が分かる。D=ℚ[[ℓ]]上の有限加群のtorsion summandは4個。
C2の長さから指数1は2個、残り2個は2以上となり、

\[
\operatorname{length}(C_3)\ge3\cdot4+1+1+2+2=18.
\]

有限baseのhomologyから推測した下界ではなく、同じ実Cの受理済みDVR rankと
right-exact C2 sourceの値からのPID計算である。

## 既知の graded cover

first21と完全なquadratic sourceの受理済み証明から、gr(C3)は

\[
M=\frac{\mathbb Q[r,\ell]}{(g_{A,3},\ell^3)}e_A
\oplus\frac{\mathbb Q[r,\ell]}{(g_{B,3},\ell^3)}e_B
\oplus\frac{\mathbb Q[r,\ell]}{(r,\ell)^2}e_H
\]

の商である。これは実際の加群のdirect-sum分解ではない。
元sourceの二つの三次式を別実装で合成し、r³係数125/24、−15125/144を確認した。
各枝でr³を逐次置換するmonomial rewritingから、全NF3の4×12、NF4の2×15を再構成した。
単にrelationsをannihilateするだけでなく、degree3の関係minorがg_r³、degree4ではその二乗という
unit triangular blockを持つことから、quotientの次元2・1とfree座標を確認した。
従ってMのHFは(3,6,6,4,2)、長さ21である。

## 同じ source に現れる d

各行 j=63,64,68,73,74 では、bB2のplane成分が0なので、tauAの完全kernelの同じfree unit列
e_jが使用できる。受理済みtauA sourceのcoherent clearingはeA成分を変えず、
三次eA成分は実際のdelta_j3である。eBの全carryは未知のまま残す。

tauB列は一次が0、delta2のnonzz部分は0である。
ただし **bA2のH成分には実際のzi*r、zi*ℓ項がある**。
これをr/ℓ倍のRHで消す。RH2のeB-plane成分は既に確認したdelta_i2-plane=0なので、
このmixed clearingはeB三次planeを変えない。zz部分はz倍のframeで消し、
beta2などのcarryもzを保つ。eH三次は既知のplane²eHにより消える。
従ってtauBの三次はeA=0、eB=delta_j3となる。

元C2のliteral Taylor dualからdelta_j3のr³、r²ℓ成分は0である。
ℓ³を法としてdelta_j3=d*rℓ²と書ける。未知のℓ³成分は0と置かない。
M3の座標を(A_r²ℓ,A_rℓ²,B_r²ℓ,B_rℓ²)とすると二つの実際のsource関係は

\[
p=(0,d,u,v),\qquad q=(0,0,0,d)
\]

となる。u,vは全てのopaqueなB carryを含む。固定2×2minorはd²。
rを掛けた実際の関係のM4座標は

\[
rp=(d,-(g_{B,r^2\ell}/g_{B,r^3})u+v),\qquad rq=(0,d)
\]

であり、こちらの固定minorもd²である。
両恒等式を新しい自由多項式環ℚ[d,u,v]で検査した。

d≠0なら、Mの次数3で2次元、次数4でさらに2次元が実際のinitial relationsにより消える。
従って実C3の長さは21−2−2=17以下となり、上の下界18に矛盾する。
よって **五行すべての実際のdelta_(r*ℓ²)は0** であり、五tauB graph cubic classesも0となる。
非零仮定の場合のquartic witnessはr倍であり、ℓ倍ではない。

全三次sourceの消滅、実C3の長さや全torsion指数、元H0長さ、nil指数は決定していない。

## 実行しなかった数値経路の保存

最初の graph5 prepared worker とsource majorantには rational c-nodeの分母を落とす不一致があり、
またfull列のmixed H項を省いたproof wordingがあった。
独立 precheck `a56c4a4d…` に32個の非整数定数と実mixed項を保存し、値の実行前に停止した。
別V2はauthor側で用意されたが、今回の純粋な証明ではV1/V2どちらのreadbackも必要ない。
両方とも未実行のまま保存し、この監査はV2数値workerのSOURCE PASSを代用しない。
