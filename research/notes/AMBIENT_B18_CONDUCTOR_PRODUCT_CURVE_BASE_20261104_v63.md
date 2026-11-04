# 二つのactual conductor CI曲線のproduct base：Torと有限embedded部分

2026-11-04。I_A,I_Bは同じactual R0九変数環内のzero odd idealである。
各idealは七つのunit一次座標と、それを零にした非零三次函数で生成されるので、各々CI8。
既受理combined105/Nakayamaの正確な等式は

```
I_A+I_B=(z1,...,z7,r²,rs,s²),  A3=R0/(I_A+I_B)=QQ[[r,s]]/(r,s)².
```

二族の座標をabstract型だけから合わせたものではない。

## 1. Combined Koszul complexを完全に小さくする

二つのCI8 Koszul resolutionのtensorは、combined16生成元のKoszul complexである。
そのhomologyがTor_i^R0(R0/I_A,R0/I_B)となる。
和のidealのminimal generator数は10。local環のunit column操作で、16生成元を
このidealの任意のminimal10生成元と六つの零へ変換できる。
実際、unit minorで10本を取り、残り六本をそれらの線形結合として引けばよい。
二つのminimal generator組の間のmatrixもresidue fieldで可逆なので、
さらに上の(z7,r²,rs,s²,0^6)へ移せる。

七regular座標のKoszulを除けば、残るcomplexは

```
K(r²,rs,s²; QQ[[r,s]]) tensor exterior(6 zero generators).
```

ここでunit操作matrixの巨大展開は必要なく、その存在と正確なideal等号が根拠である。

## 2. 三quadric部分の全homology

K3の一次cycleを

```
u=(-s,r,0), v=(0,-s,r)
```

と取る。二次boundary三本は正確にr*u、s*u+r*v、s*vである。
従ってH1は、二generatorにこの三relationを持つmoduleである。
uをdual basis eps_sへ、vを-eps_rへ送ると、これは正確にomega_A3となる。
長さ3、minimal generator数2、socle1。
H0=A3であり、grade2・三生成元Koszulの上限からH_i=0(i>=2)。
後者は三つのKoszul matrixのkernelを直接見ることでも従う。

したがってR0-moduleとして非canonicalな同型

```
Tor_i(R0/I_A,R0/I_B)
 = A3^(choose(6,i)) direct_sum omega_A3^(choose(6,i-1))
```

を得る。長さはi=0,...,8の順に

```
3,21,63,105,105,63,21,3,0.
```

とくに

```
(I_A intersect I_B)/(I_A I_B)=Tor_1=A3^6 direct_sum omega_A3
```

は長さ21、socle dimension13であり、I_A+I_Bにannihilateされる。

## 3. Product curve ringの正確な有限部分

完全列

```
0 -> R0/(I_A intersect I_B) -> R0/I_A direct_sum R0/I_B -> A3 -> 0
```

から、二つのCI curveのscheme-theoretic unionはdepth1。
各conductor idealがradicalであるとは仮定していないし、このunionを被約とは呼ばない。
さらに

```
0 -> (I_A intersect I_B)/(I_A I_B)
  -> R0/(I_A I_B) -> R0/(I_A intersect I_B) -> 0
```

を使えば、product curve ringのH0_mは上の長さ21moduleそのものとなる。

## 4. Actual normalizer比較への適用の限界

前packetのactual map D0→E4のkernel/cokernelはI_A I_Bに消されるので、
この一次元product base上で扱える。ただし**任意のそのmoduleのsocleが長さ21以下**とは従わない。
curveのparameterの任意高冪で割ったmoduleを考えれば、有限長は任意に増える。
よってこのring自体のH0と比較kernelのH0を同一視しない。
また、この積がT0全体をannihilateすると主張しない。

この結果は、actual比較kernelのfull presentationを計算するときの小さいbaseを指定するものである。
E4のsocle・比較kernelの非零性・全Dのdepth19/20は未決。
新しいsource jet、point、prime、GB、Tor comparison matrixの大きい計算は行わない。

Code `code/ambient_B18_conductor_product_curve_base_20261104_v63.py`。
Report `reports/AMBIENT_B18_CONDUCTOR_PRODUCT_CURVE_BASE_20261104_v63/report.json`。
