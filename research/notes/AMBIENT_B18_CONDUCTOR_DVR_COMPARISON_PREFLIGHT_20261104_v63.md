# Full D0 socle比較のDVR化：actual20-base、raw coker26と残る完全性gate

2026-11-04。K=ker(D0→E4)、C=coker(D0→E4)とする。
両方がI_A I_Bに消されることは既独立受理済み。
ここではDVR上の有限性と具体的なcokernel入力を証明し、kernelのTor補正を省略せずに残す。
大きいSmith、全640 second carries、新しいsource jetは実行しない。

## 1. ellは実際の両conductor曲線のparameter

ell=5*(s_CA-r_CA)は共通R0のactual CA F9函数である。
既存source2のphysical inverseを使うと、CB自身のzero7 planeへの一次制限は

```
ell_B=(35/3)r_B-(85/3)s_B+O(m²).
```

cached actual cubicでg_A3(1,1)=125/24、g_B3(17,7)=-181500。
従ってどちらのconductor CI曲線でもellはparameterかつ非零因子で、intersection multiplicityは3。
Weierstrass preparationにより両R0/I_iはD=QQ[[ell]]上free rank3。
CBでr_Bを横断変数に取ると、そのcubic leading unitは-181500/4913である。
nonlinearなCB座標変更はunitなので、このleading multiplicityを変えない。

## 2. Product curve ringのDVR module型

scheme-union R0/(I_A intersect I_B)は二つのfree rank3 moduleのsubmoduleなので、D上free rank6。
前packetの有限kernelはA3^6 direct_sum omega_A3、長さ21。
A3上でellの乗法はrank1であり、omegaでもrank1。
よって各三次元blockはD/(ell²) direct_sum D/(ell)となる。
free quotientへの拡張はD-moduleとして分裂するから

```
R0/(I_A I_B) = D^6 direct_sum (D/ell²)^7 direct_sum (D/ell)^7
```

という非canonicalなmodule同型がある。special fiberの次元は20。
これはringのmultiplication tableやcanonical20 basisを既に計算したという主張ではない。
また、その任意のmoduleのtorsion長を21で制限するものでもない。

## 3. Normalizer raw cokernelはさらに具体的

physical quotient前のraw mapは

```
A_A tensor_R0 A_B -> B_A tensor_R0 B_B=R0^4.
```

像はR0 direct_sum I_A*tauA direct_sum I_B*tauB direct_sum I_A I_B*tauA*tauB。
従ってraw cokernel C_rawは

```
R0/I_A*tauA direct_sum R0/I_B*tauB direct_sum R0/(I_A I_B)*tauA*tauB,
C_raw = D^12 direct_sum (D/ell²)^7 direct_sum (D/ell)^7.
```

D-generator数は26、intrinsic torsion relationsは14。
raw Tor2 kernelのD型は(D/ell²)^21 direct_sum (D/ell)^21で、length63、D-generator数42である。
これはその像がactual D0で生きるという主張ではない。

## 4. Actual cokerはright-exactに計算できる

十二physical差をF_j=delta_j+b_Aj*tauA-b_Bj*tauBと書く。
F_jはraw image algebraの元なので、その乗法はC_raw上のendomorphismになる。
三つのR0-generator(tauA,tauB,tauA*tauB)に関するmatrixは正確に

```
[[delta, 0, -nuB*bB],
 [0, delta,  nuA*bA],
 [-bB, bA, delta]].
```

各行はそれぞれI_A,I_B,I_A I_Bを法とする。
source annihilatorを保ち、異なるF_iのcommutatorは各行のideal内に入るため、これらは可換endomorphismである。
全coefficientはpure CA48/CB48の同R0関数であり、旧r/g微分を計算する必要はない。

右完全性からactual CはC_raw/(sum F_j*C_raw)。
actual26 DVR-generatorとR0-actionをWeierstrass sourceからbindすれば、
14 intrinsic relationsと十二actionの26列ずつ、計**26 generator/326 relation以下**のDVR表示となる。
この小表示は具体的に構成可能な次の入力だが、本preflightではactual DVR basis/actionsをまだ展開しない。

## 5. Kernelで省略できないTor gate

K,CがD上有限であることは、source matrixをell=0で割ってkernelを取ってよいことを意味しない。
例えばD --ell→ Dはkernel0だが、special fiber mapのkernelはQQ。
逆にD/ell²→D/ellのkernelはQQだが、二つのspecial fiber間のmapは同型である。

現表示P:R0^640→R0^9、Q:R0^48→R0^4、U:R0^9→R0^4に対して、UP=QBのactual relation liftが必要である。
これは28 intrinsic列では零、他612列では「68本のCB polynomial times9 C0 basis」に因数分解できる。
CB zero kernelは28 syzygy+36 products+4 graphの全68本で、target次数≤2。
そのHadamard差F(qA)-F(qB)から、48-column relationへのliftを同じR0上で作れる。
旧scalar regular cofactorを新たな0/0へ置換しない。

しかしこれだけでも足りない。L=ker[U,-Q]と書けば

```
K=L/( columns(P,B) + columns(0,Syz(Q)) ).
```

**Syz(Q)の完全性**がまだ未入力である。
E4のall-order factor/core4 reductionにより、それは五Koszul方向とactual七core函数のsyzygyを使う小さい問題へ移せるが、
現在の二次initialだけからその全syzygyを仮定してはならない。
complete K presentationが得られた後なら、I_A I_Bでtensorしてもright exactnessにより同じKとなり、20-generator DVR baseへ安全に移せる。

## 6. 有限complexが完成した際の終了criterion

H1=K,H0=C、他のhomology0を持つactual comparison complexをD上の有限free complexへ移せれば、
Dのglobal dimension1から

```
H2(derived special fiber)=K[ell],
0 -> K/ellK -> H1(derived special fiber) -> C[ell] -> 0.
```

従ってK torsionの有無はH2の有限QQ rank検査で決定できる。
またmu_D(K)>rank_D(K)と同値であり、generic rankとspecial fiberを比較する場合も上のC[ell]補正が必要。
ここでいうcomplexは完全な比較解像度であり、二つのpresentationをspecializeしただけのmatrixではない。
このcomplete complexとactual actionsは未構成なので、本packetはKの有限socleやD0深さを決めない。

Code `code/ambient_B18_conductor_DVR_comparison_preflight_20261104_v63.py`。
Report/folder `AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63`。
