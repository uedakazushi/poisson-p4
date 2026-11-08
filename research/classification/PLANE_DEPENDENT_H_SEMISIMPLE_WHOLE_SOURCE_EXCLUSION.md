# rank零plane・dependent H二源は既約quinticで空

2026-11-08、PRIMARY E、純全源証明。元全K1と全J2の完全性は既独立六源receiptを使う。以下はG像20/21という次元上界からの帰属推測ではない。新しいM1、profile、点、normal jet又はGBを計算しない。別作者読解は未実施。

## 1. 既存全kernelからの実六shear

元座標は(a,b,c,z,w)、H=a∂a+b∂b−2c∂cである。対象は既六源表の

```
SS_H_RESONANCE:  Π0=H∧(z∂z−w∂w)−5zw∂z∧∂w,
NIL_H_RESONANCE: Π0=H∧z∂w+(5/2)z²∂z∧∂w.
```

先行別作者receipt `PLANE_DEPENDENT_SIX_SOURCE_SEPARATE_REVIEW.md` は、元全350×30/9、全kernel、実六shear及び全非零minorを確認した。本稿ではそれを再計算せず、保存された完全kernelの全係数を使う。いずれもJ2は350×9のrank九であり、固定P1に対してweight二Jacobiを満たすP2は存在すれば一意である。

K1座標αを元reportのcomplete_kernelの列順に取る。実六shearをa_i→a_i+s_i z+t_i wとする。SSの全八座標に対し、shearの像は

```
(α1,α2,α3,α4,α5,α6)=(3t_a,3t_b,9t_c,−3s_a,−3s_b,−9s_c),
```

で、残るα0,α7は任意。nilの全十座標に対しては

```
(α0,α1,α2,α6,α7,α8)=((3/2)s_a,(3/2)s_b,(9/2)s_c,
                         (3/2)t_a,(3/2)t_b,(9/2)t_c),
```

で、残るα3,α4,α5,α9は任意である。これは保存full30×6のfree-coordinate行そのものであり、定数逆が全域に存在する。六次元の数え上げだけではない。

実finite shearはnormal weightを上げ、Π0はweight0/1/2、P1はweight1/2、P2はweight2のままになる。canonical P1へのEuler regauge後も、全九P2への作用は元P2にshear/residualの多項式を加えたtriangular写像で、九変数の線形部はidentityである。従って逆shearにより全sourceを以下の残sourceへ移せる。P2は零と仮定せず、全九方向を残す。各定数逆とfinite inverseは既先行全source shear argumentと同じ実変換である。

以下、三変数(c,z,w)のJacobian bivectorを

```
Jac(F)=F_w ∂c∧∂z−F_z ∂c∧∂w+F_c ∂z∧∂w
```

と書く。

## 2. SS H全源はrank二以下

E∧(z∂z−w∂w)を引くとΠ0は−3Jac(czw)になる。残る全P1は

```
α0 z²∂c∧∂w+α7 w²∂c∧∂z
   =Jac(−α0 z³/3+α7 w³/3).
```

従ってP2を加える前に、全residualは既に三変数c,z,wの厳密Poisson bivector

```
Jac(−3czw−α0 z³/3+α7 w³/3)
```

である。P2=0は元全projective Jacobiの解であり、J2のinjectivityから任意の実解のP2はこれに等しい。全tensorは射影rank二以下なので、rank四／非零actual quinticの枝には存在しない。六shearとEulerの逆を戻してもrankは変わらない。

## 3. nil Hの全sourceと明示Euler regauge

残るkernel列をp=α3、q=α4、r=α5、s=α9とする。元P1の全式は

```
p[−bw ∂a∧∂b−5cw ∂a∧∂c+zw ∂a∧∂z]
+q[ aw ∂a∧∂b−5cw ∂b∧∂c+zw ∂b∧∂z]
+r[−aw ∂a∧∂c/3−bw ∂b∧∂c/3+zw ∂c∧∂z]
+s w²∂c∧∂z.
```

E∧[z∂w+p w∂a+q w∂b−(r/3)w∂c]を引くと、全必要源は正確に

```
Π=Jac(F)+(p∂a+q∂b)∧W+P2,
F=(3/2)cz²+(r/3)zw²+(s/3)w³,
W=w(−4c∂c+2z∂z+w∂w),
P2=A(z,w)∂a∧∂b+B(z,w)∂a∧∂c+C(z,w)∂b∧∂c,
```

となる。A,B,Cはそれぞれ全binary quadraticであり、九係数を一つも削っていない。全係数はa,bに依存しない。

## 4. 一意P2とactual stabilizerがlinear Casimirを強制する

GL2(a,b)はΠ0、Jac(F)、Wを保ち、p∂a+q∂bだけを標準二次元表現で動かす。全九P2空間もこの実作用で保たれる。Euler regaugeは元projective Jacobiを保ち、J2は同じinjective mapのままである。

(p,q)=(0,0)ならΠ0+P1はGL2全体に固定される。実Poisson解P2が存在すれば、その任意のGL2像も同じweight二方程式の解である。一意性からP2は全GL2不変である。しかし九P2空間は

```
[Λ² C² ⊕ (C²⊗C∂c)] ⊗ Sym²〈z,w〉
```

であり、不変vectorを持たない。従ってP2=0、全Πは三変数Jacobianでrank二以下である。

(p,q)≠(0,0)では実GL2でp∂a+q∂bを∂aへ送る。例えばp≠0では

```
M_ab=[[1/p,0],[−q/p,1]],     inverse=[[p,0],[q,1]],
```

p=0,q≠0では[[0,1/q],[1,0]]（inverse [[0,1],[q,0]]）を使える。これで二chartと零補集合を全て保持する。

このnormal formのΠ0+P1は、∂aを固定する全実行列

```
[[1,t],[0,d]], t任意、d≠0
```

に固定される。P2の一意性からP2自身もこのstabilizerに固定される。∂a∧∂b方向はdでscaleするのでA=0、∂b∧∂c方向もd及びtで動くのでC=0である。残るのは任意のB(z,w)∂a∧∂cだけとなる。

従って全実解Πはb-rowを持たず、bはこのraw liftの線形Casimirである。Euler-Pfaffianは

```
D=b Pf(Π restricted to(a,c,z,w))
```

（座標の符号を除く）となり、零又はproperな一次因子を持つ。座標逆変換後の因子も非零一次式であり、被約既約quinticにはならない。

この論証はP1の全parameterでweight二compatibleと仮定していない。実解がある点だけで一意性を使い、その点のactual stabilizerからP2を制限した。chosen left inverseの大域equivariance、reductive averaging、schemeのradical計算、又は未証明のfull Poisson normal formを要しない。

## 5. 射程

SS H全源はrank二以下、nil H全源はrank二以下又はactual Dに一次因子を持つ。従ってどちらもrank四・被約既約quintic枝では空である。先行dependent N三源のconstant-factor部分は既essential三八成分、他八modulesは既proper因子／不整合なので、今回のH二源と合わせたdependent globalV全枝のE被覆へ使える。先行whole六源の必要還元と全multiplicityを維持する。

他のrank零plane symbol、scalar、independent globalV、co-Higgs k一、Π0零の枝は本稿だけでは閉じない。旧G像20/21証明及び旧report bytesを変更しない。
