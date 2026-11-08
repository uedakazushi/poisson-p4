# constant modular：被約nonconcurrent最高三次式の全known被覆

2026-11-08、PRIMARY `E_multiple_hyperplane`、担当 `/root/plane_collision_review`。
**全源の純係数・幾何証明。新CAS/operator/profile/点/minor/jet/GBは実行しない。**
A=0,c≠0で最高三次式f3が被約、かつ三concurrent linesではない全源を扱う。
三concurrent lines、非被約f3、A=c=0のunimodular源は本稿の残余として区別する。
producerの本稿の独立受理は別途必要である。

## 1. 必要十分入口と使う最弱補題

元の全源はactual affine frameと非零tensor尺度により

```
P=∂t∧Y+Jac(f), divY=1, Yf=f, degY≤2, degf=3
```

と書ける。全Y2/Y1/Y0とf3/f2/f1/f0を保つ。
最高条件はY2(f3)=0、divY2=0である。
この二条件の全解がJac(f3)に対するlinear Hamiltonianなら、actual shear
T=t−ℓ(x,y,z)で全fのHamiltonianを引いてY2を消せる。
逆はt=T+ℓ、det1。全divergenceとYf=fは保たれる。
以下の全syzygyはlinear Hamiltonianの厳密な三次元空間に尽くされる。

## 2. 全被約nonconcurrent plane cubicと全syzygy

smooth irreducible cubicでは三偏微分が正則列なので、degree二syzygyは
三Koszul generatorだけである。divergenceはmixed partialの相殺により零。
既約node/cuspの全normalformsと全低次係数は別稿
`CONSTANT_MODULAR_IRREDUCIBLE_SINGULAR_CUBIC_WHOLE_MPS.md` §§2–4で手比較した。
nodeは全syzygy四係数でdiv=−(c+k)z、cuspは全五係数でdiv=−py−3qzなので、
div零は両者とも三つのlinear Hamiltonianそのものである。
その補題はsmoothのKoszul完全性をsingularへ仮定したものではない。

残る全被約型はline+smooth conicの横断／接配置、及び三nonconcurrent linesである。
それぞれ可逆GL3で以下のFへ移る。全体係数c≠0は保持し、
f3=cFなら最後のHamiltonian関数を1/c倍するだけでよい。

### 三nonconcurrent lines：F=xyz

各prime factorの保存から(A,B,C)=(xL1,yL2,zL3)、L1+L2+L3=0。
div零はL1のx係数、L2のy係数、L3のz係数が零となることと同値である。
従って全解は

```
A=bxy+cxz, B=dxy−cyz, C=−dxz−byz
  = components of gradF cross grad(d*x−b*y+c*z).
```

全六係数syzygyへdivの三条件を課した全解であり、非零係数を反転しない。

### 横断line+conic：F=z(xy−z²)

prime factor zを保存するのでC=zL、他方のfactorから
 yA+xB=−Lxy+3Lz²。
mod(x,y)よりL=a x+b y、二変数Koszul syzygyから全解は

```
A=xM−bxy+3bz²,
B=−yM−axy+3az²,
C=z(ax+by), M=cx+dy+ez.
```

divはcx−dy。従ってc=d=0が必要十分で、全div零解は
 gradF cross grad(−a*x+b*y+e*z)。
全五係数syzygyの必要十分表示である。

### 接line+conic：F=x(xz−y²)

A=xL、factor q=xz−y²の条件から、全解は

```
A=xL, B=−yL/2+xM, C=−2zL+2yM,
L=ax+by+cz, M=dx+ey+fz.
```

divは(e−a/2)x+(2f−2b)y−7cz/2。
従ってe=a/2、f=b、c=0が必要十分であり、残る全解は
 gradF cross grad(d*x−a*y−b*z/2)。
全六係数syzygyを最後まで保っている。

以上でsmooth、node、cusp、line-conic横断／接、triangleの全型を処理した。
被約plane cubicの他の型は三concurrent linesだけである。
三線が同一点を通ることはその三linear formが二次元spanに入ることと同値なので、
本稿では「essential三変数の被約f3」と言い換えても同じ範囲となる。

## 3. 全affine Yを既知閉和へ送る補題

ここではYは**任意のaffine場**、divY=1、Yf=f、f3被約とする。
Y1のsemisimple/Jordan/zero-eigenvalueやf2/f1/f0の共鳴を先に分類しない。
次の実二次Poisson族を作る：

```
Pε=(1+εt)∂t∧Y+Jac(f).
```

Yはaffineなので全係数次数≤2である。
Yf=fとdivY=1から[Y,Jac(f)]=0、従って全εで元Jacobiが零。
signed Pf=(1+εt)f、actual D=h(h+εt)Fであり、Fはfの三次斉次化。
f3が被約なのでFも被約である。実際、Fのprime factorの平方があれば、
hを因子に持たないこととh=0への制限からf3も平方因子を持ってしまう。
Fの次数は三、hはFを割らず、ε≠0ではh+εtもFを割らない。

* F既約：actual Dは相異なる二Hと既約cubic。
  actual affine変更x=1+εtでPε=x∂x∧V+W、V=εY、W=Jac(f)。
  旧 `TWO_HYPERPLANES_CUBIC_SUPPORT.md` の行列Bは正確にεY1、
  traceB=εなので半単純部分は必ず非零である。
  同稿の**全nonzero semisimple B支持定理**がそのまま適用される。
  全低次c,dとW45を保持する定理なので、generic固有値を仮定しない。
  よってこの全源はMPS40の閉和へ入る。
* F=LQ、Q既約：三つの相異なるHと一つの既約quadricの全支持定理を使う。
  各Q rank、全配置、全frameを保った既known閉和への定理である。
* F=L1L2L3：五相異Hの全配置定理を使い、T01へ入る。

したがって全ε≠0は有限個の既知閉成分の和Kに入る。
全係数はε多項式で、ε=0もrank四なので、Kの閉性から元PもKに入る。
これはGrassmannian/sourceの一般平坦性からの推測でも、
一つのSNC極限から任意の未知fiberを被覆したとの推論でもない。
全εで明示されたactual Poisson曲線を使う。

## 4. 全元sourceと補集合

§2の全div-zero syzygyから実t-shearでY2を全消去し、§3を適用する。
その後の可逆shear／frame／tensor尺度を戻すことで、元の全低次項と全Y2を含む
被約nonconcurrent f3源の全点がknown有限閉和Kに入る。
Y1が零固有値を持つ場合、Y0を除けない場合、f2/f1/f0がresonanceに乗る場合も、
§3ではそれらを一切反転しないので同じ全sourceに含まれる。

本稿は三concurrent linesに同じHamiltonian消去を適用しない。
f3=xy(x+y)ではdiv-zero vertical Q2(x,y)∂zが全三係数あり、
linear Hamiltonian vertical像はgradf3の二次元spanなので一方向余剰がある。
非被約f3では更にsyzygyが増え、本稿の外に保持する。
またA=c=0のunimodular全h⁵源も別義務である。
新成分増分零、全A0の完成・全分類の完成は主張しない。
