# Twisted cubic-cone rank零の全源：T38への実認識と余枝の一次因子

2026-11-08、PRIMARY `E_reduced_nonnormal`。全保存 source／全 vertex kernelからの純証明。
**指定 twisted cubic-cone rank零・被約既約actualquintic の全点は既知 T38 に入る。**
新 matrix、M1、rank、minor、native、GB、point又はjet計算は用いない。
以下は別作者のwhole-source／全literal括弧読解の前の producer 証明であり、統合担当の受理を待つ。

## 1. 全十八源、全base分岐と実逆

accepted `CUBIC_CONE_ABELIAN_CONORMAL_NONDOMINANCE.md` と同 independent reviewを入力とする。
座標 `(x0,x1,x2,x3,c)`、idealは

`q0=x0*x2−x1², q1=x0*x3−x1*x2, q2=x1*x3−x2²`。

全三十rank零 sourceの abelian-conormal必要条件は完全十八kernelを与える。
baseの全六座標 `(A,B,C,D,E,F)` は元十八の indices `[1,6,10,12,13,15]`、
残十二は vertex block `W_i=P_{i,c}=Σ_g q_g*t_(4g+i)` そのもの。
ここでは W を raw upper-matrix entry と定義する。
旧noteの `∂c∧V` 表示では V=−Wであり、符号を混同しない。

全350 Jacobiは strict `[Pbase,Pbase]=0` と全350×12 vertex linear systemに分かれ、V²条件は無い。
全base reduced支持は

`L3={B=0,F=2A,D=−C}`、free(A,C,E)、
`Q2={C=0,E=6B,F=3A,AD=4B²}`

の和である。
base零は全tensor `∂c∧V` がdecomposableでrank≤2、対象rank四から外れる。
正しい actual SL2(cubic representation)×nonzero scalar全軌道定理 ORBITS_V2を用いる。
初回誤ラベルのORBIT reportを全被覆証明として使わない。

L3非零は discriminant E²+4AC が非零／零の二軌道、
Q2非零は一つの平方binary quadratic軌道。
各 actual transformation は元全十八sourceとvertex systemを可逆に写す。
以下は canonical pointだけの標本でなく、それぞれ全orbit／全linear fiberへ適用する。

## 2. L3 semisimple全六vertex sourceを四つの実shearで戻す

代表 `(A,B,C,D,E,F)=(0,0,0,0,1,0)` の base full tensor は

`P02=q0, P03=2q1, P13=q2`、他のbase entryは零。

保存済み完全6kernelへ全係数 k0,…,k5 を戻すと、raw vertex 全 entryは

```
W0=k0*q0+2k4*q1+k3*q2,
W1=k4*q2,
W2=(k2/2)*q0,
W3=k1*q0+k2*q1+k5*q2.
```

これは whole twelve-source affine kernelの全六方向であり、chosen magnetic particularではない。
actual PGL shear `c'=c+ell0*x0+ell1*x1+ell2*x2+ell3*x3` の下で

`W_i -> W_i+Σ_j P_ij*ell_j`。

従って全 parameter 上で

```
k0 -> k0+ell2, k4 -> k4+ell3,
k2 -> k2−2ell0, k5 -> k5−ell1,
k1,k3 は不変。
```

ell2=−k0,ell3=−k4,ell0=k2/2,ell1=k5で四係数を全て消せる。
inverse は反対のlinear shear、determinant一、唯一の分母二はconstant unit。
parameterを反転せず、k1/k3零の全補集合も保持する。
全sourceは正確に

`W=(k3*q2,0,0,k1*q0)`

へ actualに戻る。

## 3. 真のcompleted T38との全literal一致

入力は `MPS_SMOOTHING_COMPLETION_20261103_v62/T38.json` の全completed tensor。
座標を `(y0,y1,y2,y3,y4)` とし、全source parameterを

`p0=2, t2=−1, t3=1, t0=−k3, t1=−k1`

と置く。period_pfaffian_open=p0²/4=1なので sourceは全 k1,k3 上で定義される。
このT38の全三高次補正を残す。
この sourceの全非零entryは

```
P_y0,y1 = −k3*(y2*y4−y3²),
P_y0,y2 = −k1*(y1*y3−y4²),
P_y1,y2 = 2*(y1*y2−y3*y4),
P_y1,y3 = y1*y3−y4²,
P_y2,y4 = −y2*y4+y3².
```

actual coordinate permutation

`(y0,y1,y2,y3,y4)=(c,x0,x3,x2,x1)`

で、全十 quadratic entries は§2の Pbase／W とliteralに一致する。
特に raw P_{i,c}=−P_{c,i} の符号を含め W0=k3q2,W3=k1q0である。
Euler correction、tensor scale又はnormal-symbolだけの一致ではない。
従って whole L3 semisimple全sourceは既知closed PGL-stable成分T38に入る。

## 4. 全nilpotent L3 fiberを含むconstant-rank vector bundleの橋

nil代表 `(0,0,1,−1,0,0)` へ新しいpoint/profileを調べる必要は無い。
旧 V2 completekernel theoremは、L3のsemisimpleとnilpotent両非零orbitにおいて
元 full350×12 vertex systemのrankが正確に六と証明している。
actual SL2／scalar逆で全orbitを覆うので、
**irreducible base L3\{0} の全点で rank六**である。

係数がbaseにregularなこのmatrixのkernelは、したがってrank六vector bundle。
非零six-minorの開被覆上で直接linear inverseを取る標準constant-rank命題であり、
未証明flatnessを仮定するものではない。
その全totalspaceはirreducible、semisimple discriminant非零の逆像はopen dense。
semisimple部分の全imageは§3によりT38に入る。
Euler quotient後のsection空間でT38のaffine coneは閉、全sourceのtensor mapはregularなので、
全totalspace、従って**全nilpotent six-dimensional vertex fiber**のimageもT38に入る。
零sectionのprojectivization問題はaffine coneで処理し、対象非零rank四点を最後に射影化する。

rank-zero baseの十二fiberへこのvectorbundleのrankを外挿していない。
その別補集合は§1でrank≤2として処理した。

## 5. Q2全sourceはactual一次因子を持つ

全Q2非零軌道の代表 `(A,B,C,D,E,F)=(1,0,0,0,0,3)` は

`P12=q0, P13=(3/2)q1, P23=3q2`、全base P0i=0。

保存された complete three-dimensional vertex kernelを全 k0,k1,k2で書くと

```
W0=0,
W1=−(k2/3)*q0+(k1/2)*q1,
W2=(2k0/3)*q0+k1*q2,
W3=k0*q1+k2*q2.
```

従って全vertex coefficientに対し **full五変数Πのrow0は恒等零**。
全signed Euler Pfaffianは i=0項だけが残り、

`actualD=x0*Pf(Π restricted to indices1,2,3,c)`。

rank四でactualD非零ならproper homogeneous一次因子を持ち、被約既約quinticに矛盾する。
全actualSL2／source scalarのinverseを戻しても非零一次因子を持つ。
元の全Q2 point、全三vertex coefficient、全zero補集合を保持したwhole factorである。

## 6. 全branchの正確な結論

全base支持 L3∪Q2、base零、各全vertexkernel、実normalizer逆を戻すと、
指定rank零twisted cubic-cone・被約既約actualquintic sourceは全てT38に入る。
Q2／base零はその対象では空。
新成分の増分は零。
可約Dを持つsourceの全既知component分類、coneのrank二制限、他曲面、全E完成は本稿の結論ではない。
