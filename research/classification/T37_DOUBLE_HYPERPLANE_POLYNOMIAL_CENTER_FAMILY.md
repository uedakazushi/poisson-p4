# 真T37の一つのdouble-H polynomial曲線とconstant DOUBLEへの正しい向き

2026-11-08、PRIMARY E、`/root/quartic_sources`。純手の全係数導出。独立exact readback、新operator、point、jet又はprofileは未実行。以下の全actual族は新しいfinite sourceとして保存するが、constant DOUBLEの陽性親認識はまだ主張しない。

## 1. 真sourceと実frame

`reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T37.json` の真completed tensorを使う。source座標は(H,a,b,B,C)、p0=1、四smoothingは

```
(t0,t1,t2,t3)=(p,p,1/p,1).                         (1)
```

p≠0では全period/smoothingがunitで、named T37内にある。全correctionを保持し、sourceの十entryは

```
PHa=p*B², PHb=p*C², PHB=PHC=0,
Pab=a*b+p*B*C/18-H²/48,
PaB=3aB+H*b/p-C²/3, PaC=-aC+H*B/3,
PbB=bB+H*C/3, PbC=-3bC+H*a+p*B²/3,
PBC=-2BC-H²/(2p).
```

delta≠0、alpha=1/(144delta)とする。target座標(h,f,g,z,w)からsourceへの実linear mapを

```
H=h, a=alpha*h+2delta*w+p*f,
b=delta*h+p*g, B=z, C=w                           (2)
```

とする。detはp²、全逆はh=H,z=B,w=C,f=(a-alpha*H-2delta*C)/p,g=(b-delta*H)/p。tensor scalarはpである。非線形の座標変更は使わない。

## 2. 全十の正確な多項式

(1)を(2)で実pushしtensor p倍すると、全括弧は

```
Phf=p*z², Phg=p*w², Phz=Phw=0,
Pfg=3delta*h*f+alpha*h*g-4delta*w*g+z*w/18
      +(5delta/3)*z²-alpha*w²+p*f*g,
Pfz=h*g+3alpha*h*z+2delta*z*w-w²/3+3p*f*z,
Pfw=-alpha*h*w-2delta*w²+h*z/3-p*f*w,
Pgz=delta*h*z+h*w/3+p*g*z,
Pgw=alpha*h²-delta*h*w+p*(h*f-3g*w+z²/3),
Pzw=-h²/2-2p*z*w.                               (3)
```

例えばPfzのp極は `H*b/p-2delta*H²/(2p)` によってちょうどh*gとなる。Pfgの極の分子は

`a*b-H²/48-6delta*b*C+2delta*H*a`

で、p零係数は `(3alpha*delta-1/48)h²=0`。残全p係数は(3)の表示である。これらは全p,deltaの恒等式であり、genericの一parameter点列やjetの打切りではない。

したがってraw global150係数自体がp多項式である。真sourceからgeneric pでprojective Jacobiを満たすので、全pでその恒等式が延長される。p零でh行が零となり、h一restrictionから中心はstrict Poissonでもある。必要ならGL-equivariant canonical divergence-zero sectionを同じ全族に用いられるが、generic hをCasimirとは置かない。

signed Euler-Pfaffianの変換係数はscalar二乗p²とframe forward det p^-2の積一である。従って全actualDは真sourceのDへ(2)を代入したものに等しい。これは将来の同一曲線readbackで全126係数を確認出来る。

## 3. 中心の全係数と型

p零は(3)から最後のp項を全て除いた六式である。発散の規約 `N^j=sum_i partial_i Pij` により全cofactorは

```
N=(-3alpha*h+6delta*w)*partial_f+3delta*h*partial_g.   (4)
```

そのaffine linear部はrank一nilpotent、translationはimageの外にある。Pzw、Pgz、Pgwがそれぞれh²、h、hで割れ、Pfgのf依存もh*fなので、全Pfaffianはhで割れる。actualDは従ってh²で割れ、二simple H+cubicの第一二plane源には適用出来ない。このwrong-target判定には新演算が不要だった。

逆に、この全中心族はgeneric pの真T37からの正則非零中心であり、閉性によりT37の境界として使える。ただし本稿時点では手恒等式を保存した段階で、独立arith receiptとglobal台帳への登録は別である。

## 4. delta中心だけをconstant DOUBLEへ潰す向きは不可能

delta≠0に対するactual対角変更

```
(h,f,g,z,w) -> (h,delta⁶*f,delta⁴*g,delta³*z,delta²*w),
tensor scalar delta^-5
```

を全(3)へ掛けるとdelta=1、p'=p/delta⁵の同じ族となる。特に中心p=0は全delta unitで一つのactual PGL/scalar軌道である。

(4)は実非零stabilizerなので、その射影PGL軌道は高々23次元。一方、現在のconstant DOUBLEの全源は既受理で一つの23次元軌道であり、h-Casimir cofactorは−h∂t、affine linear部零である。唯一のdouble Hが内在的なので、linear部rank一と零はPGL/cofactor scalarで同値にならない。

異なる二つの同じ次元の軌道の一方が他方の閉包に入ることは出来ない。より大きい軌道の閉包にも、その最大次元以下を越える軌道は入れない。従って**この中心軌道だけをdelta→0/∞へ実PGL退化してconstant DOUBLEへ入れる方法は不可能**である。このorientationを逆に用いない。

また中心でPfgへ単にsigma*h²を追加する場合は、f'=f+sigma*h/(3delta)で除ける。他の五括弧はfに依らないので、これも有限actual shearの全式である。単純なfirst period derivativeで得るこの一項を、新しい中心modulusと呼ばない。

## 5. 残る正確な陽性義務

constant DOUBLEのcanonical六式は(t,x,y,r)順に

```
(1+x-2xy, 4y²-y, (1-6y)r+x, 0, -x², 1+2xy),
D=h²*(x²*y+x*h²+h³), cofactor=-h*partial_t.
```

全D101/M48から射影orbit23、接空間上界24が既に受理されている。25以上の親は不可能で、T37は既残存24親の一つである。

陽性を得るには、p≠0の真T37を保った**同時**period/moving actual frameの正則非零中心としてこのDOUBLEを作る必要がある。p零の中心を先に固定してからdeltaだけを動かす経路ではない。sourceHをgeneric Casimirとしたり、rawEulerの個別項のpoleを勝手に除いたりしない。全true period coupling、actual両逆、全六とPf、必要なparameter零補集合を保持する。

本稿だけではその同時curveもDOUBLE所属も証明していない。新power ansatz列、点、局所障害、jetを実行する契約は作っていない。Ga9 cuspidal等の別既知親の検討は別担当と分担する。
