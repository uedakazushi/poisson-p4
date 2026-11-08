# short nil31 r=2/3：全多項式線形化、rank二cofactor親の閉障害とT19の有限入口

2026-11-08、PRIMARY E、`/root/log_manuscript`。純手証明。
新しい点・profile・operator・rank・minor・jet・CASを実行していない。
旧source・証明書・主稿・frontier・catalogueは変更しない。
全shortnil r2/3の既知親への正の帰属は未解決である。

## 1. 全sourceと多項式automorphism

対象は `MODULAR_NIL31_SHORT_SUPPORT.md` §7の全六parameter源を、同節の
正逆affine変更で戻した全二parameter模型である。座標を `(X,Y,Z,W)` とし、
`t!=0,d` 任意を保持する。全六括弧は

```
PXY=dZ-tY+2tZW-Y²/3, PXZ=tZ-2YZ/3, PXW=d+YW/3,
PYZ=-2Z²/3, PYW=(Y+ZW)/3, PZW=2Z/3.
```

旧正規化によりtdunitは一つの22次元actual PGL/scalar軌道であり、d0は
同じ元既約多項式源内の閉境界である。本稿はd0を除外して別源に置かない。

次は開集合Zunitだけの有理座標ではなく、**A4全域の多項式automorphism**である。

```
F=X-YW+ZW²/2, Q=Y-ZW, Z=Z, W=W;
Y=Q+ZW, X=F+QW+ZW²/2.
```

Jacobianは1。全Leibniz展開から、変換後の全六式は

```
{F,Q}=-tQ, {F,Z}=tZ, {F,W}=d,
{Q,Z}=0, {Q,W}=Q/3, {Z,W}=2Z/3.                (A)
```

例えば `{F,Y}=dZ-tY+2tZW` と `{F,ZW}=tZW+dZ` を引くと
`{F,Q}=-tQ` となる。全式は分母無しであり、Z0も失わない。
(A)のJacobiは、F/Wが可換な対角derivationでQ/Zへ作用し、
`{F,W}=d`が定数であることから直ちに従う。全sourceのJacobiもこの正逆で戻る。
Pfaffianは `-tQZ` で、元sourceのactual五次は
`-t h²Z(hY-ZW)`（全体scalarは不問）となる。

このautomorphismは次数三を持ち、PGLではない。従って線形Poisson模型(A)を
既知成分へ同定するだけで、元projective tensorの所属を結論してはいけない。

## 2. 全一parameterの正則な逆向き退化

alphaを任意parameterとし、(A)から

```
Y=Q+alpha*ZW,
X=F+alpha*QW+alpha²*ZW²/2
```

という全A4多項式automorphismを用いる。その全pushは、元の座標で厳密に

```
PXY=alpha*dZ-tY+2alpha*tZW-alpha*Y²/3,
PXZ=tZ-2alpha*YZ/3, PXW=d+alpha*YW/3,
PYZ=-2alpha*Z²/3, PYW=(Y+alpha*ZW)/3, PZW=2Z/3.  (B)
```

全spatial次数は二以下、全alpha冪は非負。したがってこれはprojective quadratic
Poissonの実多項式曲線であり、全Jacobiと全Pfaffianは恒等式である。

```
Pf(B)=-tZ(Y-alpha*ZW),
D_alpha=-t h²Z(hY-alpha*ZW).
```

alphaunitではactual affine/PGL変更

```
X'=X/alpha, Y'=Y, Z'=alpha*Z, W'=W
```

により全六式は旧canonical源の `(t',d')=(t/alpha,d/alpha)` に一致する。
従ってtdunitでは(B)のalpha非零部全体が同じ既存22軌道へ戻る。
alpha0は(A)の線形模型であり、actualDは `-t h³YZ` である。
**この矢印はshortnilから線形模型への退化であり、既知親からshortnilへの正の矢印ではない。**

## 3. 線形端点はM11に入るが、唯一成分という橋は使えない

alpha0の模型を `P=partial_W wedge V+Jac(f)` と書くと、base座標(F,Q,Z)で

```
f=-tQZ,
V=-d partial_F-(Q/3)partial_Q-(2Z/3)partial_Z,
div(V)=-1, V(f)=-f.
```

wholeconstant-modular quadratic f=zw sourceIの一部そのものである。
`CONSTANT_MODULAR_SOURCE_ONE_M11_WHOLE_ACTUAL_LIMIT.md` と、正式な
`R4_MODULAR_RANK1_GZERO_M11_25_COMPONENT` の後続全成分採録を接続した
sourceI全点包含を使うと、全tunit/dでこの端点はM11に入る。
これは新しい一点評価や接空間検査を要しない。

しかし次節より、alphaunitの全shortnil軌道はM11の**全閉包外**である。
(B)の全既約像閉包はNoether性によりあるambient成分Cへ包含される。
CはM11でなく、alpha0端点は両者の交点に属する。
従ってこの端点を「既知M11の一点」と同定しても、全曲線の所属を
局所唯一成分性から戻すことは出来ない。実際この端点では少なくとも二つの
ambient成分が交わる。新しい成分の相違数や次元を主張したのではなく、
Cは未同定の既存成分でもよい。

## 4. 全rank一modular long/short親を同時に排除するproper閉条件

rank四開集合上で、標準divergence-free同次liftをbarPiとする。
marked hyperplane ellに対しell²がactualDを割る場合、Poisson tangencyから
`barPi#dell`はellで割れ、

```
barPi#dell=ell*B_ell,
tr(B_ell)=0, B_ell(ell)=0
```

というlinear vectorfieldが一意に決まる。最後の二式は、それぞれ
`div(barPi#dell)=0` とskew性から従う。ellのunit係数を一つ固定したdual-projective
chartではexact polynomial divisionがregularなので、B_ellはregularに依存する。

したがって、rank四開集合とdual P4の積内のmarked条件

```
ell² | D, tangency, rank(B_ell)<=2
```

は閉である。第一条件はproper multiplication incidence、第二は係数の割切り条件、
最後は全3minorの消失である。ellを忘れる写像はproperだから、その像
L2はrank四開集合内の閉集合になる。Pfaffian零やrank低下へは適用していない。
これは `RANK_FOUR_OUTSIDE_LOW_MODULAR_FLAG_DIMENSION_BOUND.md` のL(rank<=1)とは
別の、ここだけで使うrank<=2の閉flagであり、大域次元上界は主張しない。

全17rank一long及び全shortの**multiple-H/Pf(Q)=0**源では、
標識hのhomogeneous modular fieldはh partial_x+x partial_y、又は
h partial_w+x partial_yのrank二であり、canonical cofactorはその非零scalar倍。
よって各全sourceとそのPGL/scalar飽和閉包はL2に入る。
既知half/third/DDHとの重複した名前がある場合も、同じproperなmarked条件が
その全成分閉包に残る。単純Hだけを選んだ任意rank一modular源へは拡張しない。

一方、shortnilのactualD `h²Z(hY-ZW)` はtunitでhが唯一のmultiple hyperplaneであり、

```
N_h=Y partial_X+Z partial_Y+h partial_W
```

はrank三である。従って全tunit/dの対象はL2外。
**全17long/shortにより被覆した各rank一25成分のどの全閉包も、この対象の親にはならない。**
源のgeneric rankだけを比べたのではなく、unique repeated factorをproper incidenceで
特殊化して得る全親閉包の障害である。

(B)ではalpha0のcofactorがrank一、alphaunitではrank三となるので、
この閉条件と矛盾しない。既知親にあるのはalpha0だけであり、矢印の方向が逆である。

## 5. 真completedT19の二factor入口（未解決）

新T19のGa²曲線は別sourceを対象としており、その結果をshortnilへ外挿しない。
ただし `TWO_QUADRICS_NIL22_LAST_CORE_WHOLE_T19_CONFLUENCE.md` の真正sourceは
純幾何の入力として使える。座標(A,B,C,D,H)、period(r,s)、delta=(r-s)/2、
alpha=(r+s)/2、beta=(3r-s)/2で、smoothingを(r-s)/2,r-sに固定すると、
H-Casimir全六式は

```
PAB=rAB+delta*DH, PAC=beta*AC+2delta*BH, PAD=alpha*AD,
PBC=alpha*BC, PBD=-alpha*BD, PCD=-2rCD+2delta*B².
```

actualDは `-r(r-s)H q4 q3`、q4=AB+DH(rank四)、q3=CD-B²(rank三)。
rank四quadric Q=hY-ZWがtargetで既約なので、rank三以下のq3はQへ移れない。
proper factor flagを用いると、shortnilへの候補は次の二つに限られる。

```
outerH -> h: q3 -> hZ, q4 -> Q;
outerH -> Z: q3 -> h², q4 -> Q.
```

この二入口はtrueT17の同様の二入口と形が似るが、period関係と全tensorは異なる。
Ga²作用がq3/q4を保存しHを動かすことは、固定旧HのSS2陰性や固定SS12companionの
陰性では扱わなかったcoupled flagを与える。ただし正則な全tensor端点と
全sourceのE所属は未構成である。同じ旧profileやhigherjetは再実行しない。

## 6. 正確な結論

本稿で確定したのは全A4多項式線形化、逆向きの全parameter Poisson曲線、
その線形端点でM11局所唯一性を使えない理由、及び全rank一long/short親の
proper閉flag障害である。全shortnil r2/3の正のknown-parent曲線、全成分性、
任意MPS閉包からの非所属は未解決のままである。
