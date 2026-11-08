# Ga9 linear-w：固定退化因子の全源と唯一の磁気係数

2026-11-08、PRIMARY E、純手証明。新しい点・matrix・operator・profile・jet又はCASは実行しない。
**固定因子の全Poisson源を二軌道へ還元する定理**であり、T36/T37へのpositive membership又はGa9の成分性はまだ主張しない。

## 1. 固定した被約因子と実cotangent源

座標は(h,x,y,z,w)、baseは(h,y,w)、fiberは(x,z)。

```
F=4h³w−8h²y²+y⁴, D=−4hF,
U=(4h²−y²)∂y+(16hy−4wy)∂w,
V=−hy∂y+(−4hw+4y²)∂w.
```

FはC[h,y]上のprimitiveなw一次式なので絶対既約であり、hと互いに素。従ってDは被約である。全二次R(h,y,w)に対し

```
P_R=∂x∧U+∂z∧V+(−hx+2yz+R)∂x∧∂z                 (1)
```

と置く。直接の微分で

```
[U,V]=−hU+2yV, U(D)=−4yD, V(D)=−4hD,
det(E_base,U,V)=D                                             (2)
```

を得る。(2)の最初の恒等式により、(1)は任意の全六係数Rでstrict Poissonである。signed Euler-Pfaffianはdet(E_base,U,V)であり、Rには依存しない。R=4w²が凍結されたGa9 canonical四倍tensorと全十括弧で一致する。

## 2. この固定Dを持つ全rank四tensorの尽くし

固定Dとはpolynomialの非零定数倍までを意味する。まずactual PGLで因子を上の座標へ置き、既済BV／Cartanのcanonical liftを取る。このliftはEuler gaugeに依存せず、Π#dD=0を満たす。ここで用いるのは `ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md` §1 のcanonical入口だけで、同稿後半の**既約**quintic分類は適用しない。

baseのEuler場をE_bとし

```
A=U+(4y/5)E_b, B=V+(4h/5)E_b.
```

と置く。A(D)=B(D)=0、det(E_b,A,B)=Dである。被約Dの各素因子の一般smooth点で、logarithmic fieldのCramer分子はその因子で割れる。Dがsquarefreeなので全Dで割れ、E_b,A,Bはlogarithmic derivation moduleのbasisとなる。Euler分解により

```
AR(D)=C[h,y,w] A ⊕ C[h,y,w] B,
AR_0=AR_1=0, AR_2=span_C(A,B).                              (3)
```

これはdeterminant criterionの直接割算証明である。

canonical liftのbase-base skew blockは∇Dを殺す。三変数のskew matrixなので、その係数vectorは∇Dの有理関数倍。被約性から∇Dの成分はprimitiveであり、その倍はpolynomialとなる。しかし係数のdegree二と∇Dのdegree四が両立しないため、このblockは零である。

mixed base vectorfieldsも∇Dを殺す。fiber(x,z)で展開すると、fiber次数一／二の係数はそれぞれAR_1／AR_0に入り、(3)により零。従ってmixed部分はA,Bの**定数**二行frameで尽くされる。その2×2 frame Cがsingularならbivector squareは零となるので、rank四ではdet C≠0である。全entryの零を許したactual fiber変換C^−1=adj(C)/det Cで、このframeを(A,B)へ戻せる。

次にEuler gauge

```
P ↦ P+E∧((4y/5)∂x+(4h/5)∂z)
```

を取るとmixed部は(U,V)、h-rowは零となる。h-Casimir liftではprojective Jacobiのhを含む四添字式がh×strict Jacobiなので、strict Jacobi自体が零である。全fiber-fiber係数をKとすれば、baseを含む三添字式と(2)から

```
K_x=−h, K_z=2y.
```

よってK=−hx+2yz+R(h,y,w)。fiber二次項や未処理のinhomogeneous項は残らず、全degree二Rを保持した(1)に戻る。

以上により、**actual divisorがこのDとPGL同値な全rank四projective Poisson tensor**を、actual fiber GL2とEuler gaugeの正逆を保って(1)へ還元した。canonical frameの特定のentryを反転するchartを選んでいない。

## 3. 全磁気項を除くactual shearと明示逆

```
x′=x+l, z′=z+m,
l=a h+b y+c w, m=d h+e y+f w
```

はdet一のactual GL5であり、逆はx=x′−l,z=z′−m。base-base bracketが零なのでmixed部は不変で、fiber bracketは

```
R′−R=h l−2y m+U(m)−V(l)
 =(a+4e)h²+(2b−2d+16f)hy+5c hw−(4c+3e)y²−6f yw.        (4)
```

w²項は変化しない。所望の差を

```
Δ=r_hh h²+r_hy hy+r_hw hw+r_yy y²+r_yw yw
```

としたとき、次の**全域の線形逆**が(4)を実現する。

```
c=r_hw/5, f=−r_yw/6,
e=−(r_yy+4c)/3, a=r_hh−4e,
b=0, d=(16f−r_hy)/2.                                    (5)
```

差を−(R−r w²)と取れば全Rをr w²へ除ける。分母は非零有理定数だけであり、零の磁気係数も全て含む。shearのkernelはl=b y,m=b hで、canonical Ga stabilizer y∂x+h∂zと一致する。

さらにκ≠0についてfiber同時scale(x′,z′)=(κx,κz)とtensor scalar1/κを組み合わせると、(U,V)とfiber線形部を保ち、RをκRへ移す。逆はκ^−1。

従ってr≠0ならκ=4/rでcanonical Ga9へ戻り、r=0ならR=0へ戻る。これで**全固定D源は二つのactual PGL軌道**となる。r=0は別に扱い、unit式を零へ外挿しない。族R=4s w²は零軌道をunit軌道閉包へ置くが、その逆の包含は推論しない。

## 4. 唯一の残係数はintrinsicなplane上の階数

Dのreduced singular supportは唯一のplane

```
S={h=y=0}.
```

実際、D_w=−16h⁴からsingular点でh=0、次にD_h=−4y⁴からy=0が従い、逆包含は全partialの代入で確認出来る。Sへの(1)の制限は

```
P_R|S=r w²∂x∧∂z.                                      (6)
```

従ってr≠0とS上の一般Poisson階数二、r=0とS全体のrank零が同値である。raw一成分だけを射影階数と混同せず、w≠0 chart上の二つの独立方向∂x,∂zで(6)を読む。SはDから一意なのでこの二分は全PGL変換で保たれる。

これにより、固定divisorへの一致だけではGa9を認識出来なかった問題が、**divisorの一致＋unique S上のrank二**という正しい全源条件へ縮約される。

## 5. T36の全period実Dと今後の成功条件

completed T36の座標を(A,B,C,u,v)、p=p0≠0とし

```
q+=(p/5)Bu+t2v², q−=(p/5)Cv−t3u²
```

と置く。全correction `P_BC=pBC+10t2t3uv/p` を含む元八括弧から、Euler-Pfaffianの五項を手で展開すると

```
D_T36=3A q+ q−+t0u³ q−−t1v³ q+.                         (7)
```

すなわちA係数は `(3/25)(pBu+5t2v²)(pCv−5t3u²)` であり、残項は

```
−p t1 Buv³/5+p t0 Cu³v/5−t0t3u⁵−t1t2v⁵.
```

全t0,t1,t2,t3を保持しており、どれかを初めから零に置いた境界族だけをT36全体と呼ばない。この恒等式は純手展開で、新しいtool算術receiptとは区別する。

T36又はT37のactual sourceから、正則非零projective tensorを与える一曲線P_εを構成し、そのendpointの実DがDとPGL同値、かつunique S上でrank二であることを証明出来れば、§§2–4の全逆によりendpointはGa9軌道に入る。known成分の閉性とPGL不変性からGa9全源閉包がその親へ入る。

**tensorの正則性は依然必要**である。divisorだけのprojective limitを合わせても、tensorのleading termがrank二でPfaffian零となる可能性がある。全係数の負冪消去又は同値なglobal正則性証明を省略しない。rank零のSを得た場合もGa9 unit軌道への逆推論は出来ない。

現在はこの曲線又はglobal閉条件によるT36/T37判定は未解決。固定preconjugationの1PSは既Ga stabilizer障害により不適切であり、旧六normal二次支持、高次jet又は固定N profileの反復を再開していない。
