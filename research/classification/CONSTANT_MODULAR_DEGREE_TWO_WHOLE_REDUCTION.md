# 定数非零modular場・deg f≤2：全source有限化とCA/CBへの全族曲線

2026-11-08、PRIMARY E、producer plane_collision_review。
`MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §3、式(6)の未処理低次数枝を扱う。
**rank三quadratic等の全空性、rank一centered／degree一の全CA/CB帰属、残る三つの完全有限sourceを証明する。**
残三源の成分帰属はまだ主張しない。新operator、点、kernel、profile、GB、jetは実行しない。

## 1. 正確な全入口とnormalization

元全sourceは affine coordinates (t,u,z,w) で

```
P=∂t∧Y+Jac(f),  divY=κ≠0,  Y(f)=κf,
degY≤2, degf≤2.
```

Jacのorientationは(u,z,w)で

```
{u,z}=f_w, {u,w}=-f_z, {z,w}=f_u.
```

fの加法定数は元固有関数式を満たすものとして既に一意に固定されている。
新たにfへ定数を加えてcritical valueを捨てない。
非零tensor scale 1/κでY→Y/κ,f→f/κとしてκ=1に出来る。
以下はこの正規化を使う。

base affine変換x'=Gx+bでJac(f)はdet(G) Jac(f∘G^{-1})へ移る。
tの非零scaleと全tensorの逆scaleを同時に使えば、fの非零overall scaleも独立に調整出来る。
従って以下のquadratic正常化はactual affine/PGL変換とtensor scaleで戻る。
逆はG^{-1}、逆translation、非零scaleの逆であり、非線形Darboux写像は同値として用いない。

quadratic Hessianのrankをrとする。congruenceでqを標準化し、image内のtranslationでlinear項を消す。
radicalに非零linear項が残れば、その項をu座標に選びconstantもuへ吸収してf=q+u。
残らなければf=q+cで、cは全て保持する。r=2はq=zw、r=1はq=z²。
r=3は全linear項を除いた非退化q+c。r=0はdegree一／零へ戻る。

## 2. rank三とcritical value非零は全空

r=3のcentered形f=q+cで、Y=Y2+Y1+Y0を全homogeneous次数に分ける。
degree三からY2(q)=0、degree二からY1(q)=q。
非退化qに対する後者はlinear matrixのtraceを3/2と強制する。
一方divY=1のconstant項はtraceY1=1。divY2はhomogeneous一次なのでこの定数矛盾を補えない。
従って全r=3 sourceが空。Y2を零と仮定していない。

r=1,2でradical linear項が零ならf=q+c。qのcritical locus上でY(f)=0なので、固有関数式はc=0を強制する。
よってc≠0は全空、c=0だけ次節へ残す。
fがdegree零なら同様にf=0だけで、Jac(f)=0ゆえrank≤2。rank四sourceは無い。

## 3. 残る三つの完全有限source

以下でY=A∂u+B∂z+C∂w。各parameterは任意の複素数で、sourceの追加unitを反転しない。
表示は必要supersourceだけでなく、各全値でdivY=1,Yf=fを満たす完全sourceである。

### (I) f=zw：全十係数

```
b1=αu+βz+γw+δ,
B=z*b1, C=w*(1-b1),
A=u*(-βz+γw)+K2(z,w),
```

K2は次数高々二の全六係数。b1四係数と合わせ十。
wB+zC=zwとgcd(z,w)=1からB=z*b1,C=w*(1-b1)が全解。
divY=1はA_u=-βz+γwなので、全Aが上式となる。
K2のconstant、linear、zw、square項を一つも捨てない。

### (II) f=zw+u：全八係数

```
L=αu+γz-βw,
B=z*L+βu+λz+ρw+b,
C=-w*L+γu+σz-λw+c,
A=u+zw-βuw-γuz-ρw²-σz²-bw-cz.
```

全parameterは(α,β,γ,λ,ρ,σ,b,c)。
Yf=fからA=f-wB-zC。Aのdegree≤2はwB2+zC2=0を強制し、
gcdからB2=zL,C2=-wL、Lは全homogeneous一次となる。
B1,C1,B0,C0を全て残してdivを比較すると、Lのz係数=C1のu係数、
Lのw係数=-B1のu係数、C1のw係数=-B1のz係数を得る。
これが上の全八係数であり、α零や全lower項も含む。

### (III) f=z²+u：全十係数

```
B=βu+γz+δw+e,
A=u+z²-2z*B,
C=(2βz-γ)w+K2(u,z).
```

Bの四係数とK2の全六係数。
Yf=fはA=f-2zB。degree≤2からB2=0が必要。
divY=1はC_w=2βz-γで、これが全Cを与える。
K2のu²、uz、z²、linear、constantを全て保持する。

各式から元parameterはB,Cの指定係数又はAのu零切断として一意に読める。
係数復元に新しい分母／unitは無い。全零境界も同じaffine sourceに含む。
Yf=fとdivY=1によりLie_Y Jac(f)=0、従って全Pはstrict Poisson。
四変数Pfaffianは正確にYf=fなので、非零fのこれら全sourceは一般rank四である。
これら三源のactual既知親／成分性は本稿で未解決のまま残す。

## 4. f=z²の全十六potential sourceはCB

Yf=fからB=z/2、div条件はA_u+C_w=1/2。
C[z]係数の二変数polynomial Poincare補題により、全解は

```
Y=(u/4+H_w)∂u+(z/2)∂z+(w/4-H_u)∂w,
degH≤3, H modulo C[z].
```

Hの全二十monomialからzだけの四項を除いた十六係数であり、
元の任意A,Cの全解に同時potentialが存在する。
これはsource係数の積分表示であって、非線形座標変換でのPoisson同値ではない。

known rawCBを(source x,y,z0,w0)=(z,t,u,w)、potential K=-Hで書くと

```
P_zt=g(z), P_uw=μ*g(z),
P_tu=-g'(z)*u/2-K_w, P_tw=-g'(z)*w/2+K_u,
```

他二括弧は零。
g=-z/2, μ=4を入れると、元Pの全六括弧とrawのまま一致する。
Euler補正、除いたsource係数又は局所化は無い。

任意の同じHに対してg_ε=-z/2+εz²、μ=4とする全二次族もstrict Poisson。
元Jacobiの唯一の非自明な三添字(t,u,w)では、
`-μ*g*g'+μ*g*(g'/2-H_uw)+μ*g*(g'/2+H_uw)=0` と相殺する。
ε≠0ではg_εは二つの相異なる根を持つ非退化一変数quadratic。
affine shift/scaleでknown CBの標準形へ戻り、変換後も全cubic potentialを保持する。
従ってε≠0の全族はCB、ε=0もactual polynomial curveと閉性によりCB。
この証明は全十六potential／全元frameに同時に適用し、zero potentialも残す。

## 5. degree一f=zの全十六sourceはCA

fが非定数linearならactual affine変換でf=z。
Yf=fはB=z、div条件はA_u+C_w=0。
従って全解は

```
Y=H_w∂u+z∂z-H_u∂w,
degH≤3, H modulo C[z].
```

Jac(f)=-∂u∧∂wで、known rawCAを同じcoordinate順、potentialK=-Hで

```
P_zt=g(z), P_uw=-1,
P_tu=-K_w, P_tw=K_u
```

と書けばg=-zで全六括弧が一致する。
g_ε=-z+εz²に変えた全sourceは、Hを固定したまま全εでstrict Poisson。
ε≠0ではquadratic leading係数が非零なのでCA標準のx²/2-Kへactual affine shift/scaleで戻せる。
全族は閉CAへ入り、中心の全H sourceもCAである。

## 6. 全補集合と残余

| 元の全degf≤2 branch | 結論 |
|---|---|
| Hessian rank三 | 全空 |
| rank二、radical linear非零 | 完全全八source(II)、帰属未決 |
| rank二、radical linear零、critical value非零 | 全空 |
| 同critical value零 | 完全全十source(I)、帰属未決 |
| rank一、radical linear非零 | 完全全十source(III)、帰属未決 |
| rank一、radical linear零、critical value非零 | 全空 |
| 同critical value零 | 全十六source、全CB |
| degree一 | 全十六source、全CA |
| degree零 | 全rank四空 |

元κ、全affine normalizer、Jacのdet補償、全critical valueと全parameter零を戻してこの表が得られる。
滑らかなf3の既受理MPS theoremや特異f3の未解決枝へ本結果を拡張しない。
rank三quadratic空性から、別modular場を持つCA/CB族の一般εを排除したと誤解しない。
CA/CB変形はconstant-modular laneを離れることが許される実Poisson曲線である。

本稿は純手証明。有限pure-proof jobのgate／hash確認は数学の証明と区別する。
新算術、global files、主稿、Gitの変更は無い。全CA/CB actual familyと有限源の別作者査読を依頼する。
