# 全 R19\N27 の射影軌道分類：一変数族と六つの追加軌道

2026年11月3日。ここでは全射影交叉の成分R19からN27を除いた **全ての点** を分類する。
R∩Nの軌道や全Poisson空間の軌道分類を主張するものではない。

## 定理

全てのP∈R\Nは

    P_Q=∂y∧(4x²∂x−2xz∂z+(2xw+Q(h,z))∂w)

という標準形を持つ。逆に任意の二次式Q=A h²+Bhz+Fz²の標準形はR\Nに入る。
全PGL5軌道は次の表で厳密に尽くされる。

|Q|パラメータ・同値|軌道次元|
|---|---|---:|
|h²+κhz+z²|κ²が完全不変量、κと−κは同値|18|
|h²+hz|一軌道|18|
|h²|一軌道|17|
|hz+z²|一軌道|18|
|z²|一軌道|17|
|hz|一軌道|17|
|0|一軌道|16|

すなわち一変数族に加えて六つの単独軌道があり、これ以外にR\Nの軌道はない。
この表はscheme厚みや境界R∩Nの記述とは別である。

## 1. proper全域coverからt≠0のsourceを全ての点へ使える

全CA/CB交叉の66表では、最後のbig relative flag case以外の全source像はNへ入る。
任意のP∈R\Nは全proper旗模型へ持ち上がるので、少なくとも一つの旗対を持つ。
そのrelative orbitは他の65ケースではあり得ず、最後のbig caseに入る。
従ってPは八係数族のPGL変換として実際に表示できる。

このsourceでt=q4=0ならCA側f=0であり、その全proper像AfはNに含まれる。
よってP∉Nならt≠0である。
既存のdeterminant1の全有限変換はt≠0だけを分母条件とするから、
A=0、F=0、Q=0であっても同じ変換でP_Qへ移せる。
稠密な族のclosureを任意のsource表示だと扱ったのではなく、全proper coverと全有限表を用いている。

## 2. 任意のQの標準形はNに入らない

canonical quotient fieldをD=5Cと書けば

    D=(-8hx,12x²,-18xz,5Q+2xw).

四成分の線形結合をhx,x²,xz,xwの係数で見ると、正確な行列は

    diag(-8,12,-18,2)

であり、全Qで可逆である。
従って非零な定covector Casimirはなく、二つのCasimirsを持つAf=N_fへ入らない。

N_gへ入るためには、独立な線形形式ell,kについて

    D(ell)=lambda*ell²,
    D(k)=−alpha*ell²−4lambda*ell*k

が必要である。
ell=a h+b x+c z+d wとする。
lambda=0は上の可逆行列に反する。
lambda≠0ではw²係数からd=0、続くh²,z²係数からa=c=0となる。
この議論はQ=0でも同じである。
従ってellはxの倍数で、lambda*ell=12xとなる。
第二条件の左辺D(k)+48xkは

    40a*hx+(60b+alpha)*x²+30c*xz+d*(5Q+50xw)

であり、x²の倍数にするにはa=c=d=0が必要である。
kもxの倍数になり矛盾する。従ってP_Qは全QでN_gにも入らない。
元のsourceにあるのでP_Q∈Rであり、これで正確にR\Nを覆うことが分かった。

## 3. Q≠0の全同値と五つのzero-pattern例外

Q≠0なら、既存の全割り算証明によりintrinsic factorLは一意で、quotient上の全invariant hyperplanesはh,x,zだけである。
xはself-cofactor lineであり、h,zはcofactor比−2/3,−3/2でラベル付けされ、入れ替えられない。
したがって標準形間の同値はquotient上の対角作用だけであり、三係数の変換則は

    A'=b*d*A/a²,
    B'=b*d*B/(a*c),
    F'=b*d*F/c².

A,B,Fのzero patternは保存される。
AF≠0の部分は既存のκ²による完全分類となる。
AF=0,Q≠0では可能なzero patternsが五つだけで、各非零係数は対角scaleで1へできる。
たとえばA,B≠0,F=0ならa=√A,c=B/a,b=d=1とし、Q=h²+hzに戻す。
B,F≠0,A=0も同様である。一係数だけ非零ならその一つを1へすればよい。
これで表の五つの単独非零例外が得られ、zero patternが異なるので同値ではない。

## 4. Q=0は別の一軌道

Q=0ではD=x*diag(-8,12,-18,2)*(h,x,z,w)である。
対角の四固有値は相異なるため、全invariant linear hyperplanesは正確にh,x,z,wとなる。
xが唯一のself-cofactor lineで、他のcofactor比は−2/3,−3/2,1/6と相異なる。
したがって四つの線も全てラベル付きで固定される。

intrinsic factorの一意性も保たれる。
もし定ベクトルu∈UについてE_U∧u∧D=0なら、xを除いた多項式恒等式から
U∧u∧diag(-8,12,-18,2)U=0となる。
二つの異なる座標だけを非零にすると、uはその座標二平面に入らなければならない。
互いに素な二つの座標平面を選べばu=0である。
同じ主張を元のcanonical175×5 factor行列でもrank4、kernel=C∂yとしてQQで確認した。

Q=0の標準形は固定された一つの点であり一軌道である。
非零Qは三つ、Q=0は四つのinvariant hyperplanesを持つので、両者が同値になることもない。

## 5. 軌道次元を全zero patternsで計算する

全normal formsでyのscale一つと四つのy-shearsが自由に作用する。
quotientの対角四変数について、非零な係数の固定条件の重みは

    A:(-2,1,0,1),
    B:(-1,1,-1,1),
    F:(0,1,-2,1)

である。非零重みのrankをrとすると、projective classを保つGL5群の次元は9−r、
PGL5 stabilizerは8−r、orbit dimensionは16+rとなる。
AF≠0または隣接する二係数が非零ならr=2、一係数だけならr=1、Q=0ならr=0である。
これで表の18,17,16が全て厳密に得られる。
GL5 dimensionはprojective class stabilizerの次元であり、affine canonical tensor stabilizerとは一つ違う。

## 証明書と範囲

* `code/ambient_R19_complement_orbits_20261103_v63.py`
* `reports/AMBIENT_R19_COMPLEMENT_ORBITS_20261103_v63/report.json`

全QでのCasimir・二つのN-flag係数、Q=0のintrinsic factor、全weight ranksとorbit dimensionsをQQ上で検査した。
測定0.11秒。R∩Nの全点の軌道分類、全P4 Poissonの分類をこの表から結論しない。
