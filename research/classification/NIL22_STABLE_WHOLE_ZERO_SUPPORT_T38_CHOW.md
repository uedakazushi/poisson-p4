# Stable二全源の全零支持と、degree三Chow cycleによるT38障害

2026-11-08、PRIMARY E。全六磁気源を保つ純証明、新しいCAS・point・operator・profileは無い。
これはstableのE帰属又は成分性ではない。全T38からの以下の必要条件だけを与える。

## 1. Hyperplane h零上の全支持

既全源の座標を(h,q1,q2,x,y)、磁気を
`R=m0+m1*x+m2*y+m3*x²+m4*x*y+m5*y²`とする。
h零で全sourceの括弧は共通に

```
Pq1q2=m3*x²+m4*x*y+m5*y²,
Pq1x=-x²/3, Pq1y=-xy/3, Pq2x=-xy/3, Pq2y=-y²/3, Pxy=0.
```

Q=q2*x-q1*y+3(m3*x²+m4*x*y+m5*y²)と置く。
projective零はraw P零ではなくE wedge Pで判定する。
h零の四cubicは正確に `(xQ/3,yQ/3,0,0)`。
x=y零でもQ零なので、h零上のreduced零支持は正確にS2=(h=Q=0)。
Qの全gradientが零ならx=y=q1=q2零となるため、全磁気でS2はsmooth quadric surfaceである。
vertex line ell=(h=x=y=0)もS2に含まれる。

## 2. Affine chartの全補集合

h一ではh-Casimir代表なので、projective零は全四変数affine括弧の零そのもの。
全mixed anchor Mは二つのordinary projective vector fieldsである。
したがってM零の点は、既 `NIL22_STABLE_PROJECTIVE_ACTION.md` の全Cartan／Borel actionの共通固定点に尽くされる。

three_linesではCartanの全共通固定点は

`p12=(-3,3), p13=(-3,-6), p23=(6,3)`。

V0零なのでF=Rはqに依存しない。
各p_iでR(p_i)非零なら零点は無く、R(p_i)零なら全q1,q2 fiberが零。
従って**全磁気値**に対して

`Zero_red(P)=S2 union (R(p_i)=0 の各plane Pi)`、
`Pi=(x=x_i*h,y=y_i*h)`。

磁気商をRnormal=k0*x²+k1*x*y+k2*y²に戻した座標では

```
r12=9(k0-k1+k2), r13=9(k0+2k1+4k2), r23=9(4k0+2k1+k2).
```

この三つはdet=-3^9の全linear coordinates。
逆は `k0=(2r12-r13+2r23)/81, k1=(-5r12+r13+r23)/81, k2=(2r12+2r13-r23)/81`。
source parameterを反転せず、三つのzero hyperplanesと全交叉を保持する。
元六mの値も単にR(-3,3)、R(-3,-6)、R(6,3)であり、全translation商と整合する。

line_conicではregular nilpotent Kのprojective固定点は唯一、Borelのsemisimple場もそこを固定する。
全共通固定点はbase cubicのconic Q_baseとline Lの接点p=(6,-6)だけ。
F=R-(3/2)(q1+q2)なので、そのfiberの零集合は全磁気で一本のaffine line。
従って全reduced零支持はS2と次のlineの合併に等しい。

`x=6h,y=-6h,q1+q2=(2/3)*R(6,-6)*h`。

そのh零端点はell上にある。
とくにline_conic全源とthree_linesのr12*r13*r23非零稠密開では、唯一の二次元零支持はdegree二のS2だけ。

## 3. T38のdegree三零cycleを特殊化する

既completed T38のgeneric全源は、degree三twisted-cubic cone上でprojectively零。
これは全十entryがその三quadricのidealに入るliteral式で確認済みであり、
`CUBIC_CONE_ZERO_T38_SEPARATE_REVIEW.md` と
`GA9_LINEAR_W_WHOLE_T38_AND_A5_EXCLUSIONS.md` §2の全source式を使う。

targetがclosed T38に入ると仮定し、genericがその全sourceの稠密開にあるDVR curveを取る。
generic点のsource frameとdegree三coneは有限field extension後に選べる。
degree三effective二cycleのChow空間はprojectiveなので、cycleはspecialへ延長する。
generic cycle上でtensorが零だから、そのhorizontal closure及びspecialization cycleのsupportでもspecial tensorは零。
非reduced multiplicityは正整数としてcycleに残り、二cycleのdegree三は失われない。
normalization-flatness、zero scheme自体のflatness又はreduced limitは仮定しない。

しかしtargetの唯一の二次元零支持がS2なら、非零effective二cycleはm[S2]、m正整数。
degreeは2mであり、三にはならない。
ゆえに

* stable_line_conicの**全六磁気源**はwhole T38の外。
* stable_three_linesのT38所属は `r12*r13*r23=0` の必要条件を持つ。
  三値全unitの稠密開はwhole T38の外。

zero hyperplanes上でPiが現れればS2+Piのdegree三cycleは次数上可能だが、所属は主張しない。
既T40 graphの三値全零境界にはS2と三planeが存在し、degree五cycleを持てることとも整合する。
Chow次数だけでT40の精密所属を再証明したとは呼ばない。

## 4. 残る正の義務

stable源はP2上の `h*T(-log C3)=O(-1)^2` の全linear Lie-Poissonに、全H0(O2)磁気を加えたものと読める。
C3はtriangleならCartan、tangent line-conicならBorelのfree cubic。
このwhole algebroidと§2の全磁気値を保つactual parent deformation又は全成分性は未解決。
T21などの別whole親を本稿で排除していない。T36/T37だけが残るとは断言しない。
旧fixed-P1／gap一 normalconeの反復、又は磁気一点だけのzero-support比較へ戻らない。
新component、全E、global／main／Gitの変更は無い。
