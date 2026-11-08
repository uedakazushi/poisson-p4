# t=−3 canonical quinticのexact X(1,2,2) normalization

2026-11-08、`/root/quartic_review`。純幾何補足、新しい演算job無し。
`PLANE_TANGENT_M2_MINUS_THREE_WHOLE_ORBIT.md`のcanonical P(0,1)と実Fを対象とする。
結論はnormalizationと全projection matrixの具体化だけで、T36/T40等への帰属はまだ主張しない。

## 全plane familyと射影bundle map

u!=0上でρ=v/u、A=a−cρ²、B=b+cρ+uρ³/3と置くと

    F=3u²[u A²−vAB−3cB²].

B!=0でη=A/Bを取れば、F=0から

    3c−η²u+ηv=0,
    3a−3ηb−η³v=0.                                    (1)

を得る。第二式はA−ηBの式へ第一式を代入したもの。
η=s/tとしてhomogenizeすると、次の全linear plane projectionを得る。

    [a:b:c:u:v]
      =[s x : t x−s²z : s²y−stz : 3t²y : 3t²z].       (2)

sourceはX=P_{P1}(O(1)⊕O(2)⊕O(2))で、x列はbase degree一、y,z列はdegree二。
全三列は

    g1=(s,t,0,0,0),
    g2=(0,0,s²,3t²,0),
    g3=(0,−s²,−st,0,3t²).

t!=0ではg2,g3のu/v成分が独立し、g1も独立。
t=0では(a,b,c)の三座標に独立となる。s=0も最初の判定に含む。
従って全(s:t)でrank三、sourceにbasepointは無い。
式(2)を元Fへ代入すると恒等零であり、imageは既約hypersurface F=0に入る。

tautological Hは、三bundle summand全てpositiveなのでample。
basepoint-free map(2)はHをO(1)のpullbackに持つからcurveをcontractせず、properなのでfinite。
H³=1+2+2=5、imageの既約quintic degreeも5なのでgeneric degree一。
smoothなXからのfinite birational mapであり、(2)はF=0のnormalizationである。
u!=0/B!=0で導出した(1)から、u=0やB=0を失っていない。
全射影mapのbasepoint-free/finite性が全chartを含むことを別に証明したためである。

## cusp gluingのliteral first datum

target三座標をU=(a,b,c)、残二をW=(u/3,v/3)と書くと

    W=t²(y,z).

従ってprojectionのCP quadraticはdouble t²。
node CP=stを無断に代用せず、cusp/special caseを保持する。
s=1、t=0近傍でU=(x,tx−z,y−tz)。
U0=(x,−z,y)をfiber座標に取れば

    ∂tU|0=N U0,
    N=[[0,0,0],[1,0,0],[0,1,0]].

Nはregular J3であり、noncyclic gluingと同一視しない。
このfirst datum自体からcyclic nodal族へのactual deformation存在を推論しない。

## rankzero surfaceの別の具体式

canonical u1 chartではSing FはA=B=0であり、全affine Poisson bracketsはその上で零。
その射影閉包S5は

    [a:b:c:u:v]
      =[3CV²:−3CUV−V³:3CU²:3U³:3U²V].                 (3)

P2上のcubic systemには一つのdouble basepointがあり、blowup後H=3L−2EでH²=5。
rankzero S5のnormalizationはF1で、Dそのものの三fold normalization Xとは区別する。
S5のa1 chartではlocal source(t=U/V,v0=V/C)から

    b/a=−t−v0/3, c/a=t², u/a=t³v0, v/a=t²v0

となる。B=t+v0/3と書けばimage ringはC[B,t²,t³]、すなわちline×cusp。
このsurface/zero geometryだけで既知Poisson成分へ帰属させない。

## 正のknown-family認識の正確な残義務

別担当のX122全source routeへ、(2)の**full cusp projection**とcanonical tensorを渡した。
generic CP=stの真completed T36とのactual同一性、及びこのfull cusp mapへの
actual cyclic family／source specializationが必要である。
CPだけを動かしてもgeneric projectionがPoissonを持つとは限らない。

もしそのactual familyが得られ、かつ本Fのcanonical L_F（全Casimir＋divergence条件）kernelが
表示rank四lineだけなら、projective properness＋閉linear incidenceでpole-free limitを示せる。
これは全special sourceのrelative-descent定理より弱い十分条件だが、
本稿ではL_Fを新計算しておらず、そのline一意性やfamily存在を未証明のままにする。
ambient又はprofile系列は追加していない。
