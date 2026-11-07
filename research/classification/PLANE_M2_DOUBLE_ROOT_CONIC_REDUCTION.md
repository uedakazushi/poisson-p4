# reduced conic+line・m2 double root：横断配置は全空

2026-11-07。主義務G/E。rank二Poisson plane、実D既約被約quintic、位数二、
generic double normal quadraticを仮定する。zero cubicはreduced conic+lineとする。
新しいsource行列、minor、GB又はjetは計算しない。

**結論。** transverse conic+lineの全normal-double sourceは空。
tangent conic+lineはmoving rootの一parameter全P0と、下記constant-rootの二つの必要diagonal型へ還元される。
tangentの全higher sourceの空性、次元上界又は帰属はまだ主張しない。

## 1. 全leading形と全higher bad row

`PLANE_M2_DOUBLE_ROOT_CONNECTION_GEOMETRY.md` よりprimitive root degreeは零又は一。
degree一のtransverse配置は `PLANE_M2_MOVING_DOUBLE_ROOT_NORMAL_FORM.md` 第4節で全P0が空。
従ってconstant rootだけを新しく扱う。G2=h u^2と規格化すると、hの全可能性は
f=Lq又はL^3であり、theta-sigma=dlog(f/h)。
全off-diagonal one-form bは閉形式omega=(f/h)b=N2/hを満たすが、本証明ではbを固定しない。

canonical tracefree diagonal mixed fieldsをX,Yとする。conic+lineでは両者は同じ一dimensional
closed-scalar field spaceに属し、X(f)=Y(f)=0かつX wedge Y=0。
従ってactualD2のuv係数は3f kappa。targetG2のuv係数零からkappa=0。

全P1のu-rowのuで割れない項はv^2 C、Cはconstant base vector。
既約actualDにはC nonzeroが必要である。
全sourceのprojective Jacobiのnormalweight一係数は、base volumeでvectorに双対化すると

\[
 E\mathbin\times(DY\,C)-(2Y-X)\mathbin\times C+
 \nabla(C\mathbin\cdot\nabla f)=0,\qquad
 (4DY-DX)C=0.\tag{1}
\]

第二式は第一式のcurl。全off-diagonal b、残りのP1、全P2はこの係数へ寄与しない。
従って必要なCを第一式ごと排除すれば全higher sourceが空になる。

## 2. diagonal field差の正しいscale

raw mixed fieldとnormal connectionの対応は

\[
 i_E i_{X-Y}(dx\wedge dy\wedge dz)=-3f\,d\log(f/h).\tag{2}
\]

これはcanonical tracefree field gaugeのscaleであり、base Jacobian tensorはP_f。
以下のinteger eigenvaluesはこのscaleを用いる。

transverse f=z(xy-z^2)ではD=x partial_x-y partial_y、
tangent f=z(xz-y^2)ではN=2y partial_x+z partial_yとする。
式(2)を三成分で直接比較すると

| 配置 | h=f | h=z^3 |
|---|---|---|
| transverse | X-Y=0 | X-Y=3D |
| tangent | X-Y=0 | X-Y=-3N |

となる。hの全support choicesを尽くし、どのperiodもあらかじめ反転しない。

## 3. transverse配置の全constant-root sourceは空

まずh=f。X=Y=lambda Dと置く。curl式は3lambda D C=0。
lambda nonzeroならC=c partial_zであり、式(1)の第三成分には
partial_z(C(f))=-6czだけが残るのでc=0。
lambda=0なら式(1)はgradient(C(f))=0となる。
fがdegree三かつessential三変数なのでC(f)=0はC=0を強制する。

次にh=z^3。Y=lambda D、X=(lambda+3)D。
curl式は3(lambda-1)D C=0である。
lambda not equal oneなら同じC=c partial_zからc=0。
lambda=1の例外を落とさず、C=(a,b,c)と置いて第一式の全三成分を展開すると

\[
 (2bz-cy,\;2az-cx,\;2bx+2ay-6cz)=0.
\]

全係数からa=b=c=0。
従って全period、全off-diagonal source、全higher sourceでC=0となる。
全u-rowがu倍となるのでactualDはu因子を持ち、既約性に反する。
degree一の既排除と合わせて、**transverse conic+lineのm2 normal-double全sourceが空**。

## 4. tangent配置：二つの全必要diagonal型

h=fならX=Y=lambda N。lambda nonzeroではcurlからC=c partial_x。
式(1)はc(lambda+2)z=0なので、C nonzeroに必要なのは

\[
 Y=X=-2N,\qquad C\in\mathbb C^*\partial_x.\tag{3}
\]

lambda=0ではC(f)=0からC=0なので別枝は無い。

h=z^3ならY=lambda N、X=(lambda-3)N。
lambda not equal minus oneではcurlからC=c partial_x、第一式からlambda=-5。
lambda=-1の例外ではC=(a,b,c)に対する式(1)は

\[
 (cz,\;-4bz+2cy,\;cx+4az-4by)=0
\]

となりC=0。
従って必要なのは正確に

\[
 Y=-5N,\qquad X=-8N,\qquad C\in\mathbb C^*\partial_x.\tag{4}
\]

式(3)、(4)はfullhigher sourceの必要条件であり、off-diagonal b、normal u^2 coefficient、
他の全P1/P2を零へ置いていない。これらのsourceの存在・全支持は次義務。
moving degree一の場合には別noteの全P0(t nonzero)を保持する。

## 5. 依存と境界

canonical full30/9 filtration、root-connection geometry、全projective bad-row/curl恒等式だけを使った。
新しいsample point又はgeneric-only matrix rankから空性を推論していない。
分母はdense base chart及びfだけで、未知higher sourceを反転しない。
conicがdegenerateする境界はthree-linesの別R型であり、この定理へ混ぜない。
nonreduced R、m3、normal quadratic zeroは別義務として残す。
本稿は手証明で、全scale/例外lambdaを含む別作者の独立読解を要する。
