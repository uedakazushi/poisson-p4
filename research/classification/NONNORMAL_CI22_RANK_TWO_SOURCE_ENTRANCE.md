# non-normal integral CI22 rank二：全五surface orbitとfull150入口

2026-11-07、主義務G、producer plane_collision_review。
normal CI22は別packetで全rank二restrictionの既約quintic枝を排除した。
ここではnon-normal integral非退化CI22を扱い、normalの場合のOS(1)延長は仮定しない。

## 1. 全pencilの五surface orbit

非退化memberがあるときは既受理のQ-self-adjoint Jordan全型を使う。
各固有値のeigenspaceは二次元以下である。
normalでないことはsingular locusがcurveを持つことであり、これは
ker(M−λ)二空間上でq1が恒等零、すなわち同じ固有値の二block長が双方二以上の場合だけ。
size五では(22)+(1)又は単一(3,2)に限る。
固有値はそれぞれ二つ又は一つで、pencil PGL2により連続parameterを残さない。

全member退化のcanonical coverでは、M1+regular2の唯一integralモデルと、
common radical一のbase integral elliptic quartic C⊂P3のconeだけが残る。
後者がnon-normalならCがsingularである。
四次元regular pencilでは各固有値のblockは一つだけ（さもなければrank二quadric）。
Cはarithmetic genus一なので、二つのsingular pointは正規化genusの非負性に反する。
Jordan4は標準q1=2x0x3+2x1x2,q2=2x1x3+x2²がline{x2=x3=0}を含みintegralでない。
従って全singular Cのblock型は(2,1,1)又は(3,1)、node又はcuspの二orbitである。
三又は二固有値はpencil PGL2で固定出来る。

以上の全五モデルを次の表で固定する。

| model | q1 | q2 |
| --- | --- | --- |
| regular22_1 | 2x0x1+2x2x3+x4² | x1²+x3²+x4² |
| nil32 | 2x0x2+x1²+2x3x4 | 2x1x2+x4² |
| singularM1 | x0x1+x3² | x0x2+x4² |
| cone_node | 2x0x1+x2²+x3² | x1²+x2²+2x3² |
| cone_cusp | 2x0x2+x1²+x3² | 2x1x2+x3² |

各々integralである。regular22_1はx1非零chartでq1がx0を消し、q2がirreducible conic。
nil32はx2非零chartでx1,x0を順に消す。singularM1はx0非零chartでx1,x2を消す。
各chart補集合のsurface内dimensionは一以下なので、CIのunmixed性とgeneric reduced性からprime。
cone_nodeのbase Cはx1非零chartでq1がx0を消しq2がirreducible conic、境界は一点。
cone_cuspはx2非零chartでx1,x0を消し、境界は一点。
coneにvertex自由座標を加えてもintegralityは保たれる。

## 2. Eulerを失わない全tangent/conormal線形条件

全150 raw homogeneous quadratic係数を使う。
各q=q1,q2について全五component Π#dqを取り、

`E∧(Π#dq) mod (q1,q2)=0`

を全十pairへ課す。これはprojective tangencyそのものである。
degree四のIはq1,q2に全十五quadraticを掛けたrank二十九の空間であり、
有限QQ rowspace reductionでexactに剰余を取れる。GBを使わない。

chart x0非零でintrinsic normal generators q1/x0²,q2/x0²のbracket numeratorは

`K=x0{q1,q2}−2q2{q1,x0}−2q1{x0,q2}`。

abelian conormalの必要条件はK∈I²である。
prime complete intersectionのI²はI-primaryなので、generic chartだけから全homogeneous条件へ戻る。
degree五のI²はq1²,q1q2,q2²に全五linearを掛けたrank十五の空間。
従ってこの条件もexact有限rowspaceで課せる。
Euler24方向が全条件のkernelに入ることを直接確認し、explicit coordinate minorによりquotientする。
完全kernelと全free coordinate逆写像を保存し、conductorのpoleを仮定してrestrictionを切らない。

## 3. 最弱のG上界を読む方法と範囲

全projective sourceのoriginal350Jacobiとactual126Pfを復元する。
必要なら元quadricsのcoefficient rowspaceからlex leading monomialsを取り、
互いにsupportが交わらないmonomialsの個数だけheight下界を得る。
これはGroebner basisと主張せず、元idealのinitial idealがそのmonomialsを含むという必要方向だけを使う。
全sourceがPoissonでなければ一つ非零元式によるproper閉部分のdimension上界も使える。

固定source dimension upper dからPGL像は24+d−1−dimAutS以下。
実Aut下界はregular22_1三、nil32四、singularM1四、cone_node/cusp五である。
regular22_1はskew centralizer二と二固有値pencil変換の明示lift一を持つ。
具体的にM=N2⊕N2⊕1、H=diag(1/2,−1/2,1/2,−1/2,0)とすると、
W=M/2−Hはδq1=q2,δq2=q2を満たし、両qを固定するcentralizer二方向から独立。
nil32には既受理 `REGULAR_PENCIL_CI22_ABELIAN_CONORMAL_NONDOMINANCE.md` §5の四fieldがある。
singularM1はP2のVeronese座標(u²,−s²,−t²,us,ut)を使うrank二center投影の像で、
PGL3のそのcenter stabilizerは四次元。birational normalizationなので作用はfaithful。
cone二型はvertex変換x4→a*x4+Σbi xi（a非零）の全五parameterを明示的に持つ。
これらはsourceから既に消費したnormalizationを重ねて引かない、fixed surfaceのfree frame fiberである。

本入口はwhole五sourceの有限線形→元Jacobi診断であり、全G閉鎖は実測上界が全型で十分な場合だけとする。
不十分な型では全sourceを保存し、GB・高いjet・少数free点の走査へ自動継続しない。
