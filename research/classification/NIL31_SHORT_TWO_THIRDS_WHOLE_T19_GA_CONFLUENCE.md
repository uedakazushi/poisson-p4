# nil31 short r=2/3全源：T19 quadric-pairのmoving-H Ga曲線

2026-11-08、PRIMARY E、producer integrator `/root`。**純手導出、独立査読・有限readback前。** 新operator/point/jet/profileを実行しない。旧二次/三次局所routeは停止したまま、全実PGL曲線によって別の陽性橋を構成する。旧notes/code/reports不変、新成分増分零。

## 1. 目標の全sourceと必要な範囲

`MODULAR_NIL31_SHORT_SUPPORT.md` §7の全normalizer逆を使う。斉次座標(n,x,y,z,w)、n-Casimir liftで、全六は
Pxy=d*n*z-t*n*y+2t*z*w-y²/3,Pxz=t*n*z-2y*z/3,Pxw=d*n²+y*w/3,
Pyz=-2z²/3,Pyw=(n*y+z*w)/3,Pzw=2n*z/3。              (1)
全rank四はtunit、d任意、signed actualD=-t*n²*z*(n*y-z*w)。t/d両unitは既actualtorusで一軌道、d0も同じwholepolynomialsourceの閉境界。本稿の曲線はt=1/27,d=-sigma/108を全sigmaに対して直接与える。従って全tunit,dを元torus正逆で戻せる。特にd0を未証明の一般点から自動的に含めたとはしない。

## 2. 真T19 sourceとactual frame

真completedT19座標(a,b,c,d0,H)、period(r,v)、edges(t0,u0)。sourceのd0と目標dを混同しない。
t0=(r-v)/2,u0=r-vによりcoupledPac=2t0u0*bH/(r-v)=(r-v)bH。
H-Casimir helper全六とrawEuler差は `TWO_QUADRICS_NIL22_LAST_CORE_WHOLE_T19_CONFLUENCE.md` (3)。二quadric q3=c*d0-b²,q4=a*b+d0*H、actualD=-r(r-v)Hq3q4。

曲線parameter e、q=e²、全磁気sigmaに対し
r=e², v=3e²+e^4*sigma。                            (2)
period開-r(r-v)=e^4*(2+e²sigma)unitは全sigmaで非空。
q3/q4をliteral保存するactualGa flowのlambda=1/e,mu=0を使う。新中間座標(h,X,Y,Z,W)から**実oldsource座標**への全inverse frameを直接書く：
a=eX-Z/e,
b=eW-h/e,
c=Z-2W+h/e²,
d0=h,
H=Y+X-Z/e²。                                    (3)
これは先にa_new=eX-Z/e,b_new=eW,c_new=Z,d_new=h,H_new=Yとし、真Gaのinverseを合成したもの。source-from-target det=e²。非線形の変数変更をactualframeと呼ばない。

全forward逆は
h=d0,
W=b/e+h/e²,
Z=c+2b/e+h/e²,
X=a/e+Z/e²,
Y=H-X+Z/e²。                                    (4)
全eunitで両合成がidentity、rawEulerにも同じframeを用いる。
二quadricと親hyperplaneは正確に
q3=hZ-e²W²,
q4=hY-ZW+e²XW,
H=Y+X-Z/e²。                                    (5)
従って親Hは中心でZへ動き、rank三quadricはhZへ分裂する。中心の二重HはZであり、元chart hとは異なる。このmoving-Hと二つの非一様normal scaleが旧不成立の一様blowupとは違う点である。

## 3. 全chart六係数のpolynomial表示

h一で計算専用(U,V,I,J)=(b_new,-a_new,c_new-b_new²,H_new+a_new*b_new)を使う。
U=eW,V=(Z-qX)/e,I=Z-qW²,J=Y-ZW+qXW。
この補助は非線形だが、globalactualmapは(3)だけである。真T19の全aux括弧は既note(8)。Ga pushはU→U-1/e,V→V（I,J固定）。
K=(2+q*sigma/2)(qW-1)(Z-qX)-q(1+q*sigma/2)Jと置く。
これはpush後の{U,V}そのもの。中間targetの六括弧を次の**全polynomial**で与える：

PZW=(-2-q*sigma)(qW-1)I,
PXW=(qW-1)*[-sigma*Z/2+(2+q*sigma)W²-(2+q*sigma/2)X]-(1+q*sigma/2)J,
PYW=-(3+q*sigma)(qW-1)J-W*K,
PXZ=(4+q*sigma)(Z-qX)(Z-W)+q*(-2-q*sigma)*W*J,
PYZ=q*J*[-2Z+6W+q*(-2W²+2sigma*W)-q²*sigma*W²]
     -(Z-qX)*[2Z+q*((2W+sigma)Z-6W²)+q²*(2W³-2sigma*W²)+q³*sigma*W³],
PXY=((Z-qX)*K-PYZ)/q+(3+q*sigma)*(Z-qX)*J。       (6)

最後の表示の分子はq0で(-2Z²)-(-2Z²)=0だから、C[q,sigma,X,Y,Z,W]で厳密にqを割る。これはlocalizationではなくpolynomial quotient、負e冪無しの全identityである。PXWも元(PZW+K)/qの同じ厳密割算を上で展開した。
各全式はe/sigmaの多項式。generic eでは実PGLを掛けた真T19のchart classで、degreeが一見二を越えるchart項も捨てていない。
固定hchart restrictionはglobalquadraticprojectivebivector空間のinjectivelinearmap、imageは閉linear。従って(6)は全eでそのimageに属する。GL-equivariant canonicaldivzero sectionを採れば、全150係数のquadratic global liftはe多項式である。generic h-Casimir又は旧rawEulerpoleの個別零を仮定しない。有限readbackはこのcanonical全150を直接確認する。

## 4. 中心の全六括弧、二重HへのEuler gauge、固定PGL

(6)のq0中心は全てquadraticとなる。中間h一で
PXY=2XZ+4YZ-6WY+sigma*Z²/2,
PXZ=4Z(Z-W),PXW=2X-Y+ZW-2W²+sigma*Z/2,
PYZ=-2Z²,PYW=3Y-ZW,PZW=2Z。                       (7)
一般qの高次chart係数から中心quadraticへの相殺を省略しない。例えばPYZのq一次係数は2XZ-2YZ+6WY-sigma Z²であり、ZW²二項は正確に相殺する。
(7)をdegree二に斉次化したh-Casimir liftを取る。n=ZをCasimirへするため
L=(L_h,L_X,L_Y,L_Z,L_W)=(0,-4Z+4W,2Z,0,2h)
としてP'=P-E∧Lを用いる。新Z行は全零で、残全六(X,Y,h,W)は
P'XY=-2WY+sigma*Z²/2,
P'Xh=-4hZ+4hW,
P'XW=-hY-3ZW+2W²+sigma*hZ/2,
P'Yh=2hZ,P'YW=hY+ZW,P'hW=-2h²。                 (8)

次の固定actual linear mapを掛ける：
n=Z,
x=X/18+7Y/54,
y=-W/3,
z=h,
w=-Y/3。                                         (9)
中間順(h,X,Y,Z,W)から(n,x,y,z,w)へのdet=1/162、全逆は
h=z,Z=n,W=-3y,Y=-3w,X=18x+7w。
(8)からの全六pushは
Pyz=-2z²/3,Pzw=2nz/3,Pyw=(ny+zw)/3,
Pxz=nz/27-2yz/3,
Pxw=-sigma*n²/108+yw/3,
Pxy=-sigma*nz/108-ny/27+2zw/27-y²/3。           (10)
すなわち(1)のt=1/27,d=-sigma/108そのもの。全磁気sigmaはPxy/Pxwに残る。canonical中心は(10)からE∧N/5を引いたliftと一致し、N=y∂x+z∂y+n∂w、rank3/Jordan(3,2)である。

## 5. actualPfとwhole-source結論

(3)のdetinverse=e^-2とsourceDから全intermediateactualDは
D_e=(2+e²sigma)*(e²(X+Y)-Z)*(hZ-e²W²)*(hY-ZW+e²XW)。  (11)
固定(9)のdet1/162を掛け、中心でQ=hY-ZW=3(ny-zw)となるから
D_0=-n²*z*(ny-zw)/27。
これは(10)の全Pfの直接計算と一致し、全126係数の恒等式として保持する。
(6)とcanonicalsectionは全eで正則、中心は非零、generic eは真completed T19。従って閉T19に(10)の全sigmaが入る。元torusで任意tunit,dへ、旧全六parameternormalizer逆で元全sourceへ戻す。

**独立手査読・同一曲線exactreadbackが通れば、nil31 short r=2/3の全源はT19へ帰属する。** このproducer保存時点ではその後続段階を先取りしない。新成分数の増分零、全E/上界/主分類未完。旧二法線診断・cubic零類・固定1PS障害は変更せず、movingtrueGa/periodcurveという別の構造的証明を使う。
