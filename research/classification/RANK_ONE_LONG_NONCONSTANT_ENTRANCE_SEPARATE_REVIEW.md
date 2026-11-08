# rank一 long translation・非定数gの25次元全像：独立の純source監査

2026-11-08、PRIMARY E_multiple_hyperplane。producer
`RANK_ONE_LONG_NONCONSTANT_COMPONENT_ENTRANCE.md` (SHA256 a73cb9c7…)
の **既約25次元全Poisson像とprior129分離までをPASS** とする。
全空間成分性、smoothness、新成分の計数はまだ結論しない。
このnote作成時に新ambient点、operator、rank、minor、kernel又はjetを計算していない。

## 1. 元座標の全源・Jacobi・injectivity

X=∂x+x∂y、r=y-x²/2、q=(z,w) とする。gは任意の非定数 affine 式、
Uは affine 場で d=divU、U(g)=dg、Kは任意の次数二以下の多項式。
F=(U+d+1)K、H=gK とすると、元の四affine座標で全六括弧は

```
Pxy=-F-(1+d)r, Pxz=u, Pxw=v,
Pyz=xu+r*g_w+(gK)_w,
Pyw=xv-r*g_z-(gK)_z, Pzw=g.
```

全定数・一次・二次係数を保持する。補助変更(t,r,q)=(x,y-x²/2,q) は
Jacobian一の多項式同型であってPGLではないので、成分membershipには使わない。
この座標では P=∂t∧Y+Jac(f)、f=g(r+K)、
Y=-(F+(1+d)r)∂r+U、divY=-1、Y(f)=-f。
従って L_Y Jac(f)=Jac(Yf)-(divY)Jac(f)=0 で全Jacobiが成立する。
元Pfは -g(r+K)、column divergenceはX。任意のKでfのr係数gは非零なので全源でrank四である。

完全性は旧 `MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §5 の全逆表示(19)を使う。
g=zのchartでU=alpha*z∂z+(beta+gamma*z)∂w。
(U+1)H=Fgをz0に制限すると(1+beta∂w)H(0,w)=0であり、
有限次数の多項式上この作用素は可逆なのでH=zK。
ここで積分定数を勝手に消していない。

gはPzw、UはPxz/Pxw、FはPxyから一意に復元される。
Pyz/PywはHの両偏微分を決めるので、二候補の差H-H'は定数c。
(U+1)(H-H')=0はc=0を与える。従ってg非零からKも一意である。
Kの加法定数に隠れた正次元fiberはない。

## 2. 固定Xのsource Sの既約次元12

gの基底はaffine三係数空間から線形部分零の直線を除いた既約三次元開集合。
任意のgはqの可逆affine変換でzにでき、上のUの全解は正確に三次元。
従ってU(g)=dgの核はこの基底上のrank三vector bundleである。
KはA6全体。Sは既約12次元であり、injectivityから実際の固定X像も12次元。
例としてg=z、任意U、K=z²+w²では後述quadricはrank五。

## 3. StabXと全sourceの実際のfiber

斉次化N=h∂x+x∂yのGL5 centralizerは、順序(h,x,y,z,w)で

```
[a,0,0,0,0]
[b,a,0,0,0]
[c,b,a,d,e]
[f,0,0,g,i]
[j,0,0,k,l]
```

という十一係数の全行列、det=a³(gl-ik)である。
Nの二乗の双対像はhの線だからhを保つ条件は既にこの形に含まれる。
Aff4でh-row=hに正規化するとa=1だけを課すので、StabXは正確に十次元。
Aff4自体は4²+4=20次元である。

全stabilizerは

```
x'=x+b,
y'=y+b*x+ell(q)+c,
q'=Bq+t, detB!=0,
r'=r+ell(q)+c-b²/2.
```

を与える。qold=B^-1(q'-t)、Delta=ell(qold)+c-b²/2とすると、
全source parameterの実写像は

```
g'=det(B)*g(qold), U'=B_*U,
K'=K(qold)-Delta, d'=d.
```

である。U'(g')=d'g'、K'の次数二以下、g'の非定数性を保つ。
Pxyの変換では追加項ell(U)がF'=(U'+d'+1)K'と一致する。
Pyz/Pywは2×2 identity B*J=det(B)*J*B^-tと上のr'/K'で一致する。
従って単にPfの一致から全tensorの一致を推論したのではない。
逆も同じaffine群にあり、SはStabXで全域不変である。

map Aff4×S→Aの同じ像Pを与えるaは、divergenceのequivarianceから
StabXの一cosetに限られる。逆にその全cosetに対しa^-1_*PがS内にあり、
injectivityからsource parameterは一意である。従ってfiberは正確に十次元。
像Aは既約で dimA=20+12-10=22。

さらにa_lambda=(lambda*x,lambda²*y,q)とすると全lambda!=0で

```
lambda*P(g,U,K) = (a_lambda)_*P(lambda*g,U,K/lambda²).
```

これはg、U、r及び全六括弧で直接確認できる全source identity。
単にlambda Xの共役可能性だけでcone性を仮定せず、この式によりAが全非零tensor尺度を保つ。

## 4. marked Hを忘れても25次元

固定hのh-Casimir homogeneous liftは一意なので、affine source像をEuler商へ移す際に
隠れた正次元fiberはない。Aのprojective class像はcone性から21次元。
Hを動かすmarked sourceは PGL5/Aff4 の四次元を加えたassociated bundle像で25次元。
g=z chartでは実Dは

```
D=-h²*z*Q, Q=h*y-x²/2+K2(h,z,w).
```

K=z²+w²でQはrank五（h/yのnondegenerate2×2blockとx/z/w三つの非零square）。
従ってgeneric Qはsmoothで絶対既約、zはhと相異なる。
hが唯一の重複因子、重複度二であるので、marked Hを忘れるgeneric fiberは一点。
全PGL飽和像の既約閉包Zは正確に25次元。
rank四部分のDが重複一次因子を持つ条件はprojective multiplicationのproper像として閉であり、
Z全体のrank四部分にも保たれる。

## 5. 固定prior129の全比較

projection `rank_one_long_nonconstant_prior129_20261108.json` を、記載された実Gitcommit
`e2e718783b6811414a7f3c4b128a04a4503e5896` のCOMPONENTS blobと独立に照合した。
全catalogue SHA256は77e77007ce29dba6e9609c3e7f46bfc720c7139f69844c67d229e6490d13007eで一致し、
全129行のid/name/rank/dimension/genericD projectionが一致した。
rank四の内訳は84×23、18×22、12×24、7×25、2×28で計123、rank二が六。

* rank二六成分は閉rank条件によりZを含まない。
* rank四dim<=24の114成分はdimZ25により含まない。
* CA/CBは全Dが五一次因子積のproper閉条件内。旧普遍源では任意magnetic cubicを保持しても
  Pfが3H+H1+H2又はH+2H1+2H2であり、この閉条件は全成分へ延長する。
* nil31-DDH/half/third25も全Dが五一次因子積のproper閉条件内。
  generic rank五Qを持つZはその外である。
* positive NEW25は全DがL³*Qのproper閉条件内。Zのgeneric Dの最大重複度二と両立しない。
* 残るT01/T39/T40はdim25。もしZを含めば既約性と同次元性からZと等しい。
  しかしZは全rank四部分で反復Hを持ち、既受理MPS三成分は被約Dの非空開を持つので等しくない。

これで全prior129からの分離が閉因子条件又は次元により成立する。
未知の26次元以上の親の不存在を結論したわけではない。

## 6. 受理範囲と次の最弱判定

今回受理するのは既約25次元全Poisson閉族Zの存在と全prior129からの分離。
Zが全空間成分かは未証明である。
一つの構造的rank四・Qrank五anchorで元350×150Jacobi微分rank100が得られれば、
Euler24とscalar一を除いた全接上界25と実Z下界25が一致し、smooth/unique25と全成分性が従う。
rankが100未満なら新点列へ移らず、その値と全arraysを保存する必要がある。

g定数、g零、shortchainの全境界はこの新Zへの包含をまだ主張しない。
旧nonlinear r座標をPGL同値とせず、旧全source/normalizationと現在の元affine sourceだけを用いた。
新job、CAS又はambient rank計算は本純reviewでは行っていない。
