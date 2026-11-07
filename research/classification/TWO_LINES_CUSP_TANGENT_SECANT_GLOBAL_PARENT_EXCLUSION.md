# 尖点接線・割線の全源：二楕円成分とT21の親を除外する

2026-11-07、PRIMARY `E_reduced_nonnormal`。producerの全係数計算と下記の手証明を保存する。
独立review前であり、新成分・全分類完了・他の親の不存在は主張しない。

## 1. 全六磁気源からの忠実な一parameter判定

対象と全translationの逆は
`TWO_LINES_CUSP_TANGENT_SECANT_MAGNETIC_NORMAL_FORM.md` に固定した。
I=r5≠0の全六源は、可逆PGLと非零scalarにより
`R=3(lambda h²+y²)`、lambda=(r0+r2)/r5の全源へ移る。
固定評価点(1,1,1,0,0)のPfは全lambdaで−9、Darboux行列式は−1/9である。
faithful S61へ移る全PGLと逆、Euler差、全150係数の逆再構成、全134式、
元の全350射影JacobiをQQ[lambda]上で確認した。
parameter分母は非零有理定数だけであり、有限lambdaを除外していない。

T39全45graph残式のうち11本、T40では13本が非零である。
前者は定数1/4、後者は定数−2/3を含むから、両方の全残式idealは厳密に(1)である。
従って特殊lambdaの例外も無く、I≠0の全対象点はC39とC40の両方の外にある。
全global componentと忠実S61全graphの同値は、既受理の
`S61_GLOBAL_COMPONENT_EVALUATION_FIBER.md` とその独立reviewを用いた。
この陰性をI=0へ延長することはしていない。

公開reportは
`reports/CLASSIFICATION_TWO_LINES_CUSP_TANGENT_SECANT_WHOLE_GRAPHS_20261107.json`。
有限job `jobs/two_lines_cusp_tangent_secant_whole_graphs.json` のgate後、2.241秒で完了した。
既two-lines runnerへisolated modeを追加しただけで、旧五関数のASTは不変である。
共有Ga9 runner、旧source、旧証明書は変更していない。

## 2. 二単純超平面に付随する閉条件

一般階数四のprojective Poisson tensorと二つの相異なる単純超平面因子H,Wを考える。
Hを無限遠に選んだquadratic affine代表の発散を、射影linear vector field Z_Hとする。
もう一つの超平面はPoisson divisorなので

```
X=P#dlog(W/H)
```

も全体に正則なlinear projective vector fieldである。
H-Casimir homogeneous代表では単に `X=P#dW/W` となる。
割るWは全括弧を正確に割り、Wの零集合を捨てる操作ではない。
H-Casimir代表を作る線形Euler gauge、発散、Wでの正確な多項式除算は、
独立な二linear formsの各座標chart上でtensorとH,Wの正則な演算である。

次の二条件は閉条件である。

```
(a) Z_H+X は W 上でprojectively零、
(b) Z_H+4X は H 上でprojectively零。
```

projectively零とはEulerとのwedgeがその超平面で零という意味で、
各条件はhomogeneous係数の多項式恒等式である。
HとWを交換すると `X_new=−X`、`Z_W=Z_H+5X`（Euler差を除く）であり、
(a),(b)は交換される。この式は二つのaffine volumeの比(W/H)^(-5)からも得られる。
従って二超平面の順序に依存しない。

これをD=HWC3の因子flagとともに取る。
H,W,C3の射影係数空間の積から五次形式への乗法写像はproperである。
その像の閉包でDが依然として二つの相異なる単純Hと既約三次からなるなら、
極限の二超平面は正確にその二因子であり、既約三次へ紛れ込めない。
従って両条件を一般源が満たせば、この種類の全boundaryでもどちらかの因子順で満たす。
途中でH=Wとなるflag補集合をtargetへ無理に適用していない。

## 3. 完成されたT21全源は両条件を満たす

accepted public sourceは
`reports/MPS_S61_RAW_CATALOGUE_20261103_v62/T21.json` である。
全period p0,p1と全edge t0,t1,t2をそのまま残す。
座標をA=x0,B=x1,C=x2,H=x3,W=x4と書くと、全非零括弧は

```
{A,B}=t0 C², {A,C}=t1 B², {B,C}=p1 BC+t2 A²,
{B,H}=p1 BH/3, {C,H}=−p1 CH/3,
{B,W}=p1 BW/3, {C,W}=−p1 CW/3, {H,W}=p0 HW.
```

Eを全Euler場、L=(0,−p1 B/3,p1 C/3,0,p0 W)とし、PからE∧Lを引くとH行は零になる。
残る三変数括弧は厳密にJac(f)で、

```
f=p1 ABC/3+t0 C³/3−t1 B³/3+t2 A³/3,
P_H=Jac(f)+p0 W∂W∧(A∂A+B∂B+C∂C).
```

これは全八括弧を直接引き算した恒等式であり、Hesse normal formの外部分類は使わない。
全periodとedgeの積を保ち、パラメータを一つも反転していない。
Jac(f)の発散は零だから、lambda=p0として

```
Z_H=lambda E_ABC−3lambda W∂W, X=−lambda E_ABC,
Z_H+X=−3lambda W∂W,
Z_H+4X=−3lambda(E−H∂H).
```

最後の二式は(a),(b)そのものである。
このcompleted sourceの閉包がT21全成分であることは、
`research/notes/MPS_S61_CATALOGUE_20261103_v62.md` の受理済みsource定理への依存である。
従って§2のflag閉性から、二単純H+既約三次を持つT21の全boundaryも両条件を満たす。

## 4. 対象全磁気源は一つの条件を破る

対象のh-Casimir代表から、全六磁気係数について厳密に

```
Z_h=−4h∂z+(3h−12x)∂w,
X=P#dlog(x/h)=h∂z+(−3h+3x)∂w.
```

磁気項はbaseだけの二次式なので、これら二場に寄与しない。
ゆえに

```
Z_h+X=−3h∂z−9x∂w,
Z_h+4X=−9h∂w.
```

x=0上で第一の場は−3h∂zとなり、projectively非零である。
具体的にEulerとのwedgeの(h,z)成分は−3h²である。
第二の条件が成立しても第一の失敗は消えず、H,xの交換でも両条件が入れ替わるだけである。
従って**I=0も含む全六磁気源はT21成分の外**である。

この証明は新しい点・小行列・局所jetを使わず、全実源に一様な二linear fieldsと
accepted full T21 sourceの全parameter恒等式だけを使う。
T36/T37/T38又は既知Ga9 cuspidal24成分への帰属は未解決のまま残す。
