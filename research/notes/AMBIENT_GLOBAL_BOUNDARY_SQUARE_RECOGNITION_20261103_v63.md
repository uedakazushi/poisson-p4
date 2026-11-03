# 全域 N の階数二境界と CB endpoint の二次式が平方となる条件

2026年11月3日。これは全射影既知成分の境界比較であり、全Poisson分類の主張ではない。
全CA∩CBのgeneric-rank≤2部分を尽くしたとは、まだ主張しない。

## 定理

NをS61の18次元交叉成分の全射影27次元閉包とする。
全CBの二つのrank≤2 endpoint像Ba(a=0),Bb(b=0)において、
binary quadratic Fが平方となるproper閉像をBa,sq,Bb,sqと書く。

Nの全generic-Poisson-rank≤2境界は二つの閉像の合併であり、正確に

    N∩{rank≤2}=Ba,sq ∪ Bb,sq

と書ける。Ba,sqは既知CAのf=0境界で25次元である。
Bb,sqはNのg=0境界である。
この包含・認識は他の境界交叉点の不存在を証明するものではない。

## 1. Nの全proper旗模型

全CAの固定旗20次元空間にf∈hW、すなわちx²係数零を課すと19次元線形空間となる。
旗parabolicはhlineとWを保つので、この条件は全群作用で不変である。
既存の全20×20 Lie作用行列でも、x²係数の行がその列以外全て零となることから直接確認できる。
従ってrank19随伴束のP18から全射影空間へのproper射が定まる。

rank4開ではこの条件はS61のrA=0と同じで、全共通旗交叉の九次元q族を与える。
よってこのproper閉像は既知N27を含み、源の次元も27なので、全像はN27である。
実退化五次式はg*h³*fであるから、N内でrank≤2となるsourceの被約台はg=0またはf=0で尽くされる。
その二つのproper像が全rank≤2境界を覆う。

## 2. f=0境界はBaの平方像

CAのf=0表示は

    P=g*h²*∂zw+∂y∧X_H3.

CBのa=0表示は

    P=F(h,x)*∂zw+∂y∧X_H3.

前者は明らかにFが平方である後者に入る。
逆にF=g*ell²なら、a=0模型はhlineを忘れてWだけを使うため、hlineをellへ選び直せばCAのf=0表示となる。
F=0や衝突点も両方のproper閉像に含まれる。
従ってこの二つの閉像は等しい。CAの既存全境界定理により次元25である。

## 3. g=0境界からBb平方像への実座標変換

CAのg=0境界で

    f=h*(alpha*h+beta*x),    beta≠0

とする。H3は全十六係数のまま動かす。
座標を

    h'=alpha*h+beta*x,    x'=h,
    y'=y, z'=z, w'=w

へ変える。元の双ベクトルはP=∂y∧Yと書け、

    Y=−f*∂x+X_H3.

変換後はY(h')=−beta*x'*h'、Y(x')=0である。
projective Euler-verticalを除くためYをY+beta*x'*Eulへ変えると、h'成分は零となり、

    Y(x')=beta*x'²,
    (Y_z,Y_w)=beta*x'*(z',w')+X_H3'.

これはCBのb=0表示で

    f_CB=−beta*x'²

としたものに厳密に等しい。
生の差も正確に

    P_CA_transformed−P_CB=Eul∧(beta*x'*∂y)

である。canonical150係数では全て恒等的に零となる。

逆にBbの平方Fがmarkedhとは異なる線の平方である一般点では、その平方根をx'に取り、
逆座標変換でCAのf∈hW表示へ戻る。
両側の平方条件を持つsourceは既約で、そのproper像は閉である。
従って一般点の双方向の一致から、beta=0、F=0、平方根がmarkedhへ重なる点も含めて全閉像が一致する。

## 証明書

* `code/ambient_global_boundary_square_recognition_20261103_v63.py`
* `reports/AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63/report.json`

全H3十六係数とalpha,betaを記号のまま使い、25個の生行列成分についてEuler-verticalの全差、
canonical150係数の全差をQQ(alpha,beta)上で検査した。測定0.21秒。
全CA/CBのrank≤2交叉を尽くすには、別旗のlinear-space intersectionをさらに比較する必要がある。
