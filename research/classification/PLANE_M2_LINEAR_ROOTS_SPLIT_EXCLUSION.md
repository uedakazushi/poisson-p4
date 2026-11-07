# m2 flag(1,1)・split determinantは全空

2026-11-07。主義務G/E。新しい全source数学証明、計算jobなし。
rank二restriction plane、reduced nonzero cubic R=(f=0)、実D被約既約quintic、位数二、
normal leading quadratic generic squarefreeとする。
全root-matrix flag(1,1)を

\[
 G_2=h\ell m,\quad \binom\ell m=M\binom uv,\quad h\det M=f
\]

と書き、detMが相異二直線の積である全枝を扱う。
hはそれらと異なるline。従ってRはnonconcurrent又はconcurrent triangle。
**結論：この全Poisson sourceは空である。**
conic determinantは別noteで全MPSへ認識済みであり、本ノートと合わせflag(1,1)全体が完了する。
他flags、double normal quadratic、nonredR、m3等は別残余である。

## 1. root-line connectionと正当なraw gauge

`PLANE_M2_ROOT_LINE_CLASS_OBSTRUCTION.md` の全証明を使う。
二primitive root rowsのdegreesが同じ一、hが一つのlineなのでη=θ1−θ2=0。
従ってΘ=θIで、全conormal connectionは

\[
 \beta=\theta I-M^{-1}dM. \tag{1}
\]

θはclosed logarithmic formだが、βのnonlog部分を除いていない。
非被約Rにはこのη零の論証を適用しない。

canonical tracefree mixed gaugeを取る。
off-diagonal entry βvu=0なら、対応するhorizontal formも零、
iE iVvol=3fτのtracefree linear-field同型から対応mixed columnも零となる。
このgaugeではbase partはJacobian P_f、normal-normal bracketは全三係数を持つ
K=k20u²+k11uv+k02v²である。
Euler representativeを途中でstabilizerと誤認せず、projective Euler-Pfaffianの不変性だけを用いる。

## 2. 一つのnormal coefficientで十分なPfaffian obstruction

βvu=0であれば、全base-normal v-columnはv倍であり、v=0で零となる。
actual Euler-Pfaffianをv=0へ制限すると、base三×normal-normalの項だけが残り、

\[
 G_2|_{v=0}
 =\bigl(x_0f_{x_0}+x_1f_{x_1}+x_2f_{x_2}\bigr)k_{20}u^2
 =3f\,k_{20}u^2. \tag{2}
\]

mixed u-columnにv項があっても、mixed v-columnが零なので寄与しない。
normal coordinatesをEuler contractionする項も、u側はv-column零、v側はv=0で消える。
k20はconstant。従ってnormal u² coefficientは必ずfで割れる。
これは全higherP1/P2以前の必要条件であり、generic rank/minorや追加Jacobiを使わない。

## 3. linear image三次元の全matrix

W=im(M)⊂Mat2が三次元ならhyperplaneのannihilatorはrank一。
rank二ならdetのrestrictionはsmooth conicとなり、本split仮定に反するからである。
annihilatorをrow vector λとcolumn vector μの積で書く。
primitive両row条件からλの両成分は非零。片方零なら一rowがconstant normal directionとなり、
degree一の共通因子を持つためprimitiveでなくなる。
legitimate row rescalingとnormal右GL2により、両rowのfirst columnを同じnonzero base form xへ戻せる。
base coordinatesを取ると、全marked matrixは

\[
 M=\begin{pmatrix}x&y\\x&z\end{pmatrix},\quad
 c=z-y,\quad\det M=xc.
\]

marked rootsを混ぜるleft GL2をsourceのparameter無しで使っていない。
direct differentiationから

\[
 M^{-1}dM=
 \begin{pmatrix}
 d\log x & (z\,dy-y\,dz)/(xc)\\
 0&d\log c
 \end{pmatrix}.
\]

式(1)はβvu=0を強制する。
一方target G2=h(xu+yv)(xu+zv)のu² coefficientはh x²。
f=hxcでありx,cはindependentなので、h x²はfで割れない。
式(2)に反し、全sourceが空となる。
nonlog off-diagonal termを保持したまま、一つのPfaffian coefficientだけで閉じた。

## 4. linear image二次元の全matrix

W二次元でdetの二rootが相異なるregular matrix pencilは、strict equivalenceで

\[
 M=U\operatorname{diag}(x,y)V,\qquad U,V\in GL_2(\mathbb C)
\]

と書ける。
normal右GL2でVだけを除き、**left Uは全parameterとして保持する**。
各rowがprimitiveなのでUの全四entriesは非零。
一entry零ならそのrowの残成分にx又はyのdegree一共通因子があるからである。

Θscalarなので、Uを除かずに式(1)を計算しても

\[
 \beta=\theta I-\operatorname{diag}(d\log x,d\log y).
\]

したがってβはdiagonal、特にβvu=0。
target G2=h(U11 xu+U12 yv)(U21 xu+U22 yv)のu² coefficientは
h U11U21 x²で非零。一方f=h det(U)xyなので、これもfで割れない。
式(2)に矛盾し、この全U-familyも空となる。

## 5. 完了範囲

matrix image一は両linear primitive rowsに不可能で、detのsplit仮定にも反する。
従って以上でsplit determinantの全flag(1,1)が空となる。
codimension二のroot-row basepointを除外したり、special Uをsamplingした証明ではない。
sourceのrow共通因子がある場合は最初のprimitive factorizationで別flagへ移り、ここへ無断で含めない。

`PLANE_M2_LINEAR_ROOTS_CONIC_MPS.md` と合わせると、
reducedR、m2、generic squarefree、flag(1,1)の全sourceはempty又はMPS40に被覆される。
flags(0,1)/(0,2)、nontriangle(0,0)、normal double quadratic、nonredR、m3は未完了。
新しいsource/minor/GB/jetは一切計算していない。
正式独立receiptは別ノートへ保存する。
