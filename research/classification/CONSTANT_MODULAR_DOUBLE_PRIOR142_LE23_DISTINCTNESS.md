# constant DOUBLE：凍結 prior142 の全 dimension23以下・rank二成分との相違

2026-11-08、PRIMARY E。新算術無しの純 accepted-theorem 合成。
担当 `/root/plane_collision_review`。**producer相違証明、独立査読待ち。**
current142 の読取専用projectionに含まれる108行を全て扱う。
このprojection後の新しい成分、24以上の親、全成分数の確定は別問題である。

## 1. 入力と次元の矢印

DOUBLEのcanonicalは受理済み全highest-double sourceの一actual orbit代表であり、
既exact full51 action matrixのrank48とEuler rank24よりprojective orbit dimensionは23。
本稿は新minorを計算せず、この既証明書を入力として用いる。
そのactual divisorは

```
D=h²*C, C=x²y+xh²+h³.
```

Cは絶対既約で、唯一のdouble hyperplaneはh。
Cのintrinsic essential-plane quotient上、唯一のsingularityは
p=[h:x:y]=[0:0:1]。y=1でX=x+h²/2とすれば
C=X²+h³*(1−h/4)だから、unibranch A2 cuspである。
Dはnonreducedで、全actual orbitとその閉包はproper closed条件
「あるlinear formの平方がDを割る」に入る。

rank二以下はclosed rank条件なので、generic rank二の全六成分にrank四DOUBLEは入らない。
rank四dimension22の全十八成分にも23次元actual orbitは入らない。

残るdimension23の成分ZにDOUBLEが入るなら、PGL不変性よりorbit closure Obar⊂Z。
両者既約dimension23なので Obar=Z。
この等号を使うので、genericDがreducedである成分の境界にnonreducedDが有り得るという
一般論を誤って否定してはいない。

## 2. 八十一のgeneric reducedD成分

凍結projectionのdimension23/rank四は84行。
そのうち29行はMPS T04–T15,T18–T20,T22–T35であり、
accepted dense SNC-local fullsource/PGL theoremからgenericDはreducedである。
SNC一点だけから未証明の全source dimensionを推論しない。
別の52行はprojectionにgeneric actualD reduced=trueを持つ既受理成分で、
全証明pathと独立review statusも各行に保存した。
それぞれの元定理の一般元記述を再使用し、未読の新分類を加えていない。

これら81成分のいずれも、genericDがreducedである一方、
もしDOUBLEを含めば§1より全成分Z=Obarがmultiple-H closed locus内になる。
矛盾。従って全81成分からDOUBLEは除かれる。
全81idは報告書のtwo listsに列挙し、省略された未分類行は無い。

## 3. 残る三つのnonreducedD成分

### 5H zero-log

R4_ZERO_LOG_5Hは全actualDが五乗linear formというproper closed条件に入る。
DOUBLEのD=h²*C3はその外なので含まれない。

### nil22 image二つの nodal-cubic 成分

残りは

```
R4_NIL22_IMAGE_NODAL_CUBIC_FLEX_23,
R4_NIL22_IMAGE_NODAL_CUBIC_NODE_TANGENT_23.
```

いずれも既受理のactual dense23 orbitを持ち、genericD=h²*nodalC3。
flex側canonical cubicは

```
x³−9hx²−3hy²+108h³.
```

h=1、singular point(x,y)=(6,0)でHessian diag(18,−6)は非退化。
node-tangent側canonical cubicは

```
x³+hxy−h³.
```

y=1、(h,x)=(0,0)でHessianのcross項hxは非退化。
従って両者のplane cubic singularityは二枝ordinary node。
DOUBLEは一枝cuspなのでactual divisorのPGL同値ではない。

本質的三変数のplane quotientはgradient-annihilatorから内在的に定まり、
任意のPGL5同値はそのP2 quotientとunique double-H、prime cubicを保つ。
従って余分な二つのvertex方向を混ぜてnode/cuspを同値にすることは出来ない。

仮にDOUBLEがどちらかの23成分へ入れば、§1でそのorbit closureが成分全体となる。
同じ既約varietyの二つのdense locally closed group orbitはともにopenなので交わり、同じorbit。
これは上のnode/cusp相違に反する。
従って両whole23成分にもDOUBLEは入らない。

## 4. 範囲

これで凍結prior142の全rank二六、rank四dim22十八、rank四dim23八十四、計108成分を排除した。
同じ凍結priorの24以上は別proofで扱い、未証明の新global上界を本稿の仮定に加えない。
既ambient D101はscheme tangent24を与えるが、それだけから新23成分又は非被約厚みを推論しない。
他workerのT37親認識が陽性ならこの結果は新成分を意味せず、その境界orbitのprior23非所属となる。
