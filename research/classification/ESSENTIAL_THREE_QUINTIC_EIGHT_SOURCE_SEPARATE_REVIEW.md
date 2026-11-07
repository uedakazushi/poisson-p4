# 八つの具体本質三変数五次：全源・尖点型の独立review

**独立PASS。** 本稿は八つの**具体**curveとその全GL2 frame・全六磁気源を監査する。
必要な全curveがこの八つに尽きるという、さらに強い純証明は別receiptである。
その全リストの完全性を本稿の仮定にせず、各具体源の成分採録に必要な
複素絶対既約性、essential rank三、全Saito source、内在的尖点型の相違を独立に確認した。
全周囲・真cut・支持及びprior104監査はrootの別作者ambient reviewに委ねる。
global catalogue、原稿、Gitは変更していない。

## 1. 複素絶対既約性と具体正規化

公開inputの各三degree5 binary sectionについて、s/t chartの全gcdが一、
t=0でも三sectionが同時零でないことをQQのEuclid法で確認した。
従って全P1上basepoint-freeなmorphismである。

全degree四十五係数の代入行列21*15、全degree五二十一係数の行列26*21を
独立構成し、全行と保存minorを照合した。rankはそれぞれ十五、二十。
五次の完全kernel一とprimitive Fへの直接代入零も確認した。
これらはQ上の因数分解検査ではない。
P1の複素imageは既約であり、FがC上で可約ならimageを含むproper factorの次数は四以下。
そのfactorにmonomialを掛ければdegree四relationを作れ、rank十五に反する。
従ってFはC上絶対既約である。
basepoint-freeなpullback O(1)の次数五とimage degree五から、mapはbirational、
従ってこのP1は実際の正規化である。
各Fの三偏微分係数行列の全rank三minorも再計算した。

## 2. 全十八syzygyと全十parameter Poisson源

各Fの全gradient quadratic syzygy map28*18を組み直し、rank十六の保存minor、
二本の独立kernel U,V、全row substitution、det(E,U,V)=cF非零を確認した。
したがって既受理の自足的Saito criterionから全AR(F)の二次自由基底となり、
選んだ部分源だけではない。

| curve | det(E,U,V)/F |
|---|---:|
| (4,5) | -9/64 |
| A12 | 225/104 |
| E8+A4 | 2/9 |
| E6+A6 | -175/144 |
| E6+A4+A2 | -125/48 |
| A8+A4 | 125/72 |
| 3A4 | -5/8388608 |
| A6+3A2 | -75/3584 |

全frame四係数と全magnetic quadratic R六係数を残した
Pi=U' wedge du+V' wedge dv+(-div(V')u+div(U')v+R)du wedge dvを構成し、
全十strict Jacobi、全五divergence、全五Casimir contraction、
全五raw Pfaffianと実D=-c det(frame)Fを恒等式として照合した。
局所化はdet(frame)非零、すなわちrank四のための必要条件だけであり、Rに条件は課さない。

各固定Fでdelta(l,m)=div(lV-mU)の全6*6行列と逆も独立に作った。
全行で非零定数determinantを得たため、任意Rに対し(l,m)=-delta^-1 Rを取れる。
actual shear u'=u+l、v'=v+mでRはR+delta(l,m)へ変わるので、これは全Rの実除去である。
frame U'=aU+bV、V'=cU+dVでは、domainを
(l,m)->(dl-bm,-cl+am)へ移すのでdet(delta_frame)=det(delta)*(ad-bc)^3。
従って一般frameだけでなくGL2の全点でinverseが存在する。
vertex GL2による全frame除去も既受理のuniform lemmaを適用する。
accepted canonical essential3 theoremにより、各固定Fの全projective Poisson源は
この表示へ入り、一つのactual PGL5軌道二十三次元となる。

## 3. 全リスト分類を使わない内在的尖点型の認識

homogeneous binary differential cross product phi_s wedge phi_tの全gcdを
二つのchartで検査し、無限遠の零も保持した。
basepoint-free birational mapでは、unit affine chartとEuler関係から、このgcdの
根重複度はbranch multiplicityから一を引いた値に等しい。
gcdとその微分の反復Euclidで、複素全根の重複度を因数分解なしで確認した。

各指定点でunit chartの有理関数を十五次までFractionで展開した。
重複度二ではorder二の座標aを選び、もう一方の線形接線成分とaの全偶数冪を順に除去し、
最初の非零奇数次数を確認した。これは固定degree五curveの型の有限検算であり、
Poisson変形の高次jet探索ではない。

| curve | 独立局所確認 |
|---|---|
| (4,5) | 無限遠でorders(4,5)、delta六 |
| A12 | s=0で最初の奇数十三、係数一 |
| E8+A4 | s=0でorders(3,5)、無限遠で奇数五、係数一 |
| E6+A6 | s=0でorders(3,4)、無限遠で奇数七、係数一 |
| E6+A4+A2 | s=0でorders(3,4)、無限遠で奇数五・係数四、s/t=-1/2で奇数三・係数1024 |
| A8+A4 | s=0で奇数九・係数八、無限遠で奇数五・係数-8 |
| 3A4 | s=0で奇数五・係数-6、他の二根でも奇数五非零 |
| A6+3A2 | 無限遠で奇数七・係数一、残る三つのsimple differential rootsはmultiplicity二 |

3A4の有理representativeは、全必要curve分類から同じ型と推測したのではない。
他の二根はrho=-2+(6/5)sqrt5とその共役である。
QQ(sqrt5)での奇数五の係数は
-47271875/196608-(7046875/65536)sqrt5非零。
このfieldの二つのembeddingとも非零を保つので、三つすべてがA4である。

各行の局所delta下界の和は六、すなわち既約quinticの算術種数を使い切る。
A6+3A2の残る三点はdelta各一以上であり、和から各々正確にA2となる。
このbudgetは別の特異点、複数正規化点のimage衝突、追加の枝も排除する。
従って上表は単にramificationの個数でなく、全curveの正確な内在的configurationである。
全八signatureは相異なる。
本質三変数の最小変数空間は一意なので、PGL5による五次coneの同値は
quotient P2のcurveのPGL3同値を誘導する。この相違はambientでも保存される。

## 4. Receiptと独立範囲

checker `code/classification_eight_essential_quintic_source_review.py`、SHA256
`2e5ba640b75f31c15f2f939709b26d1752146a607dba50cd55a50472adaef8a2`。
job `research/classification/jobs/essential_three_quintic_eight_source_review.json`、SHA256
`33a8991c7de57c24e4b262aa3877e769771d7645494c9f4bfb67156fe138b92d`。
gate PASS後、全八source/cuspは0.202秒、最大RSS14520KiBでPASS。
report `reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_EIGHT_SOURCE_REVIEW_20261107.json`、SHA256
`f35ee50f85b3f07c04652e537abb95f7362f077c74221c7a8cae44afcf72a6dc`。
producer sourceはimmutable source reportのexecuted bytesへ固定し、mutable runnerをimportしない。
標準Fraction、独立sparse helperだけを使用し、CAS・producer helperをimportしていない。

最初のchecker実行はEuclid入力の末尾零をtrimし忘れて停止した。
producerへの異議や数学的失敗ではなく、保存reportも作成されなかった。
入力tupleをtrimする実装修正後に再gateし、上記実行で全八を再確認した。
必要八curveリストの完全性pureproof、全周囲componenthood/prior監査及び新countは別範囲のままである。
