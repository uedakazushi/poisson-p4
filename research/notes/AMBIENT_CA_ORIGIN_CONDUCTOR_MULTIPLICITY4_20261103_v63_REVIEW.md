# CA原点conductor multiplicity4・全initial・depthの独立検証

2026-11-03。`AMBIENT_CA_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63` を**CAだけ**の定理として承認する。
CBの新しいself-pair packet、CA/CBのactual交叉D、full Poisson schemeは対象外。
既承認の次元上界とactual self-pair geometryの全入力hashを照合し、新しい論理的帰結を監査する。
producer再実行、追加高次jet、GBは行わない。

## intrinsic conductorのgluing

generic self-pair点で全normalization fiberは正確に二つの被約点、各sourceはsmooth29かつimmersion。
受理済みの全image表示を `A=T/(I1 intersect I2)`、`Bi=T/Ii`、`D=T/(I1+I2)` とする。
すると `A=B1 times_D B2`、normalizationは `B1 directsum B2`。
conductorの任意元を `(b1,b2)` とし、normalizerの `(1,0)` と `(0,1)` を掛ければ、
conductorは必ず `ker(B1->D) directsum ker(B2->D)` に入る。逆包含も直ちに成立する。
従って `A/conductor=D` はcanonicalな同定である。旧R29と別chartのalignmentは不要。
二つのbranch接空間の和37をwhole imageのedimと同一視していない。

独立コードには `y*(y-x^m)=0` の二つのsmooth branchという小例を含める。
conductor商は `k[x]/x^m` で、m=1は被約、m=2は非被約。
これも「二つの旗」だけでは非被約性が従わないことを明確にする。

## generic非被約性と第三のsupport branch

actual G×a族の像Zは既承認のrank20によるirreducible20次元族。
上界定理でintrinsic conductor Cはdim20なのでZはsupport componentとなる。
もしZのgeneric conductor長が1なら、他成分・nilpotent supportを除いた非空openで
Cは被約、さらにchar0 finite-typeのregular openでsmooth20となる。
稠密constructibleなG×a像はこのopenへ入り、その点のedimは20となるはずである。
しかし全generic aおよびG移動におけるactual meeting edimは21。
従ってZ上のgeneric長は2以上である。

W/V critical familyは全source critical locusに入り、quadratic表示ではt=0、従ってnu=0。
他方、原点で合流する二つのself-pair source arcはR値が同じで異なるt値を持つ。
tはrank2 algebra全体を生成するので、t値が同じなら全source座標値が同じとなってしまう。
従ってformal arc上で `nu=(t1-t2)^2/4 != 0`。
このt/nuは必要ならformal algebraic recipeをarcへ引き戻して解釈すればよく、
未収束形式級数へ任意の複素数を代入する議論ではない。
Zはcritical W/Vの二branchとは異なる。

## associativityと全initial equality

代数的local conductorでmultiplicity associativityを適用すれば、W/Vの寄与が各1、
Zはgeneric長2以上・support multiplicity1以上なので、e(C)≥1+1+2=4。
completionはmultiplicityを変えない。Zがanalyticに分岐しても、その寄与の総和は減らない。
既受理のe(C)≤4からe(C)=4である。

minimal表示の初期idealに入る `I0=(WU,WV,U*(aU²+bUV+cV²))`、a≠0は、
Hilbert–BurchでCM20次元・e4である。全射 `M=S_gr/I0 -> gr(C)` のkernel Kが非零なら、
`Ass(K) subset Ass(M)`。MはunmixedなのでKはdim20、正のtop multiplicityを持つ。
同次multiplicityの加法性は `4=e(K)+4` となり矛盾する。
従ってK=0、全initial equalityが成立する。未知のb,cやhigher termを零に置かない。
Hilbert seriesは `(1+T)^2/(1-T)^20`、CはCM20である。

等号後のassociativityはgeneric長を正確に(1,1,2)へ固定し、他の20次元branchを排除する。
CMのunmixed性により低次元minimal/embedded primeもない。Zのsupport multiplicityは1。
excellentな被約equidimensional local supportのmultiplicity1から、その完成もregularかつ一枝。
従って原点のformal minimal primeは正確にW,V,Zの三つとなる。

## 全十生成元・type・actual image

初期idealは二つの二次式と一つの三次式で生成されるため、complete filtered divisionで
対応する三つのactual liftが残余idealを生成する。U³係数が非零なので三次式は
`m*(WU,WV)` に入らず、この三生成はminimal。
七linear functionを戻してmu(I)=10。
独立コードは全150係数上で三liftと七linear functionからraw十rowへの遷移を作り、
その10×10行列式1と外側140列の零を確認する。
raw rowは18,19,28,29,60,67,72,46,95,63である。

type2はgraded typeの上界だけから断定しない。actual残余は正則22変数上のheight2 CM idealで
minimal生成数3なので、そのminimal Hilbert–Burch resolutionの最後のrankは2となる。
これがactual CM type2を与える。

R29上の `0->I->R->C->0` からdepth I=21。
`A=R directsum I*t` よりdepth A=21。
`m_A=m_R directsum I*t`、`I subset m_R`、`nu subset m_R` なので、
odd product `nu*I² subset m_R³`。従ってedim A=29+mu(I)=39であり、
`m_A²=m_R*m_A` も成り立つ。m_Rがm_Aのreductionであることを確認した上で、
finite R-rank2からintrinsic multiplicity e(A)=2を得る。rankだけの推論ではない。
coneはdim29/edim39/depth21/e2、nonzero coneのscalar smooth factorを除くprojectiveは
dim28/edim38/depth20/e2。Aが一つのregular normalizerへ注入するdomainであることは変わらない。

一般Z点ではnilradicalの平方零性とrank1局所自由性、smooth20 reduced supportをopen上で取れる。
formal smoothnessでsplitするため、conductor完成は `C[[g20,eps]]/(eps²)`。
これはgeneric self-pairでの記述で、原点の未知higher analytic normal formは同定していない。

独立コードとhash証拠は
`code/ambient_CA_origin_conductor_multiplicity4_review_20261103_v63.py` と
`reports/AMBIENT_CA_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/` に保存する。
