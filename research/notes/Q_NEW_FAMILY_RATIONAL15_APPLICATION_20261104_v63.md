# Actual firstjet像B_NEWはQQ/C上で有理15次元

2026-11-04。これは三つの**FINAL独立レビュー済みgate**を合成する応用である。
対象はW16内のreduced firstjet image `B_NEW`、次元15。quadratic Poisson成分25全体では
ない。関数体・section・generic stabilizerを混同せず、下記三gateのhashを専用receiptで
固定して初めて結論を採用する。

## 1. 使用する三gate

1. TARGET_FIELD_AND_SPECTRAL_DEGREES_REVIEW:
   actual dominant-source pullbackにより
   `E_Q=QQ(B_NEW)^Sp4=QQ(xi,r,s,q,B)`、index1。
   同じ主張はC上でも成立し、QQでの定義体と幾何的base extensionも受理済み。
   別の有理表示 `QQ(xi,R,Z0,r,s)` との間の両向き回収も、この同じ体についての結果。
2. RATIONAL_ORBIT_SECTION_REVIEW:
   上の五座標上のactual rational16-vector `beta_bar(z)`、z=(xi,r,s,q,B)、
   およびそのSp4-saturationが同じB_NEWに稠密であること。
   後者は元六sourceのt=2 sectionでkernel derivationの横断性を検証した結果である。
   「一つのorbitに入るから稠密」という推論ではない。
3. SP4_CENTRALIZER_FREE_REVIEW:
   actual既存source点0/prime32003の32×15 matrixのminor17760が、QQ上の
   非零polynomial minorを与える。非空openで全GL4共通commutantがscalar、
   従ってSp4 stabilizerは±Iに含まれ、affine W16上の−Iはbetaを−betaへ送るので
   generic geometric stabilizerはidentityである。

第三gateはLie stabilizer零だけの証明ではない。trace/4を引くことで**全group**の
commutantをscalarにし、残る有限centerも落としている。char0のgroup-schemeに
対しても同じ結論が成立する。選択minor自体をSp4-invariantとは仮定しない。

## 2. 同じquotient pullback上の有理section

Sを表示された非零条件を満たすA5の有理openとする。
beta_bar:S→B_NEWはQQ上のrational mapである。補助式u²=1600rの平方根を
QQ(S)内の有理関数と仮定していない。元sourceとのorbit同一性は二次補助拡大内の
実Sp4行列で示し、その後beta_bar全16式のu奇数部が消えることを検証して、
vectorそのものをQQ(S)へ下降させている。

Sp4不変有理関数のrestrictionは、その同じorbit恒等式により元sourceのpullbackと
一致する。補助拡大後の恒等式は、もとのfieldへの包含が単射なので下降する。
gate1のindex1により、五section parametersに対応する五つの**actualE_Qの元**を
選べる。これらは有理quotient map `pi:B_NEW ⇢ S` を定義し、

`pi(beta_bar(z))=z`

がrational mapsとして成立する。有限spectral baseだけへの写像や、sourceの別の
五変数に置き換えた推論ではない。Inv14を使う関数体回収のexpanded式は不要だが、
actual field inclusionとindex1は本質的である。

各invariant rational functionのregular locusはSp4安定である。もしsection全体が
そのpole/undefined locusに入れば、稠密なSp4-saturationもその真の閉集合に入って
しまう。従ってSを縮めて、上のrestrictionとpiを同時に定義できる。

## 3. Product mapの支配性とgeneric degree1

QQ上のrational map

`Phi:Sp4×S ⇢ B_NEW, (g,z)↦rho(g)*beta_bar(z)`

を取る。gate2により支配的。sourceとtargetは共に次元15である。

gate3のfree locus Uは非空なSp4安定openである。sectionがUを全て失うなら、
その全Sp4-saturationはB_NEW\Uに含まれ、稠密性に反する。従ってsectionの
generic pointもgeometrically freeである。sectionがgeneric torsorの境界orbitだけを
選んでしまう可能性を、ここで排除している。

幾何的generic imageに対して

`g1*beta_bar(z1)=g2*beta_bar(z2)`

なら、piとそのinvarianceからz1=z2。そのgeneric zではstabilizerがidentityなので
g1=g2。従ってPhiの幾何的generic fiberは一点である。

支配的な同次元rational mapはgenerically finiteであり、char0なのでそのfield extensionは
分離的である。一点fiberをpurely inseparableなdegree>1と取り違えることは無い。
よってdegreeは1、**PhiはQQ上birational**。同じargumentはC上でも成立する。
全写像がQQ上で定義され、幾何的generic degree1であるため、単にC上で得た逆を
QQへ仮定しているのではない。

## 4. Torsorと有理性の結論

有理quotient fieldだけでは一般にtotal spaceの有理性は従わず、generic G-torsorの
問題が残る。本証明ではactual rational sectionとgeneric freenessが、そのgeneric
productを直接与えた。quotientrationalityだけからtorsor trivialityを仮定しない。
Sp4がspecialであるという別の定理を使う必要もない。

split Sp4はQQ上のbig Bruhat cell

`Uminus×T×Uplus ≅ A4×(Gm)^2×A4`

を持つので、有理10次元である。従って上のbirational productにより

`QQ(B_NEW)=QQ(Sp4)(xi,r,s,q,B) ≅ QQ(t1,...,t15)`。

**B_NEWはQQ上で有理15次元、従ってC上でも有理15次元**である。
保存されたrho(g)=g*は有限pullbackのanti-actionである。この規約では
`Phi(g*h,z)=rho(h)*Phi(g,z)`、source第一factorの**右乗法**に対応する。
左乗法equivarianceとは主張しない。この規約差はbirational productの証明を変えない。

作用mapは既存rhoと保存beta_barの有限式で定義される。今回は15変数の全展開や、
group coordinateのexpanded inverse formulaを新たに保存したという主張ではない。
birational inverseの存在は上のdegree1 argumentによる。

## 5. 範囲と保存

結論はreduced affine firstjet hypersurface B_NEWについてのgeneric birational主張。
全点stabilizer、global quotientの積同型、boundary orbit分類、projectivizationの
stabilizer/有理性、Poisson成分25全体の有理性、source26のfield index、NEW方程式の
具体式・正確な次数、全component/exhaustionを決定しない。

checkerは三つのFINAL独立レビューとauthor receiptsを読み、同じpullback・section・
freenessというgateを照合し、input hashesとproduct-map recipeを保存するだけである。
新しい点、prime、minor、source26、不変式展開、CAS計算は実行しない。
