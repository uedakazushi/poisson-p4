# rank零plane・globalVの二十一whole低parameter源はG非優越

2026-11-07、主義務G。producerの全source証明。独立review前はglobal採録を待つ。

**命題。** 二つのfirst-jet notesで列挙した全normalnil b非零source、及び
`PLANE_RANK_ZERO_GLOBAL_VECTOR_LOW_PARAMETER_WEIGHT_BOUND.md`に指定した二十一の
whole低parameter sourceについて、全higher Poisson tensorのambient projective像はdimension二十四以下。
normalss square型とnormalnil b零regular-nil型のK零を含み、common-kernelの零Kも含む。
全m二plane、generic二/三parameter源、又は全Gを尽くしたとは主張しない。

## 1. 既に閉じた部分と、五つのK零fiber

`PLANE_RANK_ZERO_GLOBAL_VECTOR_LOW_PARAMETER_WEIGHT_BOUND.md`はfull350×30/9 kernelsと
全exact parameter rootを用いて十六whole models及び二つのpunctured modelsを閉じた。
当時の五つの未閉鎖fiberと実行時native serialization障害はimmutable report内に保持する。
`PLANE_RANK_ZERO_GLOBAL_VECTOR_ONECUT_BOUND.md`はそのうち二fiberを、
fullD2を消すexplicit cokernel functionalによる非零二次条件で全像二十四以下へ閉じた。
本noteでは残る三fixed sourceを、同じ元のquadratic-bivector Jacobiの全weight三/四で閉じる。
jet次数又は未知bivector次数を追加していない。

## 2. actual有限shearとfull weight二source

base change`x→x+C(z,w)`はplaneを保つactual unipotent PGL5変換で、逆はminus。
normal weight一の全三十canonical係数に対する作用は

`Π1→Π1+S*y`

であり、Sはreportのactual30×6 shear matrixである。
weight零は変わらず、残りは全weight二へ入る。
これはnormal weight一のvector fieldの有限changeなので、first-weight作用はその一次項に等しい。
projective Euler gaugeでcanonical representativesへ戻した場合も同じ商の作用である。

Sのindependent columnをBとし、kernel D1内でcomplement Rを取り、
`ker D1=im B⊕im R`とする。
reportは全matrix、選んだbasis minor及びinverseを保持する。
任意のfirst-weight sourceはactual有限shearでRのsourceへ移せる。
したがって全higher source dimensionは、normalized sliceのdimensionにrank Sを加えた値以下。
sliceを取ったという理由だけで非線形equationsを省略しない。

ss/nil commonではrank S四、slice dimension八。
nil T-squareではrank S五、slice dimension六。
各sliceのΠ1に対し`Q2=(1/2)E∧[Π1,Π1]`、weight二は
`D2*r+Q2=0`である。
全D2 rank六のconstant rational minorで六係数を解き、**kernel D2三係数v0,v1,v2を全て戻す**。
残るcompatibility conditionはfullD2のcokernelにQ2を射影した全式。
その必要式の選択と全cokernel functionalも保存した。

## 3. ss common：四whole linear branchesとheight二

eight slice variablesをu0,…,u7とする。full weight二compatibility idealは

`u7*(2u2+u3+2u5)=0`, `u4*(-2u0+u1-2u6)=0`

という二本で生成される（非零constant multipleを省略した）。
従ってその全被約台は次の四linear spacesの和。

1. u7=u4=0。
2. u7=0、u1=2u0+2u6。
3. u3=−2u2−2u5、u4=0。
4. u3=−2u2−2u5、u1=2u0+2u6。

零枝と全intersectionを含む和であり、open chart又はgeneric branchへの制限ではない。
各linear sourceはfirst-weight六係数、weight二kernel三係数、合計九係数を持つ。
full inhomogeneous weight二solutionを上のsubstitutionへ戻して、全三百五十weight三と
全三百五十weight四Jacobi係数を計算した。
各枝でweight三のrow65とrow66の二polynomialのgcdはQQ上で1である。
全polynomialsとexact gcd chainをreportに保存した。

polynomial UFDのheight一primeはprincipalである。
二式のgcdが1なら、それらを含むheight一primeは存在しないので、全Jacobi idealはheight二以上。
**unit ideal又は空とは結論しない。** このdimension boundはQQからCへのfield extensionで保たれる。
従って各normalized sourceのdimensionは高々七。
actual shear四、symbol base八、plane六、projectivization一を戻し

`dim image≤7+4+8+6−1=24`。

## 4. nil common：三whole linear branchesと一cut

full compatibilityは

`u7*(u3+u7)=0`, `u7*(-u1+2u2+3u5)-u3*u5=0`。

u7零又は非零を式の中で分ければ、全被約台は次の三linear spacesの和である。

1. u7=u3=0。
2. u7=u5=0。
3. u3=−u7、u1=2u2+4u5。

ここでも全intersectionと全零係数を保持する。
各first-weight六係数と全weight二三係数に対し、weight三row65は非零polynomialである。
full七百係数と非零witnessをreportに保持する。
従ってnormalized sourceのdimensionは高々八。
このnormalnil sourceのsymbol baseは七であるから

`dim image≤8+4+7+6−1=24`。

## 5. nil T-square：全reduced compatibilityは一つのlinear space

slice variablesはu0,…,u5である。
七本の全compatibility generatorsをg0,…,g6とし

`L=−u2+3u3+u4`, `d=u2−u4`

と置く。literal exact identitiesは

`g1=−12u5²`, `g3=u3*L/3`,

`g5=14u0*u5/3+L*(2d+9u3)/6`。

任意のreduced pointでu5=0、u3L=0、L(2d+9u3)=0となり、

`L²=(15/2)u3L−(1/2)L(2d+9u3)=0`

からL=0を得る。
逆にu5=L=0を七generator全てへ代入すると恒等零。
よって全reduced compatibility supportは**単一のlinear四次元space**であり、
native GBのscheme-theoretic dimension上界を強くする未証明のradical仮定は不要。

u5=0、u2=3u3+u4を代入して全weight二を解き、kernel三係数もすべて戻した七係数spaceで、
weight三row65は

`u3*(2u0²−v1)`

という非零polynomialになる。
したがって全normalized source dimensionは高々六。
actual shear五、symbol base八、plane六、projectivization一を戻し

`dim image≤6+5+8+6−1=24`。

## 6. exact証明書と未解決補集合

以下の全reportはimmutableで、各実行時producer全bytesを保存する。

- `CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_LOWPARAMETER_WEIGHTS_20261107.json`：全二十一source/全linear operators/全parameter roots。
- `CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_ONECUT_20261107.json`：二つのK零fiberのcokernel obstruction。
- `CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_COMPATIBILITY_20261107.json`：全actual shear slices、full Q2、compatibility generators。T-squareのnative印字障害も保持。
- `CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_RESIDUAL_20261107.json`：common全七linear branchesのfull original Jacobi、gcd/nonzero証明。T-squareは同一QQ coreを正しく印字して一度だけ修復。
- `CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_TSQUARE_20261107.json`：T-square全reduced support、全weight二solution及び非零weight三witness。

full source geometry、base dimension、finite shearの橋は本noteと二つのfirst-jet notesにある。
実行時間は順に24.116、2.426、7.663、4.718、1.007秒。
reusable runnerにisolated modesをappendし、各回全旧FunctionDef ASTが直前immutable producerと一致することを確認した。
旧parameter certificate/report/manifestを現source bytesへ書き直していない。

これで二十一whole低parameter源の五つのK零frontierは全て閉じた。
normalss generic MIX三parameter、MIXH/MIXN/NILREG二parameter、
normalnil b零Cartan三parameter及びgeneric mixed二parameterは別の全source義務である。
全m三以上も新modular入口でnormal weightに入るが、その全復活normal formsは本命題に含めない。
全G、全E又は全分類の完了とは主張しない。
