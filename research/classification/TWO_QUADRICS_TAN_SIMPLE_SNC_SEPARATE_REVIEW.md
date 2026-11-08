# tan-simple全六磁気のcoupled SNC衝突とnamed T16/T17：独立review

2026-11-08、PRIMARY E、`/root/log_manuscript`。**全source算術及び幾何bridgeはPASS**。
対象は `contact31_tangent_at_simple_intersection` の全六磁気源であり、
既知MPS40有限閉和、さらに下記の全括弧同定によりT17へ入る。
先のchord源もT16と同定する。新成分・局所unique・全分類の結論は無い。

## 1. 一回の有限独立確認

producer入口は `inputs/two_quadrics_tan_simple_snc_confluence_manifest_20261107.json`
（SHA256 `786dc644...`、全値はmanifest自身に固定）。
実行sourceはreport内の全文、SHA256
`3c7d5eb884d27e5a2425336ab23d3dc866deff00b4ee8da1117e5ae6a61d091f` へbindし、実行又はimportしない。
`classification_two_quadrics_snc_review.py --tan-simple` は過去の独立Fraction/Laurent ringだけを用いる。
CAS・producer・S61 engineをimportしない。

旧reviewerと旧chord current jobのexact bytesは
`inputs/two_quadrics_snc_pre_tan_review_source_contract_20261108.json` に保存した。
旧全FunctionDef（parent/main）のASTは変更せず、isolated modeだけを追加した。
旧chord計算は再実行せず、current jobのcompatible input bindingのみ新sourceと旧archiveへ更新。
旧chord reportと旧manifestのbytesは不変で、その旧job/source hashは上のarchiveからexactに回収できる。

新job `jobs/two_quadrics_tan_simple_snc_separate_review.json` gate PASS後、
一回 **0.274秒**で以下を全parameterのまま確認した。

* coupled親の全五parameter、strict十式・元350式・actual Pfaffianと真SNC点。
* 全PGLと逆、det−epsilon²、全負epsilon指数の消滅。
* 全六磁気parameter matrix/逆、det−1769472/244140625、実fiber frame/逆、det−25/8。
* 全epsilon・全磁気のstrict十式/元350式/actual Pfaffian。
* 原sourceの全150端点係数、全六磁気微分h²,hx,hy,x²,xy,y²、零磁気境界。
* raw T16とcompleted T17の全五parameter、明示Euler gauge、全六括弧及びhelperへの逆binding。

新report `CLASSIFICATION_TWO_QUADRICS_TAN_SIMPLE_SNC_SEPARATE_REVIEW_20261108.json` は
正確な実測値・RSS、実行reviewer全文、全input hashesを保持する。

## 2. coupled親の全関係を保つSNC bridge

親の全五parameterはp,q,a,u,v、b=p−a、
`F=pAB+uH²`、`G=qCD+vF` で、

```
P=F∂A∧∂B+G∂C∧∂D
 +(A∂A−B∂B)∧(aC∂C+bD∂D).
```

JacobiのACD項は `vpAF+pqACD−(a+b)A(qCD+vF)`。
他の三非自明項も同じa+b=pにより零となる。
Fの各係数、vFの積を独立parameterへ緩めたsourceではない。
actualDはHFG、p=q一、a=u=v零でHABCDの真SNC rank四点。

全親parameter空間は既約なので、非零射影像の全閉包を含む同じambient componentがこのSNC点にも会う。
既受理の「SNC log-symplectic locusに会う全成分」のMPS40定理により全親像が既知閉和へ入る。
SNC一点の既知帰属を任意変形へ逆外挿せず、局所一意性も用いない。

## 3. 正則全衝突と磁気parameter inverse

実座標変換及びperiod bindingは

```
(H,A,B,C,D)=(h,y+h,x−y+h,x+epsilon z,−y+epsilon w),
p=1,u=−1,v=−1/3,
q=−1/3+epsilon q1,a=2/3+epsilon a1.
```

target tensorはepsilon倍した実pushforwardで、epsilon非零ではPGL同値。
det−epsilon²、唯一可能なnegative poleは元coupled関係を戻すと全て消える。
bare endpointは

```
U=(2x(2y−x)/3, −(hx+xy−3y²)/3),
V=((y²−2xy−3hx)/3, −y(y+h)/3),
{z,w}=−(2h+y)z/3−2yw/3−q1xy−a1(hx+2y²).
```

目標M=Σmi(h²,hx,hy,x²,xy,y²)iに対してni=−8mi/25と置く。

```
a1=−(n1+n5)/3,
q1=−n4−n3−n5/3+2n1/3,
A=3n0/2, B=n5−2n1, C=3n2−A, E=3n3/2.
```

translationはell=Ah+Bx+Cy、kappa=Ex、z'=z+ell,w'=w+kappa。
磁気差は全係数で

```
2Ah²/3+(A+C)hy/3−Bhx/3+(−B−2E)xy/3+By²/3+2Ex²/3.
```

これをbare磁気へ加え、定数fiber matrix
`T=[[-5/4,0],[25/24,5/2]]` を施すと、そのdet−25/8倍で正確にMとなる。
全6×6parameter mapのdetは−1769472/244140625、保存逆との両積は一。
translationとTの逆は全磁気係数上でregular、磁気parameterの反転は無い。
全150端点はaccepted sourceそのものと一致し、actualD端点は
`−25h(hx−y²)(hx+xy−y²)/24`。
全curveのactualDも全五Pf項から再構成して保存factorizationと係数比較した。

任意のm、特に全零mでもregular曲線の一般epsilonは既MPS閉和に属する。
その閉性により端点全源と全PGL飽和が被覆される。

## 4. 二つのnamed parentの独立全六括弧確認

`TWO_QUADRICS_CONFLUENCE_NAMED_MPS_PARENTS.md` の依存sourceと全gaugeを別に確認した。
原順(A,B,C,D,H)から(H,A,B,C,D)へ並べ替え、raw PからE∧Lを引く。

T16では
`L_B=(rho1−rho2)B/2,L_C=(rho1+rho2)C/4,L_D=−(rho1+rho2)D/4`、他零。
全六括弧はchord親と一致し、全polynomial bindingは

```
(rho0,rho1,rho2,tau0,tau1)
=(q+2k,−p−2k,p−2k,pu,qv).
```

raw sourceの全五parameterを戻してから比較したため、特定period一点の認識ではない。
従ってchordの全親と既確認の全六磁気端点はT16へ入る。

T17では **completed source** の `tau0 tau1 H²/rho1` 項を保持した。
`L_B=−rho1B,L_C=−(rho1−rho2)C/2,L_D=−(rho1+rho2)D/2`、他零。
全六括弧はcoupled親と一致し、

```
p=rho1, q=rho0+rho2, a=(rho1−rho2)/2,
u=tau0, v=tau1/rho1,
(rho0,rho1,rho2,tau0,tau1)=(q−p+2a,p,p−2a,u,pv).
```

この全正逆bindingも実際に全括弧へ代入して一致を確認した。
named source分母rho1のopenはp非零に対応し、helperのp零はpolynomial性とT17の閉性で戻す。
今回の実衝突ではp一を固定しており、全epsilonでrho1は単元。
従ってtan-simpleの全六磁気端点はT17へ入る。

各完成sourceのPGL像閉包がnamed成分であるという入力は、既受理
`MPS_S61_CATALOGUE_20261103_v62.md` §§1–2に明示されたperiod/edge全族と局所成分認識に置く。
Pfaffian型又は次元だけから番号を推測したのではない。
元のSNC全component論証はこのnamed同定に依存せず、それ自体で有効である。
