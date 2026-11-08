# contact31_chord全六磁気源のSNC実衝突：独立全源review

2026-11-08（UTC）、PRIMARY E、`/root/log_manuscript`。
`TWO_QUADRICS_CHORD_SNC_CONFLUENCE.md` のwhole source証明を別作者として読解し、
親全parameter、実変換、端点全六源と幾何的閉包の橋を独立確認した。**PASS**。
結論は `contact31_chord` の全点が既MPS40有限閉和に属することだけであり、成分数は増えない。
他のtan-simple source、特定T番号、局所一意性は結論しない。

## 1. 一回の独立Fraction検証と正しい凍結input

producer入口は `inputs/two_quadrics_chord_snc_confluence_manifest_20261107.json`
（SHA256 `9ffc5ebd94c5bd07e1a4d34155df6a236e2ad6f667746aaec19697aa44906e13`）。
原全sourceは `CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json` の `contact31_chord`、
括弧順はxy,xz,xw,yz,yw,zwである。

共有producerは後続tan-simple modeの追加により既にadvanceしたため、**mutable current runnerをこのreviewのinputにしない**。
chord report内のexecuted `producer_source` と
`inputs/two_lines_pre_tan_simple_source_contracts_20261107.json` のexact UTF8 snapshotが同じ
`2f3c3b6741f6ba8a05a9460eebc36410850c486865d0f6172ae7a19ea92411bb`
であること、当時の全保存契約hash、manifestの他全不変inputを照合した。
producerを実行又はimportせず、旧bytesも変更していない。

新checker `classification_two_quadrics_snc_review.py` は過去の独立sparse Laurent/Fraction演算だけを再利用する。
producer、S61 engine、CASはimportしない。
有限job `jobs/two_quadrics_chord_snc_separate_review.json` をgate PASS後、一回の全parameter計算で
**0.092秒**で次を完了した。正確な実測値・RSS・実行source全文・入力hashは新reportに保存する。

* 全五parameter親の全十strict Jacobi、元projective全350、Euler-signed actualDと真SNC特殊化。
* 全period/edgeの実product関係を残したconfluence、全行列と逆、det epsilon²、全負epsilon指数の消滅。
* 全六磁気の実translation/fiber scaleと正逆、6×6parameter map及び逆、非零定数det。
* 全epsilon・全磁気の元projective350/strict10とactualPf恒等式、全150端点係数の原sourceとのliteral一致。
* 六磁気微分がh²,hx,hy,x²,xy,y²を全て与えること、全磁気零を含むregular性。

この機械検証とは別に、次節以降のirreducible親とSNC-component theoremの適用範囲を手で確認した。

## 2. SNC点に会う全既約親からのMPS帰属

H,A,B,C,Dで

```
P=p(AB+uH²)∂A∧∂B+q(CD+vH²)∂C∧∂D+kN1∧N2,
N1=A∂A−B∂B, N2=C∂C−D∂D.
```

両二次元Poisson tensorは別座標に依存して互いに可換、N1,N2は両tensorを保存し互いに可換である。
従って全A5でstrict Poisson。
cross k²項はPfaffianで相殺し、actualDは
`Hpq(AB+uH²)(CD+vH²)`。
p=q一、k=u=v零は非零rank四tensorで実因子D=HABCD、五独立超平面のSNCである。

親の非零tensor openは既約で、射影像の閉包も既約。
全Poisson空間の有限個のambient componentsのうち一つCがこの全閉包を含む。
このCは上のSNC点も含む。
既受理MPS Theorem1.4の「SNC log-symplectic locusに会う全ambient component」というscopeにより、CはMPS40の一つ。
よって親の全非零像は既知有限閉和に属する。
SNC点だけのknown-union membershipから任意曲線へ逆推論したのではなく、
**全親像を含む同じambient component** がその点にも会うことを使う。
SNC点の成分一意性、component IDの同定又は親sourceのlocal smoothnessは不要である。

## 3. 実PGL衝突と全pole cancellation

親u=v負一、`p=−1+epsilon*s,q=1,k=−1+epsilon*t` とし

```
(H,A,B,C,D)=(h,x,y,x+epsilon*z,y+h+epsilon*w).
```

このtarget→parent行列Gはdet epsilon²、reportの逆行列と全両方向で一致する。
target tensorは `epsilon*G^-1 P(G target) G^-t`。
epsilon非零で実PGLと非零scalarだから既MPS閉和に属する。
唯一現れ得る負指数係数は

```
(p+q)(xy−h²)+(q+k)hx=epsilon*(s(xy−h²)+t hx)
```

で正確に消える。p,q,kの積を独立Taylor係数へ緩めていない。
全curveはepsilonについてpolynomialで、そのbare端点は

```
{x,y}=0, {x,z}=−x², {y,z}=h²,
{x,w}=−h²+hx+2xy, {y,w}=−hy−y²,
{z,w}=(h+2y)z+s(xy−h²)+t hx.
```

## 4. 全六磁気端点と可逆parameter map

M=m0h²+m1hx+m2hy+m3x²+m4xy+m5y²、c=−5/7とする。
全epsilon族へ同じregular translation `z'=z+ell,w'=w+kappa`、最後にw方向をc倍する。

```
s=−7m4/5, t=−7m1/5,
a=7m2/10−7m5/15, b=7(m0+m4)/5−a,
d=7m5/15, e=−7m3/5,
ell=a h+b x+d y, kappa=e x.
```

translationによる磁気差は全係数で
`−(a+b)h²−2(a+d)hy−3d y²+e x²`。
これをs(xy−h²)+thxへ加えてc倍するとMそのものになる。
独立に再構成したm→(s,t,a,b,d,e)の6×6matrixはdet `117649/93750`、保存逆との両積が一。
全frameのdetは−5/7、inverseはtranslationを引き、w尺度を−7/5倍する正則式である。
磁気parameterは一つも反転しない。

端点の全150係数は原contact31_chordと厳密一致し、特に

```
{x,z}=−x², {y,z}=h²,
{x,w}=−5(−h²+hx+2xy)/7, {y,w}=5y(h+y)/7,
{z,w}=M−5z(h+2y)/7.
```

actualDも全curveで直接確認した。
old z=z−ell、old w=−7w/5−kappaと書けば

```
D_epsilon=h(xy−h²)[(x+epsilon*old z)(y+h+epsilon*old w)−h²]
             *(5/7+epsilon*m4),
D_0=(5/7)h(xy−h²)(−h²+hx+xy).
```

したがって全磁気零を含む任意のmで、regular endpointは既MPS有限閉和に入る。
途中でpが零となるrank低下も親の全像に留まり、parameter補集合を落とさない。
全PGL飽和も同じ既知閉和に属する。

## 5. 保存と残余

新reportは `CLASSIFICATION_TWO_QUADRICS_CHORD_SNC_SEPARATE_REVIEW_20261108.json`。
UTC日付が変わったため新独立実行を11/08のidentityで保存し、11/07 producer測定を再現したと装わない。
新component、特定のMPS番号、局所unique、他source又は全E分類の結論は無い。
`contact31_tangent_at_simple_intersection` はこのreceiptの対象外である。
