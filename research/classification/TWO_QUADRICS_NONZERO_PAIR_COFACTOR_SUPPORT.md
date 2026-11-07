# H＋二既約quadric、非零log括弧の全支持はMPS

2026-11-07。課題E_reduced_nonnormal。対象は実退化因子
`D=h q₁q₂`、q₁,q₂は相異なる既約二次式、`Pf P=κq₁q₂`, κ≠0。
`TWO_QUADRICS_SUPPORT.md` 第2節の**全支持**から得た二つのnormal場を使う。
ここでは定数 `λ=X₁(log q₂)=−X₂(log q₁)` が非零の場合を全て閉じる。
λ=0を除外したまま全1+2+2又はGを完了したとはしない。

## 1. 三つだけのjoint-character配置への完全還元

全tensorの非零scaleでλ=1とする。五次元の線形形式空間Vにおいて
二つの可換場はhを固定し、q₁,q₂のjoint-characterはそれぞれ

`τ₁=(0,−1)`, `τ₂=(1,0)`、`tr(X₁,X₂)=τ₁+τ₂=(1,−1)`。

可換なJordan分解を取る。半単純部はhを固定し、qᵢに同じcharacter τᵢで作用する。
冪零部はqᵢを固定し、両半単純部と可換である。従ってhを第一要素とする
joint-weight基底を選べる。rank五は自身を保存する場のtraceが零となって矛盾するため、
qᵢのrankは三又は四だけである。rank四・三の順序を交換してもよい。

この時点でweightの相異性やhの非radical性を仮定しない。
weightを `t₀=0,t₁,…,t₄` とすると、全trace条件は

`Σⱼ₌₁⁴ tⱼ=τ₁+τ₂`。                                          (1)

rank rの対称行列には非零のr次principal minorがある。
そのdeterminantの非零置換項の各cycleを考えると、辺に沿うweight和は全てτである。
偶数cycleは隣接辺の組へ、奇数cycleは全頂点の固定点へ置き換えられる。
従って、あるr座標の上にinvolution matchingが存在し、各辺(i,j)について

`tᵢ+tⱼ=τ`（固定点なら `2tᵢ=τ`）                         (2)

が成立する。これは重複weight、任意のradical配置を含む必要条件である。

五座標のr部分集合と全involutionは、r=4で50個、r=3で40個。
両quadricについての(2)と(1)を全て厳密有理行消去すると、(4,4),(4,3),(3,3)で
それぞれ2500,2000,1600組のうち整合するものは24組ずつだけである。
全て解は一意であり、四座標の置換を除いて次表の一組になる。
rank不足の未定weight系や重複weight系は一つも残らない。
この有限全列挙はproducerから独立したFraction reviewerに保存した。

| rank | (tₓ,tᵧ,t_z,t_w) | q₁ | q₂ |
|---|---|---|---|
| 4,4 | (τ₁,τ₂,τ₂−τ₁,τ₁−τ₂) | hx+yw | hy+xz |
| 4,3 | (τ₁,τ₂/2,τ₁−τ₂/2,τ₂−τ₁) | hx+yz | y²+xw |
| 3,3 | (τ₁/2,τ₂/2,τ₁−τ₂/2,τ₂−τ₁/2) | x²+yz | y²+xw |

いずれもhを含む五joint-weightは全て相異なる。
従って可換冪零部は各一次元weight空間上で零であり、元の二場そのものが対角である。
このweight一覧では各quadricに許される単項式は表の二項だけ。
実際のrank三又は四から両係数は非零であり、座標と各quadricの非零scaleにより
両係数を一に出来る。従ってこの表は一つの開集合だけでなく**全λ≠0支持**を覆う。

## 2. 全90係数の正確な源と全Jacobi/Pfaffian

各表でh=1とし、Pの六括弧に次数二以下の全15単項式を入れる。
二つのnormal条件、発散、両不変条件を全て課す。

`−P·dqᵢ=qᵢXᵢ`, `div P=X₁+X₂`, `[Xᵢ,P]=0`。             (3)

第一式の符号は第一添字縮約の規約による。
これは280+20+180=480行、90列のaffine線形系である。
各型のrankは89、自由係数は一個n。全480行と右辺、全90係数の完全解、
非零89次minorを保存し、別の標準ライブラリFractionコードで再構成した。
従って許容源の部分族だけを計算したのではない。

各完全源では全四strict Jacobiが恒等零で、Pfaffianは次の通り。

| rank | 全Pfaffian | rank四の条件 |
|---|---|---|
| 4,4 | −n q₁q₂ | n≠0 |
| 4,3 | (1/2−n)q₁q₂ | n≠1/2 |
| 3,3 | (3/4−n)q₁q₂ | n≠3/4 |

零Pfaffianのparameterも全源内に保持する。対象のrank四条件だけで上記一点を除く。
別のrank、共鳴、minor零の補集合はこれ以外に存在しない。

## 3. 各全源の実SNC軌道極限

次の整数cocharacterを(h,x,y,z,w)へ作用させる。

| rank | weight(h,x,y,z,w) | 極限Pfaffian |
|---|---|---|
| 4,4 | (3,0,0,0,0) | −n xyzw |
| 4,3 | (3,0,1,0,0) | (1/2−n)xyzw |
| 3,3 | (3,0,0,−1,−1) | (3/4−n)xyzw |

各完全源の**全単項式**のbivector weightは非負で、weight零の全括弧は
定数倍の `xᵢxⱼ` だけになる。表のκが非零ならその極限は非退化logcanonicalであり、
実退化因子は五つの独立超平面 `hxyzw` のSNCとなる。
負weight項を特定parameterで零と仮定したり、一点だけを退化させたりしていない。

元の点を含む任意の既約成分は連結PGL₅で不変であり、同じ成分はこの軌道極限も含む。
MPS Theorem1.4のSNC点での全成分分類により、その成分は既知MPS40の一つである。
従って全λ≠0支持の全rank四点は既知MPS40の閉和に入り、未知成分を作らない。
これは全源を対象にした結論で、特定MPS番号又は局所唯一性は必要としない。

## 4. 証明書と残る範囲

* producer：`code/classification_two_quadrics.py --nonzero-pair-cofactor`、3.274秒。
* 独立再構成：`code/classification_two_quadrics_review.py --nonzero-pair-cofactor`、0.511秒。
* `reports/CLASSIFICATION_TWO_QUADRICS_NONZERO_PAIR_COFACTOR_20261107.json`。
* `reports/CLASSIFICATION_TWO_QUADRICS_NONZERO_PAIR_COFACTOR_REVIEW_20261107.json`。
* 実行前gate：`jobs/two_quadrics_nonzero_pair_cofactor{,_review}.json`。

reviewerはproducer/CASをimportせず、全480行、完全affine解、非零minor、
全四Jacobi、実Pfaffian、全parameterのSNC極限と、上記6100 matching組を独立再構成した。
matchingから幾何への橋及びMPSの入力はこのノートの第1,3節に明記した。

λ=0には二つの独立可換なtrace零場がh,q₁,q₂を全て保存する。
generic pair-zero半単純枝とnil22-image枝の既済結果を合わせても、
他の共鳴・冪零・pencil/hyperplane配置をまだ全被覆していない。
また平面七配置のうちtan-simpleとchordの全源帰属も未完了である。
この二源は次元≤23によりGの該当未知優越成分を作らないが、Eの認識は別に残る。
