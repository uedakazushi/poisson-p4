# 全 nonnormal integral CI22 rank 二枝の独立監査

2026-11-07。主義務 G。全五 surface、full150 tangent/conormal 源、全被約支持の
逆場合分けと actual D の零補集合を別作者として監査した。
標準ライブラリ Fraction と既済の独立多項式ライブラリだけによる再構成は
**全項目 PASS**。producer、SymPy、native CAS は import も実行もしていない。
主原稿・global 台帳・Git・旧証明書は変更しない。

## 1. 結論と範囲

S⊂P4 を integral・非退化・nonnormal な CI(2,2) とする。
S を非正規曲面支持に持ち、S 上の一般 restriction rank が二で、実退化因子 D が
既約被約 quintic である全 rank 四 Poisson 構造の PGL5 飽和射影像は二十四次元以下。
零 divisor R の既約・被約・vertical 条件を付けない。計算源は rank 零 restriction
も安全な過大評価として含む。

normal CI22 rank 二は別の受理済み空性定理である。それと合わせると全 integral
非退化 CI22 rank 二枝から未知二十五〜二十七次元 G 優越成分は出ない。
全 E 帰属、点の不存在、新成分性、全 G の尽くしは本監査の結論ではない。

独立実装は `code/classification_nonnormal_ci22_rank_two_review.py`、契約は
`jobs/nonnormal_ci22_rank_two_review.json`、結果は
`reports/CLASSIFICATION_NONNORMAL_CI22_RANK_TWO_REVIEW_20261107.json`。
実行前 gate PASS、最終実行 1.614 秒、最大 RSS 17,296 KiB。

## 2. 全五 surface orbit を尽くす幾何

外部入力は既受理の**任意の複素 symmetric pencil**の congruence canonical theorem
だけである。generic pencil の closure theorem を全 pencil に拡張していない。
本監査は既存の一次資料 receipt を背景として使い、その適用後の全 size 五場合分けを
独立に読み直した。今回新たに原典全文を取得したとは主張しない。

integral 非退化 S を含む pencil に rank 二以下の member は無い。linear factor
の積が prime I に入れば一つの linear factor が入り、非退化性に反するためである。
非退化 member Q があれば Q-self-adjoint M を使える。従って各固有値の eigenspace
は二次元以下。singular locus が curve となるのは、同じ固有値の二 Jordan block の
長さが双方二以上で、その二次元 kernel 上の Q が恒等零の場合に限る。
size 五では `(22)+(1)` 又は `(3,2)` だけである。固有値は二つ又は一つなので
pencil PGL2 により連続 parameter は残らない。

全 member 退化なら singular blocks Md の size は2d+1。M2 だけの CI は plane
を含み integral 非退化でない。M1+regular2 は distinct eigenvalues の一型だけが
残る。同一 eigenvalue の diagonal 型は rank 二 member、Jordan2 型は plane を
含む。common radical 二次元以上も三次元以下の残 pencil に rank 二 member を持つ。
M0+M1+regular1 も regular 部分を消す member が rank 二となる。
従って残りは common radical 一、regular4 の cone だけである。

cone の base C⊂P3 は integral arithmetic-genus 一 quartic。
各 eigenvalue の block は一つでなければ rank 二 member がある。
nonnormal cone では C が singular。δ の和は一以下なので二つの singular point
は無い。Jordan4 型は standard equations に line を含み integral でない。
従って node `(2,1,1)` 又は cusp `(3,1)` だけ。三又は二 eigenvalues を pencil PGL2
で固定すると、それぞれ一つの actual PGL5 surface orbit である。

| model | q1 | q2 |
|---|---|---|
| regular22_1 | `2x0x1+2x2x3+x4²` | `x1²+x3²+x4²` |
| nil32 | `2x0x2+x1²+2x3x4` | `2x1x2+x4²` |
| singularM1 | `x0x1+x3²` | `x0x2+x4²` |
| cone_node | `2x0x1+x2²+x3²` | `x1²+x2²+2x3²` |
| cone_cusp | `2x0x2+x1²+x3²` | `2x1x2+x3²` |

integrality も確認した。regular22_1 と cone_node は x1≠0 で q1 が x0 を消し、
残 q2 は irreducible rank 三 conic。nil32 と cone_cusp は x2≠0 で x1,x0 を
順次消し、singularM1 は x0≠0 で x1,x2 を消す。各補集合は surface 内で次元一以下
（cone の base では一点）。CI の unmixed 性と generic reduced 性から prime。
共通 vertex 座標を加えても integrality は保たれる。

## 3. normality を使わない full150 源

全 raw quadratic skew coefficients 百五十を保持し、両 qi について

```
E wedge (Pi#dqi) = 0 mod I,  I=(q1,q2)
```

の全 component を課す。degree 四の I は全 qi×quadratic の rank 二十九空間。
独立実装は exact QQ rowspace の剰余から全 tangency rows を再構成した。
normalization の anticanonical section を五係数へ制限する ansatz は使わない。

dense x0 chart の intrinsic generators qi/x0² の bracket numerator は

```
K=x0{q1,q2}-2q2{q1,x0}-2q1{x0,q2}.
```

rank 二の transverse Lie algebra に対し、実 D が singular surface を持つことから
generic conormal bracket は零となる。これは既済 abelian-conormal 入口の適用であり、
ordinary node 又は R の被約性を追加しない。prime CI の I² は I-primary なので、
x0∉I を使って dense chart の membership を全 homogeneous `K in I²` へ戻せる。
degree 五の I² は q1²,q1q2,q2²×全 linear の rank 十五空間。
全 conormal remainder matrix も元百五十係数から独立に再構成した。

全 Euler 二十五方向の rank は二十四で全条件の kernel に入る。
explicit coordinate gauge の非零 minor を再計算し、完全 kernel と全 free coordinate
inverse を検査した。従って conductor 上だけの section、Euler 方向又は chart 補集合を
落としていない。

| model | 全 tangent 射影 section 源 | conormal 後 affine 源 | 元 quadratic rowspace rank |
|---|---:|---:|---:|
| regular22_1 | 26 | 10 | 16 |
| nil32 | 26 | 12 | 17 |
| singularM1 | 26 | 10 | 14 |
| cone_node | 26 | 14 | 9 |
| cone_cusp | 28 | 16 | 20 |

ここで「射影 section 源」は Euler quotient の意味であり、最後の tensor scalar
projectivization の一はまだ引かない。全源の150 tensors、全350 Jacobi、
実 Euler-contracted D の全126係数、元 quadrics の完全 QQ rowspace を再構成した。
初回の coarse `all_G24_closed=false` は不変保存し、以下の全支持で補う。

## 4. regular22_1 の全零補集合

元 rowspace の p 座標を使う。`(p2-p8)²=0` より reduced support で p2=p8=t。
t≠0 では元 product 式から p3,p6,p7,p9=0、残式は p1=p0=0,p4=-t/2。
自由 t,p5 の二次元以下である。

t=0 では `p1(2p1-p4)=0`。p1≠0 なら p4=2p1、p3,p6,p7,p9=0、
自由 p0,p1,p5 の三次元以下。p1=0 なら残条件は正確に
`(p0,p4)*(p7,p9)=0` なので全支持は次の二つの五 parameter linear source に入る。

```
A: free(p0,p3,p4,p5,p6), others0;
B: free(p3,p5,p6,p7,p9), others0.
```

両方で全350 Jacobi と全126 actual D が**全係数で恒等零**となることを独立に
確認した。representative point からの外挿ではない。
従って rank 四既約 D の全 source は前二枝だけに入り、affine 次元<=3。
二乗式を scheme 上の linear equation と主張せず、reduced support に必要な方向だけを使う。

## 5. nil32 の全逆 case と九 supersource

元十七 quadrics の全 rowspace を再構成した。A=p10、P=p4、B=p8、r=p6、t=p3、
f=p5、g=p11 と書く。

A≠0 では元式から r=0,P=3A/2,t=3B/2、
`8Af=B²`、`6Ap1+(g+f)B=0`、`2Ap2+Ap9+gB+fB/2=0`。
さらに元式の p0 係数は非零定数倍 A なので p0 が決まる。
自由 A,B,g,p2,p7 の五以下。A=0,P≠0 では r=0,t=2B,g=-f、
p1=p0=0,p2=-2p9 となり、自由 P,B,f,p7,p9 の五以下。

A=P=0 では core は

```
3t²-5r²-3rB-12B²=0,
r(t+r-B)=0,
9tB-2r²-18B²=0.
```

r=0 なら t=2B。r≠0 なら t=B-r、`(3B+r)(3B+2r)=0`。
従って全 core は三 line
`(t,r,B)=(2s,0,s),(4s,-3s,s),(5s/2,-3s/2,s)`。
s≠0 の第一 line は g=-f,p1=0,p2=-2p9,p0=0 で四次元以下。
第二/第三 line は f=g=0 となり、それぞれ

```
(p2,p9,p0)=(2p1,-p1/2,p1²/(2s)),
(p2,p9,p0)=(7p1/3,-4p1/3,10p1²/(9s)).
```

自由 s,p1,p7 の三以下である。

s=0 の全 core 零では

```
f(g+f)=0,
12p1g+6p1f-8(g+f)p2-10(g+f)p9=0,
6p1²-3p1p9+4(g+f)p7=0.
```

g=-f≠0 なら p1=0、自由 f,p0,p2,p7,p9 の五以下。
f=0,g≠0 なら p2=(6p1-5p9)/4、p7=-(6p1²-3p1p9)/(4g)、四以下。
f=g=0 では p1(2p1-p9)=0 の両枝を保持し、それぞれ四以下。
これで各反転の零補集合と全交わりを尽くし、全 source<=5。

独立実装は保存された全九 Laurent maps を別 parser/多項式積で評価し、元十七式の
全残差を再構成した。反転されるのは case 名で非零とした単一 A,B,f,g だけである。
第一 map などは残る非零 compatibility を持つ**必要 supersource**であり、全 map 値が
Poisson だとは主張しない。coverage は上の逆 case 証明から得る。
map 代入が成功しただけから全支持を推測していない。

## 6. singularM1 の全 case は五次元以下

p8≠0 なら元 product 式で p0,p1,p2,p3,p7,p9=0、p4=p5p6/p8。
自由 p5,p6,p8 の三以下。
p8=0 の全 rowspace は正確に次の八式と同値であることを、独立 QQ rowspace で確認した。

```
p0p5+p3p6=0,  p1p2+p6(p3-2p5)=0,
p1(p3-3p5)=0,  p9(p3-2p5)=0,
p7(p0+2p6)=0,  p2(p0+3p6)=0,
p1p7=0,  p2p9=0.
```

p5≠0 で x=p3/p5 とすると p0=-xp6。
x≠2,3 なら p1=p9=p6=0、自由 p5,x,p2,p4,p7 の五以下。
x=2 なら p1=0、p2p6=p2p9=0。p2=0 は五、p2≠0 は四以下。
x=3 なら p9=0、p7p6=p1p7=0、p1p2+p5p6=0。
p7=0 は五自由変数上の非零一式により四以下、p7≠0 も p1=p6=0 で四以下。

p5=0,p3≠0 なら p6=p1=p9=0、p0p7=p0p2=0 なので最大四。
p5=p3=0 で p2≠0 なら p1=p9=0,p0=-3p6,p7p6=0、最大三。
p2=0,p7=0 なら自由 p0,p1,p4,p6,p9 の五。
p2=0,p7≠0 なら p1=0,p0=-2p6、最大四。
従って全零補集合を含め source<=5 である。

## 7. cone の全 actual D 零補集合

cone_node では全126 actual D の非零三十三 row が、元350 Jacobi の QQ-linear
combination となる。各 original-row functional の全係数を独立に検証した。
従って全 Poisson 点で actual D=0、rank 四枝は空である。

cone_cusp では元式 p13*p15=0。
両方零の全源では全126 actual D が恒等零である。
p13≠0 なら p15=0、元式から

```
p0=p1=p2=p3=p5=p7=p10=p15=0,
p4=2p12, p8=-p6, p14=-p9.
```

自由 p6,p9,p11,p12,p13 の五。
p15≠0 なら p13=0、

```
p0=p2=p5=p6=p8=p9=p11=p13=0,
p1=2p12/3, p4=4p12/3, p7=(-12p10+9p3)/2.
```

自由 p3,p10,p12,p14,p15 の五。
両全五 source の元350恒等零と、両零補集合の actual D 零を独立に全係数確認した。
全 open/zero case を保ち、generic cusp coefficient のみへ制限していない。

## 8. 全 actual Aut と射影像への橋

実装は full25 linear vector fields の `V(qi) in span(q1,q2)` という全係数 matrix を
別に作り、scalar を一つ除いた全 pencil-preserving Lie algebra 次元を計算した。
順に3,4,4,5,6を得た。標数零の algebraic stabilizer は smooth なので実際の群次元
を与える。最弱十分な既表示下界3,4,4,5,5を用いるだけで結論する。

幾何的にも regular22_1 の skew centralizer 二と pencil lift 一、nil32 の四 field、
singularM1 の rank 二 Veronese center stabilizer 四、cone の vertex affine 変換五が
actual group を与える。Lie 計算はこの作用下界も全係数から独立に確認したものである。

各 parameter wall が個別に Aut-invariant であるとは仮定しない。
まず全有限 cover から **stabilizer 不変な全 rank 四・既約 D source X** の最大 affine
次元 d を得て、その後 PGL frame の自由な group fiber と tensor scalar 一を引く。
Euler quotient は既に source に反映済みだが scalar projectivization は未消費なので

```
dim image <=24+(d-1)-dim Aut(S)
```

である。個別 tensor の stabilizer が零である必要はない。

| model | 全 rank 四 source d 上界 | 採用 Aut 下界 | 全射影像上界 |
|---|---:|---:|---:|
| regular22_1 | 3 | 3 | 23 |
| nil32 | 5 | 4 | 24 |
| singularM1 | 5 | 4 | 24 |
| cone_node | 空 | 5 | 空 |
| cone_cusp | 5 | 5 | 23 |

すべての source は raw full150 から得たため、R の nonreduced/vertical 境界や
conductor 上の余分な section を後で捨てる必要はない。以上で本 G 枝が閉じる。

## 9. 入力・実装履歴・再現

入口 `inputs/nonnormal_ci22_rank_two_manifest_20261107.json` の全十二入力と
各 executed producer snapshot hash を検証した。旧 fullsource/coarse report と
support-completion report は不変である。
独立 checker 初回は等しい RREF の dict 挿入順まで比較して安全停止した。
rowspace は厳密に同じであることを確認し、pivot-keyed dict 同士の比較へ直し、
自分の code hash と有限契約だけを更新して再 gate 後に全実行した。
失敗実行は証明として使用せず、producer bytes は一切変更していない。

```
python scripts/classification.py gate research/classification/jobs/nonnormal_ci22_rank_two_review.json
python code/classification_nonnormal_ci22_rank_two_review.py
```

全 source の有限算術、五 surface の網羅性、全逆 case と自由 frame の幾何を合わせた
独立受理である。gate 単独が数学の証明であるとはしない。
