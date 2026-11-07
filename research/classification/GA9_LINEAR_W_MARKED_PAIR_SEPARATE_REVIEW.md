# Ga9 linear-w と新103成分の marked-pair 非包含：独立監査

2026-11-07。producer `quartic_sources`、reviewer `log_manuscript`。
`GA9_LINEAR_W_MARKED_PAIR_OBSTRUCTION.md`
（SHA256 `63d93d812eedb7acb7b7b37d9247aa17ec5a607ca444e035fe613b08def79320`）を独立に吟味し、
以下の **marked-pair 閉包非包含**として受理する。
新成分又は全 E 被覆は主張しない。
監査契約は `jobs/ga9_linear_w_marked_pair_review.json`、結果は
`reports/CLASSIFICATION_GA9_LINEAR_W_MARKED_PAIR_REVIEW_20261107.json` である。

## 1. 正確な対象と結論

generic pair は rational integral plane quartic F と line H で、F は A3 tacnode と A2 cusp を一つずつ持ち、
H は tacnode の二 branch に各 order 二で接する。全交点重複度四をこの一点で使い切る。
既受理の103番目の一般退化因子がこの pair であることを入力とする。
target pair は
\[
 H_0=(h=0),\qquad F_0=4h^3w-8h^2y^2+y^4.
\]
F0 は w-linear primitive（\(h^3\) と定数項の gcd 一）なので絶対既約。
\(\partial_wF_0=4h^3\) から唯一の singular point は \([0:0:1]\)。
例えば正規化
\[
 [h:y:w]=[4s^4:4s^3t:8s^2t^2-t^4]
\]
は basepoint-free、generic birational であり、そこで vanishing sequence は \((4,3,0)\)。
局所 leading 型 \(h^3+y^4\) の E6 cusp と、unique inverse point も確認出来る。

受理するのは target **marked pair** が generic pair locus の PGL-stable closure に無いこと。
unmarked quartic の adjacency、103成分と Ga9 磁気零 boundary の交わり又は他親成分への帰属は否定しない。

## 2. 全 DVR 族を同じ平面へ移せる範囲

仮に closure に入れば、dense constructible marked-pair source と target を結ぶ curve を取り、
normalization・completion・有限 extension 後、\(R=\mathbb C[[t]]\) 上の族を得る。
factorization は projective factor-space の proper multiplication により延長する。
special の linear factor は H0 唯一、quartic は irreducible なので、factor の special value もこの pair である。
各係数を primitive に取り、product の scalar は R の unit で調整出来る。

Poisson の全 P4 族から来る場合にも、固定平面を追加仮定していない。
generic と special の D はともに essential-variable rank 三。
gradient coefficient matrix の全4minorは generic に零なので R で零、special のある3minorは unit。
従って R 上で rank 三、二次元 kernel と三次元 quotient は free direct summand である。
その正則 basis で全 D を三変数へ移せる。
H と F もその quotient のみに依存する。余分な変数について polynomial degree は product で加法的なので、
D の degree 零なら両 factor の degree も零だからである。
target frame は GL5(R) の unit determinant で選べ、generic PGL parameter に pole があった可能性もこの議論で吸収する。

primitive quartic は generic に geometrically integral、special に integral であり、
total \(\mathcal C\subset P^2_R\) も integral、flat である。
非normal total space を排除していない。

## 3. 同時 P1 正規化を仮定せず証明する

excellent R 上の integral surface \(\mathcal C\) の normalization
\(N\to\mathcal C\) は finite。
N は integral なので R-torsionfree、従って flat、かつ proper である。
normal surface は CM なので Cartier special fiber N0 は一次元 CM、embedded point を持たない。
\(\mathcal C_0\) の generic point で total \(\mathcal C\) は regular なので normalization は同型。
したがって N0 は generic reduced・generic degree 一である。
finite map は追加 curve を点へ縮めず、第二 component が全 \(\mathcal C_0\) を覆うことも generic 一対一性に反する。
ゆえに N0 は reduced integral curve。

generic rational quartic の正規化は geometrically P1。
有限 extension で P1 と split してから normalize してよい。
もし後で extension を追加する場合は normalization を取り直して同じ論証を用いる。
normalization が任意の ramified base change と可換だとは仮定しない。
proper flatness から \(\chi(O_{N_0})=\chi(O_{N_K})=1\)。
N0 は rational F0 と birational なので、その normalization は P1、
\(0\to O_{N_0}\to\mu_*O_{P^1}\to Q\to0\) より \(1=1-\operatorname{length}Q\)。
従って Q零、N0=P1 である。
generic/special の両 geometric fiber が smooth だから N/R は smooth。
Henselian complete R では三つの異なる special point を section へ lift して N=P1_R と trivialize 出来る。
\(\nu^*O(1)\) は fiber degree 四なので \(O_{P^1_R}(4)\) である。

この手証明は全 integral rational-quartic family を扱い、normal total、simultaneous normalization の flatness、
normalization の base-change 交換又は δ一定定理を別に仮定しない。

## 4. marked section と全 unit 座標の標準化

generic tacnode の二 inverse point と cusp inverse point は有限 extension 後に section となる。
singular section の極限は special singular point、N0 上の inverse は一つなので、三者は同じ source point へ近づく。
そこを affine coordinate 零とし、無限遠を異なる special point に取れば、全三 section は R の有限値を持つ。
midpoint による正則 translation で tacnode inverse を \(b,-b\)、\(b\in tR\setminus\{0\}\) と出来る。
cusp を \(\zeta\in tR\) と書き、generic では \(\zeta\ne\pm b\)。

marked flag（tacnode point on H）は R 上の flag であり、unit target basis により
point を \([0:0:1]\)、Hを h零へ正則に直せる。
H の pullback は degree 四を二点各二で使うため、\(h=\lambda(r^2-b^2)^2\)。
special h の零 divisor が \(4[0]\) であることから \(r^4\) 係数 \(\lambda\) は R-unit。
unit scaling で \(h=(r^2-b^2)^2\) とする。

y は二点で消え、monic \(q=r^2-b^2\) の polynomial division は R 上で行える。
generic に remainder 零なので、R 上でも
\(y=q(Ar^2+Br+C)\)、\(A,B,C\in R\)。
\(y-Ah=q(Br+C+Ab^2)\) は正則な target shear である。
special tangent h と異なる line y は E6 で order 三なので、Bは unit。
その逆で scaling すると \(y=q(r+c)\)、\(c\in tR\) となる。
最後の \(w\in R[r]\)、degree 四以下は \(w_0(0)\ne0\) を満たす。
全 source/target 座標変更は unit determinant で、tの pole 又は vanishing scalar の除算は無い。

## 5. generic cusp 条件の恒等式を R へ戻す

cusp の projective derivative は零で、\(h(\zeta)\ne0\) だから二式は
\(h'y-hy'=0\)、\(hw'-wh'=0\) である。
手で全式を展開すると
\[
 h'y-hy'=q^2(r^2+2cr+b^2).
\]
従って \(\zeta^2+2c\zeta+b^2=0\)、
\(q(\zeta)=2\zeta(\zeta+c)\)。
\(b\ne0\) in K なので \(\zeta\ne0\)、\(q(\zeta)\ne0\) も保持している。
第二式は \(q(\zeta)w'(\zeta)=4\zeta w(\zeta)\) となり、
**generic field 上でだけ**非零 q と \(2\zeta\) を消して
\[
 2w(\zeta)=(\zeta+c)w'(\zeta)
\]
を得る。ここからの両辺は R の元であるから、R の K への単射性により R の恒等式である。
specialize すれば右辺は零、左辺は \(2w_0(0)\ne0\) となり矛盾する。
generic の非unit を割った座標を special に延長したのではない。

## 6. Poisson closure への適用と負の結果の限界

103番目の既約24成分は、上の generic marked-pair source の PGL 飽和閉包。
actual degeneracy-divisor map は非零 Pfaffian open 上で regular、Ga9 canonical はその open にある。
もし canonical point が成分に属せば、dense constructible source に generic point を置く curve と
上の proper factorization により排除した DVR marked-pair family が得られる。
したがって Ga9 linear-w の canonical と、その dense23 orbit は同成分に属さない。
accepted whole11 source bound/orbit23 の橋により、全 linear-w source の同成分への包含も否定される。

これは全 source の **全点が外にある**という主張ではない。
磁気零 boundary などの交わり、他の24/25親成分、全 E membership 又は新104成分性は未決定である。
本 obstruction は一つの偽の親族認識 route を閉じる最弱の負の結果として保存する。
164頁 TeX/PDF、global 台帳、Git は今回変更しない。
