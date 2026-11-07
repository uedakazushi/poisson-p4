# quartic cone と高次数 cone：全 Euler 源による独立監査

2026-11-07、主義務 G。**PASS**。
S(0,4) の全中心 parameter を覆う既約被約 quintic 排除、及び全 integral 非退化 degree四
rank零支持の G 次元≤24への接続を確認した。
同じ証明は degree五以上の全 integral 非退化 cone の既約 quintic 枝も排除する。
Jacobi、conormal、Aut の追加計算はこの最弱証明に不要である。

## 1. 全 projective section の完全な vertex 分解

点 vertex を持つ任意の integral 非退化 cone surface S⊂P4 を考える。
base hyperplane の座標を x=(x0,…,x3)、vertex を c、base curve を C⊂P3 とする。
C は integral かつ非退化。全 global bivector は外積 Euler と H1(TP4)=0 により
homogeneous 二次係数の raw lift を持つ。

全係数を落とさず

\[
 \Pi=B_0(x)+cB_1(x)+c^2B_2
 +\partial_c\wedge\{V(x)+cW(x)+c^2v\}
\]

と書く。B0,V は二次、B1,W は一次、B2,v は定数。
E=E_x+c∂c とすると projective restriction の零性は E∧Π|S=0 と正確に同値であり、

\[
 E\wedge\Pi=E_x\wedge(B_0+cB_1+c^2B_2)
 +\partial_c\wedge\{cB_0+c^2B_1+c^3B_2
          -E_x\wedge V-cE_x\wedge W-c^2E_x\wedge v\}.
\]

各 x が C の affine cone を動き c は全複素数を動くため、c の係数比較が全点で出来る。
まず c³ は B2=0。
c² は B1=E_x∧v を C 上で与え、これは一次係数なので I_C(1)=0 から ambient 全体で成り立つ。
c は B0−E_x∧W の各二次係数が I_C(2) に入ることを与える。
実際の linear field W(x)+cv に対する Euler-vertical gauge
E∧(W+cv) を引くと、これらの base 部分と mixed c 項が**同時に**消える。
従って全源の代表は

\[
 \boxed{\Pi=\sum_j Q_j(x)A_j+\partial_c\wedge V(x)},\qquad
 E_x\wedge V|_C=0,
\]

となる。Qj は I_C(2) の basis、Aj は constant base skew bivector。
V は任意の quadratic vector field で、C 上で radial という最後の条件だけを持つ。
global な linear radial coefficient の存在を仮定せず、V を Q×constant vector に制限しない。
従って非CM・非linearly-normalな base/cone の extra section も全て残っている。
これは cohomology を q=0/ACM として縮める ansatz ではない。

## 2. 最弱の fixed-factor 結論

I_C(2)=0 なら全代表は ∂c∧V、従って projective rank≤2。
I_C(2)=CQ なら全代表は QA+∂c∧V であり、恒等的に

\[
 \Pi\wedge\Pi
 =Q\{Q A\wedge A+2A\wedge\partial_c\wedge V\}.
\]

従って signed Euler-Pfaffian D も Q を因子に持つ。
E∧(Π+E∧L)²=E∧Π² なので、これは選んだ Euler 代表に依存しない actual D の因子。
係数は全 V、全 A を残しており、Jacobi 前の恒等式である。
D が零又は二次因子付きになるので、rank4・既約被約 quintic の場合は空。
Q の rank3/4、可逆性、base の smoothness、conormal bracket はこの結論の追加仮定ではない。

## 3. S(0,4) の全中心への適用

minimal-degree cone Y=S(0,4)⊂P5 の vertex を v、base を rational normal quartic C4⊂P4 とする。
中心 p は Y の外にあるので、その base 成分は非零。
vertex shear c→c+ell(x) により中心の vertex 成分を零に出来、これは Y を保つ実際の線形変換。
従って base P4 の中心について全場合分けすればよい。

closed Sec(C4) の内側なら、Brodmann–Park Theorem6.2(a) の全中心 criterion により
image は ACM・degree四・codimension二、従って既済 CI22 枝へ入る。
外側なら、線形 subsystem が C4 の全 length二 scheme を分離するので、射影は閉埋込み。
image C⊂P3 は smooth rational quartic、O_C(1)=O_P1(4)、非退化であり、S はその cone。
secant/tangent の衝突を閉 secant 条件に含めている。

h0(O_C(2))=9 に対し ambient quadrics は10次元なので I_C(2)≠0。
独立な二 quadric は共通一次因子を持ち得ない：素性・非退化性に反するためである。
従って二つあれば degree4 CI curve を定め、C と degree が等しいので C 自身となるが、
CI(2,2) の arithmetic genus1 は smooth rational C と矛盾する。
ゆえに I_C(2)=CQ。§2 が全中心一括で既約 quintic 枝を排除する。

producer の unique Q が smooth という補足も、quadric cone の Cl=ZR、H=2R による
degree4 divisor=4R=2H、Cartier quadric section→CI22という標準計算と整合する。
ただし本監査の主証明はこの補足に依存せず、仮に rank3 Q を残しても同じ因子結論を与える。
binary-quartic stabilizer、3+1中心、二根2+2の rank4 generator、projective conormal 補正を含む
producer §§3–5 の強い A=0 導出は、この accepted G 結論の依存に含めない。
特に、可約 D を持つ一般の rank4 Poisson tensor まで不存在と読み替えない。

## 4. degree五・六を含む全高次数 cone の系

任意の integral 非退化 C⊂P3 で degC≥5 なら、独立な二 quadric があった場合、
そのどちらも reducible ではない（素 ideal が一次式の積を含めば C は一平面に入る）。
従って二 quadric は coprime、degree4 CI curve を定める。
degree≥5 の C を含むことは Bézout の degree 比較に反する。
よって dim I_C(2)≤1、§2 の二つの場合で全源を排除する。

これは `RANK_ZERO_HIGH_DEGREE_CONE_EXCLUSION.md` の全 source 証明を独立に確認したものである。
normal/non-normal、singular base、nonCM extra の有無に依らない。
未知の degree五・六 rank零支持のうち cone 部分は、既約 quintic の一般元から全て排除出来る。
non-cone 支持及び rank2 restriction は引き続き補集合として残る。

## 5. 任意の integral 非退化 degree四支持の網羅性

`DEGREE_FOUR_RANK_ZERO_SURFACE_EXHAUSTION.md` の分類入力を独立に直接照合した。
Brodmann–Park [arXiv:0808.0090v3](https://arxiv.org/pdf/0808.0090v3) の
Introduction 印刷1–2頁は、integral 非退化 almost-minimal-degree variety を
normal+linearly normal な normal del Pezzo と、一つ上の空間の minimal-degree variety の外点射影に分ける。
Remark6.1(A) 印刷16頁は linearly normal del Pezzo が maximal、すなわち ACM であることを明記する。
これは引用する既知分類定理であり、今回 Fujita の本全体の証明を再精読・再証明したという主張ではない。

surface degree4 in P4 は degree=codim+2 なので正確にこの scope に入る。
normal+linearly normal は ACM→CI22。
その他の全場合は degree4 minimal surface in P5 の外点射影。
同論文 Notation/Remark2.1 印刷4頁の minimal-degree 分類から、元曲面は
Veronese、S13、S22、S04 の四つだけである。
より高次元の Veronese cone を surface の場合に混ぜない。

|元の全場合|既済又は今回の矢印|
|---|---|
|ACM image|Hilbert function(1,2,1)→CI22→全 PGL 次元≤24|
|Veronese、中心が closed secant 内|一次 Remark6.3(A)→ACM→CI22|
|Veronese、closed secant 外|既済 projected-Veronese 全源 rank4 排除|
|S13/S22|全中心は ACM 枝又は独立監査済み四 orbit、後者は全 D が Q 因子|
|S04、base 中心が closed secant 内|一次 Theorem6.2(a)→ACM→CI22|
|S04、closed secant 外|§1–3 の全 Euler 源と fixed Q 因子|

全 CI22 は regular pencil と all-member-degenerate pencil の二つの受理済み独立証明が覆い、
全 rank零・既約 quintic 源の PGL 飽和次元≤24を与える。
その他の分岐は既約 quintic の場合に空。
従って任意の integral 非退化 degree四 rank零支持について G≤24 が従う。
これは既知成分への全 E 帰属を追加して要求するものではない。

## 6. 保存と範囲

今回の新しい核心は全 vertex-weight の手証明と exterior 恒等式であり、新たな CAS/点計算は不要だった。
結果は `reports/CLASSIFICATION_PROJECTED_QUARTIC_CONE_REVIEW_20261107.json` に、
producer 三ノートの確定 snapshot/hash と既済各証拠への依存を記録する。
四モデルの全150/350算術、PV、両 CI22 の既済証明書を理由なく再実行しない。
共有台帳、Git、原稿、producer の確定 bytes は編集していない。
