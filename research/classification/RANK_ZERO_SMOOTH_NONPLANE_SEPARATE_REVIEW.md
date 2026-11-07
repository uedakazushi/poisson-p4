# smooth非plane rank零枝：全quadric sourceと統合の独立確認

2026-11-07。producer は root、reviewer は plane_collision_review。
`RANK_ZERO_SMOOTH_NONPLANE_SUPPORT.md` を独立に読んだ。
**判定：記載された全smooth非plane枝のG還元としてPASS。**
既受理の個別証明を入力とする幾何学的統合であり、CAS・個別sourceの重複計算は行っていない。

## 1. smooth quadricの全sectionを落としていない

Sのspanは一つのP3、I_S=(ell,q)、qはその四変数で非退化と取れる。
E=Λ²TP4についてKoszul列は

`0→E(-3)→E(-1)⊕E(-2)→E I_S→0`。

通常Eulerと外積EulerからE(-3)はacyclicである。
例えばT(-3)とΛ²C5⊗O(-1)がともにacyclicなので直ちに従う。
またH0(E(-2))=Λ²C5、E(-1)の全sectionは全linear-coefficient bivector Lでliftできる。
従って全sourceは `Π=ell L+q A` で表せる。
Lの五十lift係数に五つのEuler冗長性があるが、これを条件として源を縮めていない。
定数skew Aの十係数も全て含む。

smooth非planeでrank零なので、generic normal multiplicity二とFrobenius first jetから
abelian conormalが必要となる。
prime CIのI²はI-primaryで、generic symbolic-square条件を全I²へ移せる。
local projective generatorsell/h,q/h²への変更の補正も、Π自体がIにあるためI²に入る。
従って実三次多項式 `{ell,q}` はI²に属する。

degree三のI²の項はell²linear及びell q constantだけであり、q²の項は無い。
ell=0への制限は `q*v(q)=0`、vはAのell-rowが与えるconstant base vectorとなる。
四変数qの非退化性からv=0。従ってΠのell-rowは全てellで割れ、hyperplane HがPoisson invariant。
H上ではπ∈Λ²TH、dimH=3なのでπ²=0。ゆえにellが実Dを割り、既約被約quinticに反する。

ここで普通二重点、clean blowup、vertical divisor不存在は一度も仮定しない。
quadric coneではvがradicalに残り得るため、同じ排除を主張しない。

## 2. degree boundと個別入力の適用範囲

degree≤6はclean仮定を要しない。
一般planeによるDの切断はintegral平面quinticで算術種数六。
Sとのdegree個の交点はそれぞれsingular pointでdelta≥1なのでdegreeS≤6である。

S⊂P3なら、degree≥3はprime CI squareの次数論でDにそのhyperplane因子を生む。
degree二は§1、degree一は本稿から除外したplaneである。

非退化smooth Sでq>0の場合、既受理の小次数irregular曲面補題はsmoothnessとq>0だけを使う。
そのdouble-point/ruled-surface証明を読み直し、clean residueからq>0を得る前段を
ここへ仮定として持ち込んでいないことを確認した。
結論はdegree五の楕円scroll。全stable bundleでのresolution/Bottとq52存在から
H0(E I_S)が一次元となり、今回の非零πはT40に属する。

q=0、degree≥5は既受理の全Euler源補題とlinear-normalityにより空。
二quadricを含めない次数論から源は零又はq2Aだけで、実Dは零又はq2²linearとなる。

degree三のsmooth Sは唯一のcubic scroll。
その完全abelian-conormal源と独立ideal-lift/leading-cover証明は、全PGL像dimension≤23を与える。
この入力を全E帰属又はrank四不存在へ強めていない。

## 3. degree四のscheme-theoreticな場合分け

linearly normalでないsmooth integral非退化surface in P4は既受理Severi補題の
projected Veroneseだけであり、全rank零Poisson源排除の独立証明を適用できる。

linearly normalかつq=0なら全源はI_S(2)⊗Λ²C5。
r=h0(I_S(2))≤1は既約Dを与えない。
r≥2では二つの独立quadricを選ぶ。素idealと非退化性により各quadricは既約、互いに素なのでCIを作る。
そのdegree四CIがdegree四のSを含むため、唯一のminimal componentはSでgeneric multiplicityは一。
CIのunmixed性によりembedded成分も無く、schemeとしてSと一致する。
cycleの次数だけをscheme同一視の理由にはしていない。

したがって全smooth CI22の既受理parity/Vandermonde証明が適用され、rankは二以下。
ここでは、新しい全singular-pencil CI22定理を先取りしなくてもsmooth枝を閉じられる。

## 4. 統合結果と保持する例外

以上で、smooth非plane rank零曲面を持つ今回のπはT40又はdimension≤23のcubic-scroll像に入る。
従って既知T40以外の未知25〜27次元G成分の一般元からこの全smooth枝を排除できる。
node/非node、vertical有無による補集合は、この結論には残らない。

Sの特異点、非normal性、quadric cone、rank零plane、rank二restrictionは残る。
全非正規quinticのG/E又は成分分類全体の完了を結論しない。

確認記録は `reports/CLASSIFICATION_RANK_ZERO_SMOOTH_NONPLANE_REVIEW_20261107.json`。
各採用済み個別定理の算術を今回再実行したという意味ではなく、全仮定と矢印の独立監査である。
