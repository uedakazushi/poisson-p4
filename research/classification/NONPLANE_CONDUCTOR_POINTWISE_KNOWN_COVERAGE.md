# 非平面の導手支持を持つ既約非正規 quintic：全次数・全 rank の点ごとの被覆

2026-11-08、PRIMARY E。受理済み全源定理の純集約。
新しい算術・operator・kernel・profile・minor・GB・jet はない。
G の像次元上界を E の包含へ読み替えず、各次数の最終 pointwise 排除又は actual family を使う。
新成分増分は零。主稿・global 台帳・Git の登録は integrator の別作業である。

**定理。** 複素射影四次元空間上の一般 rank 四 Poisson tensor \(\pi\) の実退化因子
\(D\) が被約既約な非正規 quintic とする。
正規化商 \(\nu_*O_X/O_D\) の integral 二次元支持成分のうち、少なくとも一つが plane でなければ、
\[
 [\pi]\in T37\cup T38\cup T39\cup T40.
\]
非平面支持がただ一つであること、他の plane 支持が無いこと、clean、normal node、
導手の被約性、normalization-flat 又は generic \(\beta\) を仮定しない。

## 1. 元の全既約非正規入口

`IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE.md` §§1–2 を使う。
\(\nu:X\to D\)、\(L=\nu^*O_D(1)\)、正規化商の integral 二次元支持を \(S_i\)、
その一般点での正の length を \(\delta_i\) とすると、一般平面切断の正規化種数 \(g\) に対し
\[
 \sum_i\delta_i\deg S_i=6-g,\quad 0\le g\le5,\qquad
 L^3=5,\quad K_X+C_X\sim0,\quad L^2C_X=12-2g.
\]
正規化商は純二次元であり、nonnormal hypersurface の \(S_2\) 性と \(R_1\) の失敗から
少なくとも一つのこのような支持がある。
各 degree は一から六であるが、支持の全てを degree 六又は genus 零と仮定しない。

Hamiltonian 導分は標数零で正規化へ延長し、
\(\partial(a)b=\partial(ab)-a\partial(b)\) により導手を保つ。
従って各 reduced integral 支持は genuine Poisson 部分多様体であり、
一般 restriction rank は零又は二である。
rank 零だけを課したり、rank 二だけを課したり、単に coisotropic であることに置き換えたりしない。
また \(S_i\subset\operatorname{Sing}D\) なので actual equation は \(F\in I_{S_i}^{(2)}\) を満たす。

degree 一の integral surface は plane である。
本定理では、少なくとも一つ存在する**非平面** \(S\) を選ぶので degree は二から六となる。
これは all-plane の場合を排除したのではなく、本定理の対象外として残す操作である。

## 2. span が退化した曲面を先に戻す

\(S\) が hyperplane に入り、degree \(d\ge3\) なら
`IRREDUCIBLE_NONNORMAL_DOUBLE_SURFACE_DEGREE_RESTRICTIONS.md` §1 を使う。
prime CI ideal \(I=(h,q_d)\) は \(I^{(2)}=I^2\) を満たす。
degree 五では \(2d>5\) なので
\[
 F=h^2A_3+hq_dB_{4-d},
\]
従って \(h\mid F\) で既約性に反する。負次数の \(B\) は零と解釈する。
normality、restriction rank 又は特異点の一般性は不要である。

この議論は degree 二の quadric に流用しない。
degree 二はその全 smooth／cone の独立全源定理で扱う。
degree 一もこの次数比較で消していない。
従って degree 三以上では非退化 integral surface だけを以下の表へ入れればよい。

## 3. 全次数・全 restriction rank の最終表

| 非平面支持の degree | rank 零の全源 | rank 二の全源 |
|---|---|---|
| 二 | smooth quadric／quadric cone とも actual 既約 quintic では空 | smooth quadric／quadric cone とも空 |
| 三 | smooth cubic scroll は空、twisted-cubic cone は T38 | smooth cubic scroll は T37、twisted-cubic cone は T38 又は T39 |
| 四 | 全 integral 曲面で空 | 全 integral 曲面で空 |
| 五 | 全 genus 零／一 source は T40 又は空 | genus 零と elliptic genus 一は空、canonical index 二 node は T37、cusp 等は空 |
| 六 | 全 integral 曲面で空 | 全 integral 曲面で空 |

以下はこの表に必要な exact proof path と scope である。

### Degree 二：quadric の二 rank と全 base

smooth quadric には `SMOOTH_QUADRIC_CONORMAL_FACTOR_SEPARATE_REVIEW.md` の
全 multiplicity／Euler lift／conormal proof を使う。
quadric cone の rank 零は `QUADRIC_CONE_ZERO_ACTUAL_FACTOR.md` と
`QUADRIC_CONE_CI22_AND_TRIPLE_Y_ROOT_REVIEW.md`。
元 full source の semisimple／nilpotent 全 family、全 \(Y\)、零／単位補集合を戻し、
actual invariant hyperplane から proper factor を得る。

cone rank 二は最終の
`QUADRIC_CONE_RANK_TWO_FULL_BASE_REDUCTION_V2.md`、
`QUADRIC_CONE_RANK_TWO_IRREDUCIBLE_QUINTIC_EXCLUSION_V2.md` と
`QUADRIC_CONE_RANK_TWO_EXCLUSION_V2_SEPARATE_REVIEW.md` を使う。
元 full 三十三係数に含まれる追加 base 方向、全 \(F_u\)、\(C=0\)、SS／nil、
全 trace／pole wall と零／単位補集合を保持した後の pole／valuation 排除である。
旧六係数入口のまま述べた V1 whole 結論は撤回済みであり、本表へ使わない。
修正履歴は `QUADRIC_CONE_RANK_TWO_ENTRANCE_ERRATUM.md` に残る。

### Degree 三：四つの全 source

`DEGREE_THREE_CONDUCTOR_POINTWISE_KNOWN_COVERAGE.md` を使う。
surface 分類、二 restriction rank、smooth rank 零の正しい二-source V2、
smooth rank 二の whole A25 residue 比較と genuine seed の actual T37 曲線、
cone の全 rank 零 T38 と rank 二 T38／T39 がそれぞれ独立確認されている。
旧 scroll 四-source 全被覆の反例・撤回を迂回しており、seed の帰属だけから未知 higher を推論しない。

### Degree 四：全 integral surface、両 rank

`DEGREE_FOUR_CONDUCTOR_IRREDUCIBLE_QUINTIC_ALL_EMPTY.md` を使う。
almost-minimal-degree の全 integral 分類により、normal／nonnormal CI22、
非 CI の rank 四 quadric 上 quartic、projected Veronese、
S04 の secant 外全連続 centre を戻す。
rank 零には `CI22_NIL32_AND_DEGREE_FOUR_RANKZERO_ROOT_REVIEW.md`、
rank 二には integral CI22 全排除、quadric Weil-class 因子、
`PROJECTED_VERONESE_TANGENT_SL2_SEPARATE_REVIEW.md` の全十五源／Borel 排除を使う。
old projected-Veronese G 上界だけを pointwise 空へ読み替えない。

### Degree 五：全 normalization と actual family

length 式により \(\delta_S=1\)、\(g=0\) 又は一であり、
\(g=0\) のときは残る length 一の plane も保持する。
元全還元は rank ごとに独立した旧 source 定理を使い、以下の最終 E の矢印を付ける。

rank 零・genus 零は `GENUS_ZERO_RANKZERO_ACTUAL_T40_CLOSURE.md` と
`GENUS_ZERO_RANKZERO_QUOT_CONDUCTOR_SEPARATE_REVIEW.md`。
113／122 全 projection 五空間と target frame、plane 導手の二点又は全 first-jet 貼合せを保持する。
指定された quotient \(O^5\to E\) の relative Quot 変形と MCM 部分導手の base change により、
中心が元 tensor そのものである ambient Poisson 曲線を作る。
bundle だけの閉包、又は未証明の descent／正規化 flatness を仮定しない。

rank 零・genus 一は `RANK_ZERO_QUINTIC_SURFACE_GENUS_ONE_SEPARATE_REVIEW.md` の
全 canonical／cusp／zero-summand／nef 非 ample 排除と stable elliptic T40 認識を使う。
残る全 unstable \(A_2\oplus F_{2,3}\) は
`UNSTABLE_ELLIPTIC_RANKZERO_ACTUAL_T40_CLOSURE.md` と
`UNSTABLE_ELLIPTIC_RANKZERO_T40_SEPARATE_REVIEW.md`。
逆 extension、完全五 section の任意 frame、finite image と MCM conductor を同時に変形し、
元 actual tensor を保つ T40 への曲線を得る。

rank 二・genus 零は `QUINTIC_RANKTWO_GENUSZERO_WHOLE_EMPTY_SEPARATE_REVIEW.md`。
113／023／014、全 marks・collision・small resolution、node／cusp／non-split の旧入口を戻す。
rank 零となる残余 plane で actual inverse form が強制する非零二-residue wedge と、
同じ dualizing generator の反 deck 性を両立出来ず、全空である。
genus 一の elliptic 全源は
`RANK_TWO_QUINTIC_SURFACE_GENUS_ONE_ELLIPTIC_SEPARATE_REVIEW.md` で空。
canonical Gorenstein index 二の全十一係数・非 terminal・nonreduced vertex curve 等を含む源は、
`CANONICAL_INDEX_TWO_NODE_WHOLE_ORBIT_ENTRANCE.md` で全零／従属補集合を戻す。
node は `CANONICAL_INDEX_TWO_NODE_WHOLE_T37_COVERAGE.md` と
`CANONICAL_INDEX_TWO_NODE_ROOT_SEPARATE_REVIEW.md` の actual 全百五十係数同定により T37、
cusp は `RANK_TWO_QUINTIC_CANONICAL_CUSP_EMPTY_SEPARATE_REVIEW.md` で空である。
古い index 二 G 上界だけに依存した包含ではない。

### Degree 六：両 rank の全排除

rank 零は `RANK_ZERO_SEXTIC_CONDUCTOR_EXCLUSION.md` と
`RANK_ZERO_SEXTIC_CONDUCTOR_SEPARATE_REVIEW.md`。
全 normal 正規化・sectional genus 零定理から 113／122／014／023／005 を尽くし、
cone、cusp parity と node での Poisson／導手の必要条件を用いる。
支持 \(S\) 自身の normality／smoothness や cone 性は仮定しない。

rank 二は `RANK_TWO_SEXTIC_TWO_CURVE_RESIDUE_EXCLUSION.md` と
`RANK_TWO_SEXTIC_TWO_CURVE_SEPARATE_REVIEW.md`。
同じ五 scroll model、generic cusp／node、single／split prime、全 finite projection を戻す。
閉一形式の非零定数 residues は一般二 hyperplane 曲線と一般 ruling line の両方で留数和零を満たす。
split 二 prime の ruling 次数の和は三なので二つの次数は異なり、二関係は独立となる。
codimension 二の bad locus と small exceptional locus を曲線から避けるだけで十分で、
全域 log-smoothness 又は clean は仮定しない。旧 023 の proper locus G 上界もこれで全空へ強化される。

## 4. 結論の量化と残す範囲

任意の元 \(\pi\) に対し、仮定を満たす一つの非平面支持 \(S\) を選ぶ。
§1–2 から表のいずれかへ入る。空の枝は actual 被約既約 quintic 仮定と矛盾し、
残る枝は全 actual source の \(T37,T38,T39,T40\) への包含である。
各既知成分は閉かつ PGL 不変なので、元の座標変換を全て戻した \(\pi\) についても結論する。
支持が複数ある場合にも、選んだ一成分についてこの矢印があれば十分である。

本定理は **all-plane conductor**、normal な \(D\)、被約可約又は非被約 actual divisor、
rank 二以下の \(\pi\) を分類しない。それらを本表の空枝から除去したと解釈してはならない。
既知閉成分の境界を含む pointwise membership であり、全点が generic member と同じ
導手型・軌道・scheme 構造を持つとは述べない。分類全体の完了、総成分上界又は新 componenthood も主張しない。

同名 manifest は、この全源連鎖の proof／独立 receipt／既存公開 manifest の exact path と SHA を固定する。
新しい数学計算を再実行した記録ではなく、旧 report の埋込実行 bytes と current compatible code を同一視もしない。
