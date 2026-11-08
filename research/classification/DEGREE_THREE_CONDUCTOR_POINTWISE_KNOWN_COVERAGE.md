# 次数三の非正規導手曲面：全 restriction rank の点ごとの既知成分被覆

2026-11-08、PRIMARY E。独立受理された全源定理の純粋な集約。
新しい計算、kernel、operator、標本、minor、GB 又は jet はない。
元の証明書・失敗記録・主稿・global 台帳・Git は変更しない。新成分増分は零。

**定理。** 複素射影四次元空間上の一般 rank 四 Poisson tensor \(\pi\) の実退化因子
\(D\) が被約既約 quintic とする。正規化商の二次元支持に integral degree 三曲面
\(S\) があれば、\([\pi]\) は既知閉成分 \(T37\cup T38\cup T39\) に入る。
他の導手曲面、特に追加の plane の存在を除かない。
\(S\) の normality、一般 normal-cone の node 型、clean、vertical 不存在、
\(R=\operatorname{div}(\pi|_S)\) の被約性又は normalization-flat は仮定しない。

## 1. 全幾何入口と二つの restriction rank

`IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE.md` §§1–2 により、Hamiltonian 導分は
正規化と導手を保つ。従って各 integral 支持 \(S\) は Poisson 部分多様体で、
その一般 restriction rank は零又は二である。coisotropic 条件だけで代用しない。
また actual quintic 方程式 \(F\) は \(I_S^{(2)}\) に入る。

\(S\) が hyperplane に入る場合、`IRREDUCIBLE_NONNORMAL_DOUBLE_SURFACE_DEGREE_RESTRICTIONS.md`
§1 の prime complete-intersection square 補題を使う。
\(I_S=(h,q_3)\)、\(I_S^{(2)}=I_S^2\) なので、degree 五では
\(F=h^2A_3+hq_3B_1\)。従って \(h\mid F\) となり空である。
この排除は restriction rank と \(S\) の特異点に依存しない。

残る非退化 integral degree 三曲面は minimal-degree 分類により、
smooth cubic scroll \(S(1,2)\) と twisted-cubic cone \(S(0,3)\) の二型だけである。
各型を固定した後も全 actual ambient 座標変換を戻す。

## 2. 四つの全 source と結論

| 曲面と一般 restriction rank | 全点での結論 | 最終全源証明と独立確認 |
|---|---|---|
| smooth \(S(1,2)\)、rank 零 | 被約既約 actual \(D\) では空 | `CUBIC_SCROLL_ZERO_TWO_WHOLE_SOURCES_V2.md`、`CUBIC_SCROLL_ZERO_ACTUAL_FACTOR_V2.md`、`CUBIC_SCROLL_SEGRE_SUPPORT_SEPARATE_REVIEW.md` |
| smooth \(S(1,2)\)、rank 二 | 全 \(T37\) | `CUBIC_SCROLL_RANK_TWO_POINTWISE_T37_RECOGNITION.md`、`CUBIC_SCROLL_A25_WHOLE_RESIDUE_SEPARATE_REVIEW.md`、`PLANE_M2_CONCURRENT_T37_SEPARATE_REVIEW.md` |
| cone \(S(0,3)\)、rank 零 | 全 \(T38\) 又は actual proper factor／rank 不足 | `CUBIC_CONE_ZERO_T38_WHOLE_RECOGNITION.md`、`CUBIC_CONE_ZERO_T38_SEPARATE_REVIEW.md` |
| cone \(S(0,3)\)、rank 二 | 全 \(T38\cup T39\) 又は actual proper factor／rank 不足 | `CUBIC_CONE_RANK_TWO_POINTWISE_T38_T39_COVERAGE.md`、`CUBIC_CONE_RANK_TWO_NODAL_FISHER_SEPARATE_REVIEW.md` |

**Smooth rank 零。** 最終 V2 は元二十四関係の被約支持を全 \(L_6\) と
一つの \(3\times3\) rank 一 Segre cone の合併へ両方向に戻す。
rank 一行列の全 outer-product 表示は零行・零列も含み、全元 tensor／Jacobi／actual Pfaffian を保持する。
二つの全像で actual quintic はそれぞれ proper factor を持つ。
旧「四 source が全支持を覆う」という証明は反例を伴って撤回済みであり、
この結論には使わない。旧四族が個別には有効であることと、全被覆とは区別する。

**Smooth rank 二。** 全 normal-root／directrix 還元と canonical-root 排除の後、
全 source は A25 node／cusp に入る。同じ flat-fiber constant を二つの divisor で比較し、
node は全空、cusp は全 actual normalizer と presentation の逆を戻して一つの実軌道又は空となる。
primitive denominator のみを使用し、非 primitive な第二 row、全 \(\Omega\) 重複、
零 parameter、全四 full-P2 special fiber、全 higher source を保持する。

一つの genuine cusp-scroll tensor の存在は別の exact source により確認されている。
その actual 五次は絶対既約、smooth scroll 上の intrinsic bracket は非零で、
normal quadratic は一般 rank 一である。
同 tensor からの全 polynomial Poisson 曲線は非零 parameter で一つの actual 軌道に入り、
その一点は真の completed T37 tensor と全百五十係数で正逆 PGL／Euler 同定される。
閉性から中心 tensor も T37 に入る。
全 A25 が同じ実軌道へ戻る定理と合成するので、一点から未知の higher source を外挿していない。

この最後の実在／actual T37 同定は、当初は独立 review 待ちであった。
**現在は** `PLANE_M2_CONCURRENT_T37_SEPARATE_REVIEW.md`
（SHA256 `1ceb7cf180db41ddb679108bba2a0c2d91cb23b2ff02c0b22051165dad4947f8`）で
genuine seed、全 polynomial 曲線、全正逆百五十係数と百二十六 Pfaffian が独立 PASS になっている。
従って本集約にはその待機条件を残さない。

**Cone rank 零。** 元十八係数全源を戻す。
semisimple linear base の全六 vertex 係数は actual shear と真の completed T38 source に一致する。
nilpotent linear base は rank 一定六の全 kernel bundle の稠密 semisimple 部分から閉性で T38 へ入る。
quadratic nilpotent base は全 actual factor、base 零は rank 不足である。
残 Jacobi 条件を省かず、generic 分母の nil fiber への代入で被覆を主張しない。

**Cone rank 二。** 元二十七係数からの全五 whole family、全十九 wall、無限遠、
base 零を戻す。A の \(\lambda=2/5\) は全六 vertex source ごと rank 零 cone source と一致して T38、
\(\lambda=3/5\) は全零／単位補集合と shear／Cartan の逆により一つの実代表又は proper factor になる。
その代表は nodal genus 一 Pfaffian model の polynomial Fisher tensor と literal に一致する。
smooth generic model からの実多項式特殊化なので T39 に入る。
対応する nil fiber も全 rank 六 kernel bundle の閉性で覆う。
他の全 A、全 B、W-only と base 零は actual proper factor／rank 不足である。
singular model へ smooth-only theorem を直接適用しておらず、
局所 tangent 次元又は G 上界から component membership を推論していない。

## 3. 全点への合成と exact 入力

任意の元 \(\pi\) について §1 で曲面型と rank を一つ選べば §2 の該当する全源定理が適用できる。
各正の結論は PGL 不変な既知閉成分への actual membership なので、全座標変換を戻せる。
factor／rank 不足の枝は冒頭の actual 被約既約 quintic 仮定と矛盾する。
従って冒頭の結論が従う。複数の導手成分のうちどれを選んだかによる一意性は不要である。

Exact proof／receipt の公開 manifest は次の通りであり、本集約の同名 manifest に現 bytes を固定する。

* `inputs/cubic_scroll_zero_two_sources_v2_manifest_20261108.json`
* `inputs/cubic_scroll_segre_support_separate_review_manifest_20261108.json`
* `inputs/cubic_scroll_a25_whole_residue_constant_manifest_20261108.json`
* `inputs/cubic_scroll_a25_whole_residue_separate_review_manifest_20261108.json`
* `inputs/plane_m2_concurrent_t37_separate_review_manifest_20261108.json`
* `inputs/cubic_cone_zero_t38_whole_recognition_manifest_20261108.json`
* `inputs/cubic_cone_zero_t38_separate_review_manifest_20261108.json`
* `inputs/cubic_cone_rank_two_pointwise_t38_t39_manifest_20261108.json`
* `inputs/cubic_cone_rank_two_nodal_fisher_separate_review_manifest_20261108.json`

歴史的実行 report の埋込 producer／archived job と、後の compatible current runner は別の identity である。
本集約は古い実行を current code で再実行したとも、既存 manifest の歴史 SHA が全て current path と一致するとも主張しない。
既知各成分の componenthood／相違を再証明するものではなく、他導手次数又は all-plane 支持をこの定理だけで閉じない。
