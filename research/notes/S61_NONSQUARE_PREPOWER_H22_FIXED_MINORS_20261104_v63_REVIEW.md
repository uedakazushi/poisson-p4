# pre-POWER h22 の固定 pure3/pure4 小行列式：独立レビュー

2026-11-04。二つの一次証明書を一つの独立 checker で検査する。

- `S61_NONSQUARE_PREPOWER_H22_PURE3_FIXED_MINOR_20261104_v63`、
  report `56140959c7a6b5852f1da7dbaa2cbf9cb2306c6c59d84fb4c7fb13cc7e99b6a9`。
- `S61_NONSQUARE_PREPOWER_H22_PURE4_NAMED_STRIP_20261104_v63`、
  report `164ed080bee0c89f86c5fe7a1ad6d2533e45d46aa2e4d547da964e72223f2b3d`。

## 実 source binding

使用する段階は受理済み **pre-POWER15×22** だけである。三行は current15
の1,3,6、original25の10,13,16、actual135ラベル18,39,68に対応する。
固定 pure targets は original32 の11,13,25、Q9の q=22 を保持する。
H-adjugate や新しい B 分母を導入していない。

108個の保存済み scalar DAG 値を、元 checkpoint・既存 pair-dual・GL3 unit
atoms から別表現で検証する。さらに cached Q9 全九係数373項をその元の
pre-POWER row0 roots に直接照合するため、同じ凍結 DAG の87 ancestorを
追加して195 nodesを処理する。POWERの次数207 cofactorsは読まず、再評価しない。

Q9には選択した九係数以外にも八つの pure tail があり、これらを捨てたり
零と置いたりしていない。元22列の全 row から、その八つの label/root を
lazy source binding として保持する。新しい八係数の評価はしない。
選んだ三行と Q9 の残存 heavy5係数だけが親 DAG の構造的零である。
受理済み ZW25 と GL3 source/target unit change の参照を hash で固定し、
新しい source query は行わない。

## 独立算術

著者の Laurent-minimum-a normalization を再利用せず、独立 native は
\(P/a^k\), \(k\geq0\) を保持し、分母を揃えて比較する。
多項式 storage も著者の LEX ではなく DEGREVLEX とする。
全 atomic division が非零 QQ scalar と a の冪だけであることを確認し、
九成分および determinant の保存済み QQ integer-height envelope も確認する。

\[
M_{ij}=h_{22}P_i[j]-P_i[22]Q9[j]
\]

の九個の実係数を固定 source roots から生成して照合する。
三個の2×2 cofactor、第一行の cofactor expansion と六置換 expansion を
独立に検証し、非零1895項・全次数40の固定 determinant に一致することを
確認する。別の列・小行列式を探していない。

## 四行小行列式と strip

普遍多項式環で h=Q9[q] を一時的に反転すると、第一行 Q9 で第一列を
消す操作から

\[
\det(M)=h^2\det([Q9;P_1;P_2;P_3]_{[22,11,13,25]})
\]

が得られる。両辺が多項式なので、これは反転前にも成立する恒等式である。
保存された818項の pure4 quotient と h22² を独立に掛け戻し、実 determinant
との全係数一致を確認する。h22を全 FULLU の unit としていない。

既知 factor a,d,F,B,C,J のみの全九段の quotient/remainder と再乗算を検査する。
五個の不割り切れ remainder については、factor の LEX leading monomial で
割れる項が無いことも確認し、保存された単一多項式除算の refusal を検証する。
新しい一般 gcd・factorization・saturation は行わない。

結果は

\[
\det(\mathrm{actual\ pure4})=a^4C^3R,
\qquad R\ne0,\quad \#\operatorname{terms}(R)=434,\quad\deg R=20.
\]

R は未知の非定数因子のまま保持する。D(R) でこの四行 block の独立性を
得るが、R の零点集合を full matrix の defect locus と呼んではいけない。
全域 pure4 unit、全九核の独立性、全 FULLU cubic57、quartic239 は未証明。
既存 h22/他 chart の補集合も除外していない。

review code は
`code/s61_nonsquare_prePOWER_h22_fixed_minors_review_20261104_v63.py/.cpp`。
実測結果・全入力 hash・算術 audit は同名 `_REVIEW` report に保存する。
著者のコード・凍結入力は変更せず、producerを再実行しない。

native 実行前の checker preflight で、reviewer が Q9の列集合外13係数を
全て零と仮定した誤りを検出した（0.218秒）。上記八 tail を正しく保持する
よう修正し、失敗時の checker/hash/frontier を保存した。数学の native
計算に到達しておらず、小行列計算の再実行はしていない。
