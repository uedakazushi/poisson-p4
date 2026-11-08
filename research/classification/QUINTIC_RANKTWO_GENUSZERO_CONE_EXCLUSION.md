# 次数五 rank 二・種数零：全023／014の排除

2026-11-08、PRIMARY E。純幾何 producer proof、独立レビュー前。
新しい算術、GB、minor、jet 又は選択点は無い。

**定理。** `RANK_TWO_QUINTIC_SURFACE_GENUS_ZERO_SEPARATE_REVIEW.md` の
残る023／014全源は、actual 被約既約 quintic を持つ rank 四 Poisson tensor を与えない。
`QUINTIC_RANKTWO_GENUSZERO_113_EXCLUSION.md` の全113排除と合わせ、
rank 二 degree 五・正規化 sectional genus 零の全分岐は空となる。
small exceptional curve、全射影、全 marks とその特別な交差を保持する。

## 1. 全源の共通記法と局所必要条件

旧 producer／別作者 review の全五 scroll・全 conductor class 還元を使う。
残る smooth scroll resolution は
\[
 Y=\mathbb P_{\mathbb P^1}(E),\qquad
 E=O\oplus O(2)\oplus O(3)\quad\hbox{又は}\quad O\oplus O(1)\oplus O(4).
\]
small map \(Y\to X\) は zero-summand section だけを点に縮め、exceptional divisor は無い。
\(X\to D\) は finite birational、\(D\) は実際の reduced irreducible quintic。
degree 五 rank 二支持 \(S\) の conductor は
\(C_S=\Gamma_1+\Gamma_2\)、\(\Gamma_i=(g_i)\in|L|\) は二 prime。
残余 plane \(P\) は generic node、\(\pi|_P=0\)、closed inverse-residue form は
\[
 \alpha=\lambda\,d\log(g_1/g_2),\qquad \lambda\ne0.                \tag{1}
\]
これらは旧全 source の結論で、normal な \(S\) 又は全域 involution を仮定しない。

113 note §2 の局所 inverse 計算を共通入力とする。
actual \(D=(uv=0)\)、\(\pi|_P=0\) なら、二 branch の regular restrictions
\(\theta_+,\theta_-\) は
\[
 \theta_+\wedge\theta_-\ne0                                  \tag{2}
\]
を満たす。理由は closed inverse が logarithmic、\(\sigma^2\) の log-volume 係数が actual Pfaffian の unit の逆、
その係数が \(\pi|_P=0\) の下では \(\pm2\theta_+\wedge\theta_-\) となることである。
全ての以下の矛盾は (2) への矛盾で、normal form 点の計算ではない。

## 2. 023 の全 quadric-cone marks

023では \(C_P=Q=\mathbb P(O\oplus O(2))\simeq\mathbb F_2\)。
\(Q\) の負 section を \(C_0\)、ruling fiber を \(f\) とすると
\[
 L_Q=C_0+2f,\qquad K_Q=-2L_Q.
\]
\(|L_Q|\) は \(C_0\) を vertex に縮める quadric-cone 埋込みであり、
\(Q\to P\) はその rank 三 quadric cone を cone 外の点から射影する degree 二 map。
有限性は \(X\to D\) に対して使い、\(Q\) の負 section を有限に写すと誤主張しない。
その deck involution \(\iota\) は cone の vertex と ruling を保ち、minimal resolution \(Q\) に延長する。
したがって ruling base \(\mathbb P^1\) 上の automorphism を誘導する。

\(\Omega_Y\) を同じ quintic dualizing generator の rational pullback とする。
small exceptional curve は divisor を追加せず、
\(\operatorname{div}\Omega_Y=-Q-\Gamma_1-\Gamma_2\)。
\(\beta=\operatorname{Res}_Q\Omega_Y\) は
\[
 \operatorname{div}\beta=-(\Gamma_1+\Gamma_2)|_Q,
 \qquad\iota^*\beta=-\beta.                                  \tag{3}
\]
後者は generic node の同じ dualizing generator の trace-zero residue。
ambient unit 自身の deck 不変性は置かない。
従って二 marked divisors の和は \(\iota\)-不変である。

各 \(g_i|_Q\) は非零 \(L_Q\)-section で、vertex coordinate \(x\) と high coordinate \(y\) を使い
\[
 a_i x+B_i(s,t)y,\qquad a_i\in\mathbb C,quad \deg B_i=2           \tag{4}
\]
と書ける。
非零性は prime \(\Gamma_i\) が異なる class の \(Q\) を含まないことによる。
全型は次の三場合に分かれる。

* \(a_1a_2\ne0\)：各 divisor は vertex を避ける irreducible plane conic の strict transform。
  (3) より deck は二つの irreducible divisors を固定又は交換する。
  従って \(\iota^*d\log(g_1/g_2)=\pm d\log(g_1/g_2)\)。
* 一方だけ \(a_i=0\)：その divisor は \(C_0\) と二 ruling fibers の和であり、
  他方の唯一の horizontal irreducible conic は (3) により固定される。
  残りの divisor も固定されるので、ratio の deck pullback は定数倍の同じ ratio。
* \(a_1=a_2=0\)：共通の \(y\) を消すと ratio は \(B_1(s,t)/B_2(s,t)\)。
  元の一形式も deck pullback も ruling base 一変数だけの一形式であり、その wedge は零。

repeated binary roots、共通 fibers、比例した restrictions、二 conic の一致もこの三場合に含む。
projective \(Q\) 上で同じ divisor を持つ二 section の比は定数なので、上の divisor 議論からの ratio 結論に余分な unit は無い。
全場合で
\[
 \alpha_Q\wedge\iota^*\alpha_Q=0,
\]
となり、(2) に反する。従って023の全 source は空。
smooth quadric の factor-swap 補題を singular cone に流用していない。

## 3. 014 の二 conductor sheets と同じ第三 line

014では
\[
 C_P=Q+F_b,\quad Q=\mathbb P(O\oplus O(1))\simeq\mathbb F_1,
 \quad Q\sim L-4F.
\]
\(Q\) の負 section を \(C_0\)、ruling fiber を \(f\) とすれば
\(L_Q=C_0+f=H\)。\(Q\to P\) は \(|H|\) の blowdown、\(F_b\to P\) は線形同型である。
どちらも generic に degree 一であり、二つの node sheets を与える。

重要なのは、\(Q\cap F_b\) が \(Y\) 上の**同じ一つの curve**であること。
その image は両方の maps の下で同じ line \(\ell\subset P\) となる。
\(Q\) 上では ruling fiber、\(F_b\) 上では一つの line で、双方に \(L\)-degree 一を持つ。
二つの独立な第三 line を選ぶ freedom は無い。

\(g_i\) の各 restriction は非零。
\(Q\) に恒等的零なら \(\Gamma_i\) が \(Q\) を、\(F_b\) に恒等的零なら \(F_b\) を含み、prime class \(L\) に反する。
したがって \(\Gamma_i|_Q\) は plane line \(A_i\subset P\) の total transform、
\(\Gamma_i|_{F_b}\) は plane line \(B_i\subset P\) である。
\(A_i\) が blowdown point を通る場合には total transform が \(C_0\) を含むが、
plane 上の divisorial image は同じ一つの line \(A_i\) のままである。

\(\operatorname{div}\Omega_Y=-Q-F_b-\Gamma_1-\Gamma_2\) から、二つの residues の divisors は
\[
 \operatorname{div}(\operatorname{Res}_Q\Omega_Y)
       =-(F_b+\Gamma_1+\Gamma_2)|_Q,
 \quad
 \operatorname{div}(\operatorname{Res}_{F_b}\Omega_Y)
       =-(Q+\Gamma_1+\Gamma_2)|_{F_b}.                           \tag{5}
\]
これらを同じ plane \(P\) の rational 二形式と見れば、generic node の trace-zero 関係により互いに反対。
従って their plane divisors は等しく、(5) は
\[
 A_1+A_2+\ell=B_1+B_2+\ell,
 \qquad\hbox{従って}\qquad A_1+A_2=B_1+B_2                      \tag{6}
\]
を与える。
blowdown の exceptional curve は plane 上 codimension 二に写るので、rational 二形式の plane divisor に別の成分を加えない。
反復 line や \(A_i=\ell\)、\(B_i=\ell\) の場合も (6) では multiplicity を保つ。

plane polynomial 環の UFD から、unordered pair \(\{A_1,A_2\}\) と \(\{B_1,B_2\}\) は一致する。
従って \(Q\) と \(F_b\) から来る (1) の二 regular restrictions は
\[
 \theta_Q=\lambda d\log(A_1/A_2),\qquad
 \theta_F=\lambda d\log(B_1/B_2)=\pm\theta_Q.
\]
比例 line の場合も両者の wedge は零である。
これは (2) に反し、014の全 source も空となる。
この証明は plane matching を generic Grassmann 閉包から仮定せず、
**同じ dualizing generator と同じ交差 curve の image**で直接処理する。

## 4. 全 degree 五・種数零 rank 二の aggregate

旧全五型の補集合は次の通りである。

| 全正規化／導手型 | disposition |
|---|---|
|005|vertex line による二次元非正規支持は plane のみ、degree 五 \(S\) 不可：既受理|
|122及び全型の \(C_P=2F\)|split \(C_S\) の prime class が存在しない：既受理|
|014、\(C_P=2Q\)|二 prime の正 \(L\)-係数の和に反する：既受理|
|generic cusp \(S\)|五 cusp／Gauss degree bound：既受理|
|non-split node \(S\)|closed residue と trace-zero 非零性の矛盾：既受理|
|113、\(C_P=Q\)|smooth quadric の全 bilinear-factor lemma と actual-node wedge：113 note|
|023、\(C_P=Q\)|本稿第2節の全 cone mark 補題|
|014、\(C_P=Q+F_b\)|本稿第3節の共通第三 line と全二 line 比|

実 projection の intersection 三／四 Schubert 補集合は113 note §6により全型で空。
exact-intersection 二の全13次元 source と全 marked sections は上記全場合で処理した。
したがって rank 二 degree 五・正規化 sectional genus 零の全 actual 被約既約 quintic 枝は空。

genus 一 elliptic／canonical index 二、他次数、他 actual factor 型、全 E と総 component 数は別義務である。
新 component 増分は零。旧 source／reports／global／主稿／Git は変更しない。
同名 manifest に入力を固定し、独立純幾何レビュー前の producer theorem として保存する。
