# 正規化 E4 の nilideal 二生成と完全関係：独立レビュー

2026-11-04。対象は正規化比較環だけであり、元 D0 の nilideal の生成数や
principality は扱わない。一次証明書は
`AMBIENT_B18_NORMALIZER_NIL_GENERATORS_20261104_v63`、report SHA-256
`9dccf13ca56ed2f5d8676680d9e023b093ebfa740de8bbaa1dcafc5fde7a1af3`、
evidence SHA-256
`3ec81afe22c28db5dea27eb1f8834b7383d8668fd3dea63d2acb65fa8b6fbc0e`。

## 独立算術

著者の SymPy builder を import・再実行せず、四変数 sparse Fraction 演算で
元 Q0,Q1,f、GL4 の両合成と determinant1/3、全 source factorization、
二つの元 nil generator、全六関係とその Q0/Q1/f cofactors を直接確認する。
源多項式は11項・27項・8項の実 core source に一致し、形式的に各五直線へ
制限して零となる。これは新しい source 点の評価ではない。

単位の扱いで重要なのは、Delta は**24項の非定数多項式**であって
Delta(0)=-2880ということである。逆を -1/2880 に置き換えていない。
保存行列 N,E,T について両側の literal polynomial identities

\[
N[E,T]=[E,T]N=\Delta I_4
\]

を確認する。したがって完備局所環で実際に
\(N^{-1}=[E,T]/\Delta\)。これが Q0,Q1 を除く unit elimination を保証する。

## radical と module の範囲

受理済み実 core の全五直線支持を依存 gate として使用する。
N,P は d=0 内、C,R,L は d を含まない。実 linear prime generators と
その独立性も照合し、C,R,L が GL4 でそれぞれ c,r,l 軸となることを確認する。

\[
c=(d+x)/2,\quad r=(d-x-2a)/6,\quad l=(d-x+a)/3,
\quad h=a-d+x-2t,
\quad I_3=(h,cr,cl,rl).
\]

三軸 prime の monomial intersection を lcm と divisibility だけで再構成した。
商は三軸の被約環で、\(d=c+r+l\) は各軸で非零の parameter だから regular。
従って、N/P の embedded components を minimal intersection に重複して
入れることなく、

\[
\sqrt{dJ}=(d)\cap(u_5+I_3)=d(u_5+I_3).
\]

ambient regular S9 で d は regular なので、multiplication by d により

\[
\operatorname{Nil}(E_4)
 =d(u_5+I_3)/d(u_5+J_{\rm core})\simeq I_3/J_{\rm core}
\]

が S9-module として得られる。u5 はこの nil module に作用すると零であり、
E4 自身から u5 を消去してはいない。全 d*u_j 関係を保持する。

## 完全 presentation と極小性

monomial ideal `(h,cr,cl,rl)` の六個の pair-lcm syzygies は Taylor first
presentation を与える。保存された五本の syzygy が六本全てを生成することを
literal に確認した。最後の pair は保存した後二本の差である。
多項式環での完全性は平坦な completion により core S4 にも保たれる。

この五本と、Q0,Q1,f の実三列を合わせれば I3/Jcore の完全4×8 presentation。
上記 N による可逆な generator change で Q0,Q1 の二列を除くと、残る
**完全2×6**が `E^t*[syzygies5,V_f]`（共通単位 Delta を掛けた形）になる。
六列の像が Q0,Q1,f の保存 cofactor 結合に等しいことも全係数で検査する。
full S9 または E4 上では、五 u_j と二 generator の十関係を追加した
**完全2×16**を用いる。六関係・十六関係の最小性は主張しない。

実 source class matrix は

\[
\begin{pmatrix}0&0&0\\32&0&0\\11&0&0\\-9&-90&0\end{pmatrix}.
\]

rank2で f の class は零。I3/mI3 の四生成元から独立な二 class が残り、
同時に全2×16 relation entries は極大 ideal 内なので、nil の最小生成数は
厳密に2である。元座標の生成元は

\[
n_1=d(a-d+x-2t),\qquad
n_2=d(d+x)(a+d-x)/6=dcl.
\]

どちらも ambient 極大 ideal の二乗内に入る。しかしこれは m*Nil に入る
という意味ではない。nil-cotangent image が零でも nilideal は非零・二生成で
あり、これを D0 の nil module に外挿することもできない。

review code は
`code/ambient_B18_normalizer_nil_generators_review_20261104_v63.py`。
独立算術・proof・固定 report は同名 `_REVIEW` reports directory に保存する。
新しい jet、source 点、prime、一般体、GB、係数探索は実行しない。
