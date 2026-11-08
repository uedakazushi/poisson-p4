# rank一 long・g零の七 whole 源：正確な25次元像と全既知・相互分離

2026-11-08、PRIMARY E_multiple_hyperplane、別作者 quartic_sources の純幾何証明。
入力は `RANK_ONE_TRANSLATION_WHOLE_FINITE_EXHAUSTION.md` の B1–B6 と V1。
新しい点、Poisson operator、Jacobi 行列、rank／minor、kernel、jet 又は CAS は作成・実行しない。
以下の結論は**七つの相異なる25次元像閉包と既知族からの分離**であり、
全空間成分性又は成分数追加ではない。

## 1. 元の全 source と七つの係数支持

\(X=\partial_x+x\partial_y\)、\(r=y-x^2/2\)、\(q=(z,w)\) を固定する。
旧 full source `MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §5 において \(g=0\) とすると
\[
\begin{array}{lll}
 P_{xy}=-F-(1+d)r,&P_{xz}=u,&P_{xw}=v,\\
 P_{yz}=xu+H_w,&P_{yw}=xv-H_z,&P_{zw}=0,
\end{array}\quad
 U=u\partial_z+v\partial_w,\ d=\operatorname{div}U,
\]
\[
 (U+1)H=0,\qquad \deg U\le1,\quad\deg H\le3,\quad F\text{ は全二次以下六係数}.
\]
全 Jacobi／divergence／Pfaffian は元式で保たれ、\(\operatorname{Pf}P=-H\)。
\(U,F\) と \(H\) の全偏微分は tensor から復元される。
潜在的な定数差 \(c\) も \((U+1)c=c=0\) で消えるので係数 map は injective。
\(F\) の自由度を Pfaffian から見えないという理由で除いていない。

全 affine frame／中心を戻した後、七つの開 source は次表で定める。
\(A,B\) はともに非零。

| ID | \(U_1\) の固有値 \((\lambda,\mu)\) | 中心座標での \(H\) |
|---|---|---|
| B1 | \((-1,-1/2)\) | \(Az+Bw^2\) |
| B2 | \((-1,-1/3)\) | \(Az+Bw^3\) |
| B3 | \((-1,1)\) | \(Az+Bz^2w\) |
| B4 | \((-2/3,-1/3)\) | \(Azw+Bw^3\) |
| B5 | \((-1/2,-1/3)\) | \(Az^2+Bw^3\) |
| B6 | \((-1/2,-1/4)\) | \(Az^2+Bzw^2\) |
| V1 | \((-1,0)\) | \(z(Aw^2+Bw+C)\)、\(A(B^2-4AC)\ne0\) |

前六行は degree 三以下の全 monomial eigenvalue 方程式 \(i\lambda+j\mu=-1\) の
ちょうど二つの解を使う。二 support の零補集合は、元全 routing の M1/M2/M11/M3/M21 へ戻す。
V1 の leading／discriminant 零も同じ元 routing に保持される。
本稿でこれらの別枝を七開の中へ混ぜない。

## 2. 七つとも全固定X像は既約12次元

B1–B6 では固有値は固定された相異なる非零値。
\(U_1\) の全 GL2 conjugacy orbit は既約二次元、全 translation は任意二次元で中心は一意。
その上の \(-1\)-eigenspace は表の二 monomial line の和で、両係数 unit の開は既約二次元。
全 \(F\in\mathbb A^6\) を戻して、各固定X source は正確に \(4+2+6=12\) 次元。
係数逆があるため像も正確に12次元である。

**V1 の中心は U 単独では一意ではない。**
\(U_1\) の conjugacy orbit は二次元、translation はその image の一次元だけなので、
\(U\) 自体の affine 族は三次元である。
各 \(U\) の全 eigenpolynomial は \(z(Aw^2+Bw+C)\) の三次元開であるから、
直接 \(3+3+6=12\) と数えられる。

同じ数え方を「中心二＋linear orbit二＋H二係数」にする場合は、
\(H\) が定める二つの root の**中点**を zero-eigendirection の中心に使う。
\(U\) の zero line と二つの simple root line の交点の中点は affine 幾何的に一意で、
\((U,H)\) から復元出来る。
その中心で \(H=z(Aw^2+\widetilde C)\)、\(A\widetilde C\ne0\) となる。
このように \(H\)-marked center を使えば全基底四次元と二 unit 係数の表示が injective である。
二つの root の順序選択は高々有限であり、余分な連続 fiber を生まない。
V1 source も既約12次元。U の zero-eigenvalue translation を誤って二重計数していない。

## 3. 全 stabilizer と scalar、唯一重複Hを戻した正確な25次元

既独立 `RANK_ONE_LONG_RESONANT_BINARY_CUBIC_ENTRANCE_SEPARATE_REVIEW.md` §§1,4 の
全 StabX 十次元・全 coefficient inverse を使う。
実 stabilizer は
\(x'=x+b,y'=y+bx+\ell(q)+c,q'=Bq+t\)、\(\det B\ne0\)。
\(\Delta=\ell(q)+c-b^2/2\) として、g零の全 parameter は pullback を省略して
\[
 U'=B_*U,\quad H'=\det(B)H,\quad
 F'=F-U(\Delta)-(1+d)\Delta.
\]
七つの source は全域でこの変換とその逆により保たれる。
V1 の中点も affine 変換で移るので、source の中心表示による欠落は無い。
\(a_\lambda:(x,y,q)\mapsto(\lambda x,\lambda^2y,q)\) と tensor scalar \(1/\lambda\) は
\(U'=U,H'=\lambda H,F'=\lambda^2F\) なので各 source を保つ。

\(\operatorname{Aff}_4\times S\) の作用写像の同じ tensor を与える fiber は、
divergence と source の全不変性・係数逆により正確に StabX の一 coset。
従って affine 飽和像は既約 \(20+12-10=22\) 次元で、scalar に関し cone である。

次表は全 source の actual divisor を中心・固有 frame で記す。
各式に任意の元 target frame を戻す。

| ID | actual divisor（非零定数倍） | 唯一の重複 H |
|---|---|---|
| B1 | \(h^3(Ahz+Bw^2)\) | multiplicity 三 |
| B2 | \(h^2(Ah^2z+Bw^3)\) | multiplicity 二 |
| B3 | \(h^2z(Ah^2+Bzw)\) | multiplicity 二 |
| B4 | \(h^2w(Ahz+Bw^2)\) | multiplicity 二 |
| B5 | \(h^2(Ahz^2+Bw^3)\) | multiplicity 二 |
| B6 | \(h^2z(Ahz+Bw^2)\) | multiplicity 二 |
| V1 | \(h^2z(Aw^2+Bhw+Ch^2)\)、\(A\Delta\ne0\) | multiplicity 二 |

B1 の quadric は rank 三で絶対既約。
B2／B5 の cubic はそれぞれ z／h に関する primitive 一次式なので絶対既約。
B3／B4／B6 の quadric も rank 三で絶対既約。
V1 は h と異なる相異三 simple linear factor を持つ。
従って全表示開で h が唯一の重複因子である。

H を動かす四次元と scalar 一を戻し、marked H を忘れる一般 fiber が一点なので、
七つの PGL 飽和射影像閉包 \(Z_{B1},\ldots,Z_{B6},Z_V\) は全て**正確に既約25次元**。
一点の orbit／tangent dimension をこの数え方に使用しない。

## 4. 七族の相互識別

全 actual D は essential-variable 数三以下という closed determinantal 条件に入る。
以下の標識は unique repeated H と irreducible factorization から generic に復元されるので、
任意の PGL 同値に対する不変量である。
各 closed 条件は、該当 factor と Grassmannian の標識を含む proper incidence から全閉包へ延ばす。

**B1 対他六族。** B1 は triple linear factor を持つ閉条件内。
他の generic D の最大 multiplicity は二であり、同次元閉既約像が一致することはない。

**B2 対 B5。** どちらも \(2H+C\)、C は essential 三変数の cuspidal cubic。
B2 の C の全 singular support は唯一の plane \(\{h=w=0\}\) で H に含まれる。
B5 の全 singular support は唯一の plane \(\{z=w=0\}\) で、H に含まれない。
「H 内に C の singular P2 が存在する」は Grassmannian の proper incidence による閉条件。
B5 にはそのような plane が無いので一致しない。
P4 の cubic cone は vertex を常に持つため、単なる \(\operatorname{Sing}C\cap H\ne\varnothing\) を
この識別の代わりに使っていない。

**B3／B4／B6。** \(D=H^2LQ\)、Q は rank 三、L は H と異なる simple linear factor。
marked pair の restriction rank は
\[
 (\operatorname{rank}(Q|_H),\operatorname{rank}(Q|_L))
 = (2,1),\ (1,2),\ (1,1)
\]
の順である。quotient essential plane ではそれぞれ H/L が secant/tangent、tangent/secant、tangent/tangent。
\(Q|_H\) 又は \(Q|_L\) の rank 一以下は、全因子・hyperplane 標識を戻しても閉条件である。
generic に rank 二となる別族はその閉条件の外なので、三つは互いに異なる。
H と L は multiplicity で標識されているため交換して同一視出来ない。

**残る組合せ。** B2／B5 の irreducible cubic は、B3／B4／B6 の reducible residual cubic と異なる。
V1 は all-linear の閉条件内、前六族の generic D は irreducible quadric 又は cubic を持つ。
同次元性によりこれらの像閉包も一致しない。
以上で七者の相互相違を得る。

## 5. 全 prior129、非定数g新成分、三共鳴成分との分離

固定 prior129 projection とその Git blob の全照合は先行独立 receipt の入力を使う。
rank二六成分と rank四 dim二十四以下百十四成分は、rank 又は正確な25次元により七像を含まない。

CA／CB の二十八次元成分の全 actual D は all-linear 閉条件内。
従って前六族を含まない。V1 については、CA の三乗 linear factor、
CB の二つの二乗 linear factor という proper 閉 incidence が、
generic V1 の multiplicity \((2,1,1,1)\) と両立しない。

旧 nil31-DDH／half／third の25次元成分も同じ all-linear／反復 factor 条件で除外出来る。
positive NEW25 は \(3H+Q\)、generic Q rank五。
B2–B6／V1 は triple linear factor を持たないので入らない。
B1 が入れば同次元で等しいはずだが、B1 は quadric rank三以下の標識 proper 閉条件内、
positive NEW25 の generic Q rank五はその外なので等しくない。
残る T01／T39／T40 の25次元成分には reduced actual D の非空開があり、
七像の全 rank四部分は repeated H の閉条件内なので、同次元一致は起こらない。

既受理非定数g成分 \(Z_{nc}\) は25次元、generic D は \(H^2LQ\) で Q rank五。
B1 の triple H、B2／B5 の irreducible residual cubic はその generic 型と異なる。
B3／B4／B6 の全 Q rank三以下の closed incidence は generic \(Z_{nc}\) を含まず、
V1 の all-linear closed incidence も generic \(Z_{nc}\) を含まない。
同次元なので、七像は \(Z_{nc}\) に含まれない。

新 constant-g 二共鳴像 \(Z_1,Z_{1/3}\) と scalar-cubic 像 \(Z_{SC}\) の
25次元下界／相互分離は先行独立 pure receipt、成分性の算術は
`RANK_ONE_LONG_THREE_AMBIENT_SEPARATE_REVIEW.md` の別入力である。
七像との比較自体には、その componenthood は必要なく同じ25次元だけで十分。

\(Z_1,Z_{1/3}\) の generic actual D は essential 五変数。
これは新 anchor を選ぶことなく、旧全 source の係数から次の微分 minor で分かる。
中心固有 frame で
\(C=g_0h^2y-g_0hx^2/2+t z^2w+hH_2+h^2H_1+h^3H_0\) のとき、
\((C_y,C_x,C_z,C_w,C_h)\) と monomial \((h^2,hx,zw,z^2,hy)\) の minor は
\(-4g_0^3t^2\ne0\)。
\(H_3=t z^3\) のときは \(H_2\) の \(w^2\) 係数を c とし、
monomial \((h^2,hx,z^2,hw,hy)\) の minor が \(-12g_0^3tc\ne0\) となる非空開を使う。
\(h^2C\) についても、constant direction による零微分を mod h で見ればまず h方向が零、
次に C の零方向が零となるため essential 五を保つ。
七像の全 D は essential 三以下の閉条件内なので同次元一致は出来ない。

\(Z_{SC}\) は all-linear なので前六族と異なる。
V1 との比較では simple 三 hyperplane の coefficient span を使う。
SC では全三 factor が \(\operatorname{span}\{z-z_0h,w-w_0h\}\) 内、span 二以下の閉条件。
V1 では \(z,\ w-\rho_1h,\ w-\rho_2h\)、\(\rho_1\ne\rho_2\) なので span 三。
proper marked multiplication と rank 条件から、この generic 相違は全閉包の不一致を与える。

## 6. M1／M11 の仮入口及び全零境界との関係

M1 の全 D は \(h^4L\)、M11 の全 D は \(h^3L_1L_2\)。
両像閉包は all-linear の proper 閉条件内なので、前六族の generic nonlinear factor を含まない。
V1 は最大 multiplicity 二であり、M1 の四乗又は M11 の三乗 factor 閉条件の外。
よって M1／M11 の帰属・成分性が後でどう決まっても、七像全体をそれらの既存源へ含める矢印にはならない。
交点となる特別な係数境界を排除したとの意味ではない。

各 B の A 又は B 零、V1 の leading／discriminant 零、F の全零・特殊係数、
全 short counterpart は元有限全 routing に保持する。
七 source の係数開だけから、その別枝の全 membership 又は全分類完成を結論しない。

## 7. 最終受理範囲

七像は exact whole source と actual stabilizer／scalar／唯一 H により、相異なる既約25次元閉族である。
全 prior129、既受理非定数g25成分、新三共鳴像及び M1／M11 仮入口とは上の意味で分離される。
これは**成分性の証明ではない**。未発見の高次元 component に含まれる可能性はまだ残る。
次に一 fixed fullambient rank100 を使うなら、この25次元下界と合わせる最弱の局所矢印で十分だが、
本稿では新しい点もその operator も作らない。
