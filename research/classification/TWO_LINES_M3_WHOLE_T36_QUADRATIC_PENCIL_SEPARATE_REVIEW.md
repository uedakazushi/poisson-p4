# node-flex／cusp-flex全六磁気源のT36非所属：独立純査読

2026-11-08、PRIMARY E、別作者 `/root/log_manuscript`。
`TWO_LINES_M3_WHOLE_T36_QUADRATIC_PENCIL_EXCLUSION.md` の§7を含む全域強化を独立に手導出した。
**判定PASS：第一・第二源は全六磁気、正規化J零を含む全点でwhole T36の外。**
新operator、point、minor、profile、局所変形jetの延長又はCASは無い。
陽性T37親、新成分性、T38の例外lambda又は全E完成を追加しない。

## 1. 既受理の全源normalizerとproperなT36条件

全域actual shearの全逆は既別作者
`TWO_LINES_PLANE_CUBIC_THREE_QUOTIENTS_SEPARATE_REVIEW.md` を入力にする。
第一・第二の全R6は、source parameterを反転せずR=Ih²+Jhyへ戻る。
本査読ではJを一にせず、全I,Jを残す。
全U/V/dはproducer§1の式をそのまま使用する。

whole T36については独立に既確認したproper marked-zero条件を使う。
order二以上のzero qのquadratic tensor Qに対する
\(\alpha_Q=\iota_{E\wedge Q}\mathrm{vol}\) が、
Gr(2,Sym²)のproper像 \([G_2dF_2-F_2dG_2]\) の閉cone C22に入る点qが少なくとも一つある。
全completed T36 sourceのq+/q-式、零／一次jet incidence及びP4 propernessによる
全成分への閉拡張は
`GA9_LINEAR_W_WHOLE_T36_QUADRATIC_PENCIL_EXCLUSION.md` §§1–3の既定理である。
固定frameだけの必要条件へ弱めていない。

## 2. 全I,Jでorder二zeroの候補をvertex lineに尽くす

h一のanchor共通零は第一(0,0),(1,0)、第二(-1,0),(0,1),(0,-1)。
第一(0,0)と第二(-1,0)はw=-Iの全fiber line。
Pxzのx微分は-2なので、その全点はsimpleである。

第一(1,0)にはI零の時だけ追加zero planeがあり、Pxzのx微分は+2。
第二(0,±1)にはI±J零の時だけ追加planeがあり、同じ微分は+2。
全I=J零も含め、追加plane全点をsimpleと直接確認した。
generic陰性をそのまま例外parameterへ代入した議論ではない。

h零、x非零のraw Pfは非零定数倍x四乗なので、rank四raw tensorはEulerとの外積にはならずprojective zeroでない。
h=x零、y非零では、raw tensorがそれぞれ∓3y² partial_y∧partial_wであり、
projective zeroは正確にz零のline。
y一chart、line上w=w0で、projective coefficientは

\[
 P_{zw}^{\rm aff}=P_{zw}-zP_{yw}-wP_{zy}.
\]

その全first jetは第一Jh+2w0x+3z、第二Jh+2w0x-3z。
したがって両者ともz微分が非零で、全I,J、全w0でsimpleである。
このlineのy零端点はvertex line \(L=\{h=x=y=0\}\) の一部として別に保持する。
以上より全order二以上zero候補は厳密にLだけである。

## 3. 両fiber chartと全Lの二次jet

q=[1:kappa]のz一chart、r=w-kappaで比の括弧を使うと、
base-base partは \(E_b\wedge U\) というdegree三項、mixed partは

\[
 V-(\kappa+r)U-E_b[-d(\kappa+r)+R].
\]

よって全quadratic partは

\[
 Q_q=(V-\kappa U+\kappa dE_b)\wedge\partial_r.
\]

無限遠q=[0:1]のw一chartでは、base-base partはdegree三 \(E_b\wedge V\)、
全quadratic partは
\(Q_\infty=(U-dE_b)\wedge\partial_r\)。
Rの全係数はE_bRによりdegree三以上に入り、二次jetには出ない。

従って全q=[zeta:omega]でalphaは、baseからpullbackした
\(\iota_{E_b\wedge(\zeta V-\omega U)}\mathrm{vol}_b\) の非零scalar倍。
Saito determinant \(\det(E_b,U,V)\ne0\) から、全qでalphaは非零。

producer§4の「base同士の括弧は零」は、全projective affine tensorではなく
**そのquadratic partが零**と読むのが正確である。
実際には上記degree三項が残る。このscope補正は二次jet式と結論を変えず、旧bytesは変更しない。

## 4. 全qでR22を除くdegree三の論証

第一の \(C=h y^2-x(x-h)^2\)、第二の \(C=h y^2-(x+h)^3\) は絶対既約。
h一のquadratic-in-y式の右辺は、それぞれx(x-1)²、(x+1)³であり有理関数体上のsquareでない。
全anchorに対する手微分は

\[
 \begin{array}{c|cc}
 &U(C)&V(C)\\\hline
 \text{第一}&2(3x-h)C&-6yC\\
 \text{第二}&6xC&6yC
 \end{array}
\]

を与える。従って全qでcubic cone C=0はalphaのfoliationのinvariant hypersurface。
alphaはcubicであるが、その全係数がCに割られるならalpha=C*constant one-formとなり、
Euler水平性からconstant一形式が零となる。これは§3の非零性に反する。
よってC一般点でfoliationは非零で、Cの接平面そのものである。

独立quadrics F2,G2のpencilなら、その有理第一積分はC一般点で定義され、C上で定数となる。
degree二の分母はirreducible cubicを含められず、char零よりF2-cG2がCに割られる。
degree二対三から零polynomialとなり独立性に反する。
したがって全q∈LでalphaはC22外。

§2は全I,Jで他のorder二zeroを除き、§3–4はL全体を除く。
whole T36の必要なproper marked-zero条件を満たさないため、
**第一・第二の全R6 sourceは点ごとにT36外**である。
全actual normalizer逆、source frame零entry、J零及び全parameter例外を保持した結論である。

本査読で新しい有限計算を再実行していない。
別のT38 zero-surface障害と併用しても、T38例外parameterやT37親認識は別義務である。
