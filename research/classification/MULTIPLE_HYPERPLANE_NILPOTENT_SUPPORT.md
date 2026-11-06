# 重複超平面：rank3 冪零モジュラー場の全支持源と非優越

2026-11-07。課題 E_multiple_hyperplane と G。
入口は `REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` §3–5。
このノートは rank3 nilpotent Jordan4 の全 source を尽くす。
この型の primitive β 像の次元は高々13である。
s=1 の全源は既存22次元成分へ認識した。
残る四源は相異なる四つの新22次元成分であり、全周囲の算術・全源包含・相違を独立確認した。

## 1. 全モジュラー標準形

`P=C+L+Q` は四変数の二次以下の真の Poisson bivector で、
最高部は `Pf(Q)=0` を満たす。div は `div(P)^j=Σ_i ∂_i P^{ij}` とする。
モジュラー線形部が rank3 の冪零行列なら Jordan 型は (4) に限られる。
線形座標を選び

\[
 N=\begin{pmatrix}0&1&0&0\\0&0&1&0\\0&0&0&1\\0&0&0&0\end{pmatrix},
 \quad X=Nx=x_2\partial_1+x_3\partial_2+x_4\partial_3
\]

とする。平行移動でモジュラー定数の im N 成分を消し、
`Z=X+γ∂4`、γ=0 または1とする。γ≠0 を1とする変換は四座標の homothety である。
従って Z の零点の有無を落としていない。

rank3 の Koszul 表示（`MULTIPLE_HYPERPLANE_MODULAR_SUPPORT.md` 補題1）により
`Q=X∧Bx` と書けるが、この冪零型では B と N の可換性を仮定できない。
div Q=X は

\[
 BN-NB=sN,\qquad s=1+\operatorname{tr}B
\]

を与える。N の中心化環が `span(I,N,N²,N³)` であることから、全 B は

\[
 B=D_s+rN^2+tN^3,\qquad
 D_s=\frac14\operatorname{diag}(7s-1,3s-1,-s-1,-5s-1)
 \tag{1}
\]

とできる。N の倍数を B から除く操作は Q を変えない。
(1) は最高部の全三変数源であり、`[Q,Q]=Pf(Q)=0` と div Q=X を恒等的に満たす。

次に全 C6、L24 を導入して、`[Z,P]=0` と `div P=Z` の全係数を比較する。
これは (s,r,t,C6,L24) の33未知数に関する有理係数の affine 線形方程式である。
γ=0 では32行・階数24、γ=1では33行・階数26となる。
以下の源はその全解であり、残る Jacobi を別に全て課す。

## 2. モジュラー零点を持たない場合は rank4 解がない

γ=1 の全線形解を七変数 `(a,b,u,v,w,e,f)` で書くと

\[
 s=-1/9,\quad r=(u-w)/2,\quad t=-v,
 \qquad(C^{12},C^{13},C^{14},C^{23},C^{24},C^{34})
 =(a,e,-b+f,b,v,w)
\]

であり、L の六成分は同じ順序で

\[
\begin{split}
 (&(u+w)x_1-vx_2+fx_3+ex_4,\quad
   (3u+w)x_2/2+fx_4,\\
  &4x_1/9+wx_3+vx_4,\quad ux_3+vx_4,\quad
   x_2/3+wx_4,\quad 2x_3/9).
\end{split} \tag{2}
\]

全四つの Jacobi 多項式の係数は、ちょうど次の四条件と同値である。

\[
 u+w=0,\quad 5b-9uw-9w^2=0,\quad
 9v(u+w)-4e=0,\quad
 7a+9b(u-w)-9f(u+w)=0.
\]

従って全解は `u=−w`、`a=b=e=0`。
これを全六括弧へ代入すると `Pf(P)=0` が恒等的に成り立つ。
よってモジュラー零点を持たない γ=1 の枝には一般 rank4 の点はない。

## 3. 零点を持つ全 rank4 支持は五つの四変数源に限られる

γ=0 の全線形解は九変数 `(a,b,c,d,e,f,r,s,t)` により、(1) と

\[
 C=a\partial_{12}+b(\partial_{23}-\partial_{14}),
\]

および次の六成分の L で尽くされる。

\[
 \left(\frac c3x_2+2dx_1+ex_4+fx_3,\quad
 \frac c3x_3+2dx_2+fx_4,\quad
 -\frac{2c}3x_4+dx_3,\quad cx_4+dx_3,\quad dx_4,\quad0\right).
 \tag{3}
\]

全 Jacobi 係数は、次の六式と同値である。

\[
\begin{split}
 &d(3s-1)=0,\qquad c(7s-1)=0,\qquad b(s-1)+4d^2=0,\\
 &8cr+12dt+(1-15s)e=0,\\
 &(15s-3)a-12br+8c^2+12df=0,\qquad bc=0.
\end{split} \tag{4}
\]

`s∉{1/3,1/7,1,1/5,1/15}` なら (4) は a=b=c=d=e=0 を強制し、
残る全 f,r,t に対して Pf(P)=0。従って rank4 の支持は次の五つで尽くされる。

| s | 全解の代入 | 自由な四変数 | rank4 条件 |
|---|---|---|---|
| 1/3 | `c=0,b=6d²,e=3dt,a=36d²r−6df` | d,f,r,t | d≠0 |
| 1/7 | `d=b=0,e=7cr,a=28c²/3` | c,f,r,t | c≠0 |
| 1 | `c=d=e=0,a=br` | b,f,r,t | b≠0 |
| 1/5 | `b=c=d=e=0` | a,f,r,t | a≠0 |
| 1/15 | `a=b=c=d=0` | e,f,r,t | e≠0 |

各行は四次元既約 affine 源の一つの principal open であり、
全ての値で Jacobi が恒等零となる。
表の非零条件は必要十分である。実際の affine Pfaffian は順に

\[
\begin{split}
 &-\frac d3(108d^3+54drx_4^2+18dx_2x_4-9dx_3^2
 -9fx_4^2+3tx_4^3+3x_1x_4^2-3x_2x_3x_4+x_3^3),\\
 &-\frac c3x_4(14cx_4+7rx_4^2+2x_2x_4-x_3^2),\\
 &-\frac b2(2b+rx_4^2+2x_2x_4-x_3^2),\qquad
 -\frac a2x_4^2,\qquad -\frac e3x_4^3.
\end{split} \tag{5}
\]

したがって零係数の境界まで保持し、単一の標本で rank4 を判定していない。
最高項・全低次項・実 Pfaffian のいずれも落としていない。

## 4. 固定 symplectic 基点の primitive β 像は次元高々13

**命題。** 重複超平面を無限遠へ取った代表のモジュラー線形部が
rank3 nilpotent Jordan4 となる全点について、固定した symplectic 基点と値 J における
primitive β 像の次元は高々13である。従ってこの全型は W=A16 を優越しない。

**証明。** γ=1 の rank4 支持は第2節で空。
γ=0 の標準場 Z=Nx の affine 変換群 `Aff4=GL4⋉A4` における安定化群は、
N の可逆中心化群（次元4）と ker N の平行移動（次元1）からなる。
従って場 Z の軌道は `20−5=15` 次元である。
この場を固定した Poisson source は、第3節の五つの四次元源で全て覆われる。
モジュラー場への写像は Aff4 同変なので、その全 affine 支持の次元は高々 `15+4=19`。
これは source の中心化群作用を余分に引いた計算ではなく、場の軌道上の各 fiber の上界である。

基点0で P(0) が非退化となる開を取り、評価写像 `P↦P(0)` を考える。
GL4 はこの支持を保ち、非退化な定数 bivector の六次元空間に推移的に作用する。
従って全評価 fiber は同型であり、固定値 J の fiber の次元は高々 `19−6=13`。
primitive β はこの fiber 上で一次部 L の線形像だから、その像次元も高々13である。

元の重複超平面を無限遠へ送る perspective は固定基点と J を保ち、primitive β を変えない。
したがって入口ノート命題2の全代表への橋を通じて、求める β 非優越が得られる。□

この主張は固定 J の支持上界であり、全空間の Poisson 成分の次元や既知成分への所属を
上界13と呼び替えていない。表の全五源の成分帰属は以下で独立に証明する。
その認識を、G の非優越を成立させる追加前提にはしない。

## 5. s=1 の全源は既存の三重超平面・二次曲面成分へ入る

`legacy/v05/triple_hyperplane.tex` の既存点は、座標 `(h,u1,u2,u3,u4)` で

\[
 P_{\rm old}=h^2J+H\wedge N_0,\quad
 J=\partial_{u_1u_2}+\partial_{u_3u_4},
\]

\[
 H=3u_1\partial_{u_1}-3u_2\partial_{u_2}+u_3\partial_{u_3}-u_4\partial_{u_4},
 \quad N_0=u_4\partial_{u_2}-u_1\partial_{u_3}+u_3\partial_{u_4}.
\]

線形座標変換

\[
 (x_1,x_2,x_3,x_4)=(u_2,2u_4,4u_3,-8u_1)
\]

により、`2N0→N`、`−H/2→D1`、`J→8∂14−8∂23` となる。
従って既存点は、本ノートの s=1 源の `b=−8,f=r=t=0` に厳密に一致する。
これは退化因子だけの一致ではなく、全 bivector の線形変換による一致である。

旧定理はこの点で全350×150 Jacobi 微分の階数103、
vertical・scale・PGL 軌道方向の階数47を証明しており、
この点を通る全空間成分は唯一で、射影次元22、既存軌道閉包
`R4_TRIPLE_HYPERPLANE_QUADRIC` である。
s=1 の全 rank4 源 `(b≠0,f,r,t)` は既約で、この既存の滑らかな点を含む。
従ってその像閉包を含む全空間成分は上記唯一の成分であり、**全源**がこの既存成分に入る。
この認識と全源包含の論理は別担当が独立に読解確認した。

## 6. 残る四源の成分性と区別

各源で指定された非零パラメータを1、他の三パラメータを0とした点を取る。
四変数の各係数を h で次数2へ斉次化した全五変数 bivector に対して、
元の150係数と全350本の射影 Jacobi 条件を使う。
exact QQ 計算は次を与える。

| s | 全 Jacobi 微分の階数 | vertical・scale・軌道の階数 | 実四源方向を加えた階数 |
|---|---:|---:|---:|
| 1/3 | 103 | 47 | 47 |
| 1/7 | 103 | 46 | 47 |
| 1/5 | 103 | 47 | 47 |
| 1/15 | 103 | 47 | 47 |

各行で射影全接空間の次元は `150−103−24−1=22`。
一方、PGL と全四源の既約な実像の微分は `47−24−1=22` 次元である。
よって局所次元と embedding dimension が22で一致し、その点は滑らかで、
全空間成分は唯一となる。
全四変数源の PGL 飽和の閉包は既約で、この成分に含まれ、次元も22なので成分そのものに等しい。
従って各行は、選んだ点だけでなく**全源の成分帰属**を与える。

s=1/3,1/5,1/15 の三つの代表は射影軌道次元22なので、各成分の稠密開軌道を与える。
s=1/7 の代表の軌道次元は21である。
この一点の軌道次元から一般軌道次元や有理商の次元を推論しない。

### 6.1 実退化因子による区別

四成分の一般退化因子はそれぞれ

\[
 2H+C_3,\qquad 2H+H'+Q_3,\qquad 3H+2H',\qquad 2H+3H'
 \tag{6}
\]

である。C3 は既約三次、Q3 は rank3 の既約二次、H≠H'。
第1型の三次式は x1 について一次で、その係数は非零定数倍の x4²。
x4=0 に制限した残りの式は `x3³−9dh x3²+108d³h³` で零ではないため、
係数の gcd は1である。Gauss の補題により三次式は既約。
第2型の二次式は、(5) の斉次化で
`14ch x4+7r x4²+2x2x4−x3²` となり rank3、従って既約である。

非被約五次因子を持つ条件は、射影的な積写像の像の有限和として閉じている。
各源の全点がこの条件を満たすので、その PGL 飽和閉包である各成分全体も満たす。
従って SNC 点を持つ MPS の40成分のいずれとも一致しない。
一般 rank2 の成分とも異なる。

既存の同次元の非MPS成分である滑らかな二重二次曲面型、三重超平面・rank4二次曲面型、
rank3二重二次曲面の混合型とは、(6) の因子次数・重複度・二次因子の rank が異なる。
他の既知成分とは成分の次元で区別できる。
第1・第2型は相互にも最後の二型にも因子次数・重複度で区別される。

### 6.2 同じ因子重複度を持つ最後の二成分の区別

最後の二型は超平面を交換すれば同じ因子重複度なので、(6) だけでは区別しない。
全五変数で `E=h∂h+Σx_i∂i` とし、斉次 lift の canonical divergence-free representative を

\[
 \Pi=P-\frac15E\wedge N
\]

とする。この lift は線形座標変換に同変である。
第3・第4源では b=c=d=0 なので、全 f,r,t に対し

\[
 \Pi^\sharp dh=-\frac15hN,\qquad
 \Pi^\sharp dx_4=\frac{25s+1}{20}x_4N
 \tag{7}
\]

となる。ここでも第一添字での縮約規約を用いる。
実退化因子の**唯一の三重超平面**の一次式を ℓ、その cofactor を
`X_ℓ=Π♯dℓ/ℓ` とする。この組は射影座標変換で同変で、Π の定数倍は X_ℓ を定数倍するだけである。
従って `ℓ∈im(X_ℓ^*:V*→V*)` という条件は射影同型不変である。

s=1/5 では ℓ=h、`X_ℓ=−N/5`。
一次式上の N の像は `span(x2,x3,x4)` なので h はその像に入らない。
s=1/15 では ℓ=x4、`X_ℓ=2N/15`。
この場合 x4 は `im((N*)³)=span(x4)` に入る。
従って二つの代表は射影同型ではない。

両代表の軌道は、それぞれの22次元成分の稠密開軌道である。
同じ既約成分に異なる二つの稠密開軌道は存在できないから、二成分は相異なる。
これにより、四源は既知成分とは別の四つの22次元成分を与える。

この節の全周囲係数計算は独立実装でも再構成された。
全源包含の論理、canonical cofactor と固有の三重因子による最後の二区別は、
独立担当および統合担当が読解確認した。四成分の下界台帳への反映は統合担当が行う。

## 7. 正確検証と独立確認の射程

`code/classification_multiple_hyperplane_nilpotent.py` は全33未知数を導入し、
二つの γ chart の全線形方程式から(2),(3)を導く。
その後、全四 Jacobi 多項式、全六残条件、五つの全パラメータ代入、
実 Pfaffian (5) と各非零係数を QQ 上で検証する。
`jobs/multiple_hyperplane_nilpotent.json` の有限契約で実行し、実測約1.8秒。

最初の線形還元だけの証明書
`reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_NILPOTENT_20261107.json` は不変に保持した。
全支持源の最終証明書は
`reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_NILPOTENT_SUPPORT_20261107.json` で、
両証明書とも対応する producer の全 bytes と hash を埋め込む。
前者の hash と後者の更新コードの hash を混同しない。
独立 checker の
`reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_NILPOTENT_REVIEW_20261107.json` は、
全 C6/L24 を含む別実装の150×33線形系の階数24/26、全源、全 Jacobi・Pfaffian、
γ=1 の rank4 空性、五源の rank4 条件を約0.065秒で確認した。
行の重複・零係数の保持方式による150行と一次 producer の32/33行の違いは解集合を変えない。
第4節の場の軌道15・固定Jでの次元減少6・primitive β の橋は統合担当が独立に読解確認した。

第6節の有限な全周囲判定は `code/classification_nilpotent_source_components.py` と
`jobs/nilpotent_source_components.json` により実行した。
`reports/CLASSIFICATION_NILPOTENT_SOURCE_COMPONENTS_20261107.json` は、全四代表、
103次非零 minor と全 kernel basis、実 source・軌道方向の消滅と47次非零 minor、
producer bytes・入力 hash を保持する。実測約7.1秒。
`reports/CLASSIFICATION_NILPOTENT_SOURCE_COMPONENTS_REVIEW_20261107.json` は
一次 producer の helper を import せず、全四源の普遍 Jacobi、全Dの階数103、
vertical の階数24、全source方向を含む階数47、全55方向のDによる消滅を約6.8秒で独立確認した。
s=1/7 の軌道方向だけが階数46となる点も一致した。
