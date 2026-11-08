# rank一 long 三共鳴源：三つの25次元像と全既知分離の独立純監査

2026-11-08、別作者 quartic_sources、PRIMARY E_multiple_hyperplane。**純全源入口 PASS**。
対象は `RANK_ONE_LONG_RESONANT_AND_BINARY_CUBIC_ENTRANCE.md`
（SHA256 `7d5699db0a7d5ba249aab27b4ee1140aee618b9a7e566c984c3d2ff5bfef8369`）。
新しい点、Poisson operator、Jacobi 行列、kernel、minor、jet 又は CAS は作成・実行していない。
受理するのは三つの相異なる既約25次元 Poisson 像閉包と全 prior129／新非定数g25成分からの分離まで。
三像の全空間成分性・smoothness・新しい成分数は、ここでは未証明として保持する。

## 1. 元の全係数と復元

既全長鎖 source は \(X=\partial_x+x\partial_y\)、\(r=y-x^2/2\)、\(q=(z,w)\) として
\[
\begin{array}{lll}
 P_{xy}=-F-(1+d)r,&P_{xz}=u,&P_{xw}=v,\\
 P_{yz}=xu+rg_w+H_w,&P_{yw}=xv-rg_z-H_z,&P_{zw}=g,
\end{array}
\]
\(U=u\partial_z+v\partial_w\)、\(d=\operatorname{div}U\)、
\(U(g)=dg\)、\((U+1)H=Fg\) である。
\(\deg g\le1,\deg U\le1,\deg F\le2,\deg H\le3\)。
旧全源定理 `MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §5 の同値を入力として使い、
highest symbol だけの必要条件へ置き換えない。

補助座標 \((t,r,q)=(x,y-x^2/2,q)\) では
\(P=\partial_t\wedge Y+\operatorname{Jac}f\)、
\(Y=-(F+(1+d)r)\partial_r+U\)、\(f=rg+H\)。
従って \(\operatorname{div}Y=-1\)、\(Yf=-f\) で全 Jacobi が成立し、
元 volume を保つ逆変更後も \(\operatorname{Pf}P=-f\)、\(\operatorname{div}P=X\)。
非線形補助変更はこの恒等式確認だけに使い、PGL membership に使わない。

係数逆は全域である。\(P_{zw}\) から \(g\)、\(P_{xz},P_{xw}\) から \(U,d\)、
\(P_{xy}\) から \(F\)、他二式から \(H\) の両偏微分を得る。
二候補の差が定数 \(c\) なら \((U+1)c=c=0\) なので差も零。
g零、T零又は lower coefficient 零でも隠れた加法定数の fiber は生じない。

## 2. 全 constant-g 三次共鳴の二つの12次元 source

\(g=g_0\ne0\) なら \(d=0\)。\(H_3\ne0\) に対し
\((U_1+1)H_3=0\)、\(\operatorname{tr}U_1=0\)。
重複固有値なら両方零で、三次多項式への作用も nilpotent なので eigenvalue −1 を持たない。
相異なる固有値 \(\pm\mu\) では三次重みは \(\pm\mu,\pm3\mu\)。
従って符号・固有方向交換を除いて \(\mu=1\) 又は \(1/3\) だけである。

各 \(U_1\) は trace零・determinant \(-\mu^2\) の全二次元 conjugacy orbit。
任意 translation 二次元を戻して \(U\) の全基底は既約四次元である。
\(U_1\) 可逆なので中心 \(q_0\) は一意。
中心座標の三次 \(-1\)-eigenspace は一次元で、
\(\mu=1\) は \(z^2w\)、\(\mu=1/3\) は \(z^3\)。
非零 \(H_3\) はこの固有線 bundle の零 section を除く全体である。
\(H_{\le2}\) は六係数全て任意、\(F=(U+1)H/g_0\) のみが一意に決まる。
translation が三次項へ作用する低次寄与もこの \(F\) に含まれる。

従って両 \(S_\mu\) は既約で、次元は正確に \(4+1+6+1=12\)。
係数逆が injective なので実固定X像も12次元である。
\(H_3=0\) は除かれた零 section として保持し、既 NEW25 の別全源へ戻す。
nilpotent／zero \(U_1\) をこの非零三次 source に含めていない。

## 3. g零 scalar-cubic の12次元全 source

\(U=-\operatorname{Euler}_{q-q_0}/3\)、\(q_0\) 任意、
\(H=T(q-q_0)\) を任意 binary homogeneous cubic、\(F\) を任意二次以下とする。
\((U+1)H=0\) なので全源に入る。
中心は \(U\) から復元出来るため、\((q_0,T,F)\in\mathbb A^2\times\mathbb A^4\times\mathbb A^6\)
の coefficient map は injective。source \(S_B\) と固定X像は正確に既約12次元である。
squarefree T は非空開で、重根と T零は同じ源の全閉境界として残す。
g零の任意 \(U,H\) を scalar へ共役出来るという主張ではない。

## 4. 全 stabilizer・尺度と正確な飽和次元

既全 centralizer から、\(\operatorname{Stab}_{\operatorname{Aff}_4}(X)\) は正確に十次元で
\[
 x'=x+b,\quad y'=y+bx+\ell(q)+c,\quad q'=Bq+t,\qquad\det B\ne0.
\]
\(\Delta=\ell(q)+c-b^2/2\) とすれば \(r'=r+\Delta\)。
全 tensor への実変換を独立に戻すと、old \(q=B^{-1}(q'-t)\) の pullback を省略して
\[
 g'=\det(B)g,\quad U'=B_*U,\quad d'=d,
\]
\[
 H'=\det(B)(H-g\Delta),\qquad
 F'=F-U(\Delta)-(1+d)\Delta.
\]
\(2\times2\) Jacobian の恒等式と \(P_{x'y'}\) を直接使った式であり、
Pfaffian の変換則だけから全 tensor を推測していない。
逆は同じ affine 群に属する。

constant-g では \(H\) への補正は affine で、\(H_3\) 非零と \(U_1\) の固有値集合を保つ。
g零では \(H'\) は同じ scalar 中心での全 cubic であり、\(F'\) は依然次数二以下の任意式。
従って三 source は全域で StabX 不変。

さらに \(a_\lambda:(x,y,q)\mapsto(\lambda x,\lambda^2y,q)\) の actual pushforward と
tensor scalar \(1/\lambda\) の合成は、全 \(\lambda\ne0\) で
\[
 (g,U,H,F)\mapsto(g/\lambda,U,\lambda H,\lambda^2F)
\]
を与える。全六係数の次数・\(r'=\lambda^2r\) から直接確認出来る。
三 source をそれぞれ保ち、affine 飽和像は全非零 tensor scalar 倍で安定である。

作用写像 \(\operatorname{Aff}_4\times S\to A\) の同じ tensor を与える affine 変換は、
divergence の equivariance から StabX の一 coset に限る。
逆にその全 coset は source の不変性と係数逆により実際に fiber を与える。
従って fiber は正確に十次元、\(\dim A=20+12-10=22\)。
tensor 自身の automorphism の次元をこの fiber へ二重に足し引きしていない。

\(S_\mu\) の actual divisor は
\[
 -h^2\{g_0h^2y-g_0hx^2/2+H_3(z,w)+hH_2+h^2H_1+h^3H_0\}.
\]
cubic は \(y\) について一次で、leading coefficient \(g_0h^2\) と定数部の gcd は一。
\(H_3\ne0\) が定数部の h-divisibility を防ぐためである。
Gauss により絶対既約で、h は唯一 double hyperplane となる。
\(S_B\) の squarefree 開では divisor は \(h^2\) と相異三線の積で、
三線は \(q-q_0\) の二次元空間にあり h と異なる。従ってここでも唯一 double H。

marked H を動かす四次元を戻し、scalar 一を除き、一般点で H を忘れる fiber が一点なので
三つの PGL 飽和射影像閉包 \(Z_1,Z_{1/3},Z_B\) はいずれも**正確に既約25次元**。
固定 point の tangent 又は orbit 次元は使用していない。

## 5. 相互相違

\(Z_B\) は全 actual D が五一次式の積という proper multiplication の閉条件内。
\(Z_\mu\) の全 source には絶対既約 cubic があるから、両者は一致しない。
同次元閉既約像なので一方が他方へ含まれることもない。

二つの \(\mu\) について、同じ tensor class が両 actual source の PGL 像にあれば、
唯一 double H を保存するのでその同値は affine 変換と非零 scalar である。
modular field X を戻す全群は、上の StabX と \(a_\lambda\)／tensor 尺度で生成される。
前者は \(U_1\) を共役し、後者は \(U\) 自体を変えない。
従って固有値集合 \(\{-1,1\}\) と \(\{-1/3,1/3\}\) は同一視出来ない。
これは一般点だけでなく両非零 \(H_3\) source の全 actual 像の交わりを排除する。
もし閉包が等しければ、二つの稠密 constructible 像は同じ既約像内の非空開を含んで交わるので矛盾する。

## 6. 固定 prior129 と既受理非定数g成分の全比較

`inputs/rank_one_long_nonconstant_prior129_20261108.json` の全129行を、
その exact Git HEAD `e2e718783b6811414a7f3c4b128a04a4503e5896` の
`research/COMPONENTS.json` と直接照合した。
元 blob SHA256 は `77e77007ce29dba6e9609c3e7f46bfc720c7139f69844c67d229e6490d13007e`。
個別 id、name、rank、dimension、generic actual D projection を全て保持する。
rank四内訳は \(84\times23,18\times22,12\times24,7\times25,2\times28\)、rank二は六。
この metadata 照合と数学的因子条件の読解は別に記録する。

* 六 rank二成分は全 rank 条件の閉性により三つの rank四像を含まない。
* rank四 dim二十四以下の全百十四成分は次元により含まない。
* CA／CB の二十八次元成分では、全 actual D はそれぞれ三乗線形因子、又は二つの二乗線形因子を持つ。
  この proper な因子 incidence は全閉包に保たれる。\(Z_\mu\) は cubic irreducible でその外。
  generic \(Z_B\) も multiplicity \((2,1,1,1)\) なので、両条件の外。
* nil31-DDH は二つの二乗線形因子、half／third と positive NEW25 は三乗線形因子を必要とする閉条件内。
  同じ理由で三像の generic D は全てその外。
* 残る T01／T39／T40 は25次元。三像が入れば同次元性から一致するはずだが、
  三像は反復 H の閉条件内、各既知成分には actual reduced D の非空開があり一致しない。

既受理の非定数g25成分 \(Z_{nc}\) も別に比較する。
`RANK_ONE_LONG_NONCONSTANT_25_COMPONENT.md` と
`RANK_ONE_LONG_NONCONSTANT_AMBIENT_SEPARATE_REVIEW.md` により、その次元と componenthood は既証明。
全 rank四部分で \(D=H^2LQ\) の proper 閉 incidence に属すので、irreducible cubic を持つ両 \(Z_\mu\) を含まない。
\(Z_B\) が入れば同次元性で一致するが、\(Z_B\) は all-linear 閉条件内、
\(Z_{nc}\) の generic Q は rank五で絶対既約だから一致しない。
別の見方では、唯一 double H の開上で上の全 normalizer は g の degree 型
（零／非零定数／非定数）も保つため、各 actual source 像は互いに識別される。
因子条件の比較だけでも結論に十分である。

## 7. 受理範囲と次の有限義務

受理したのは全 source・全 stabilizer／scalar・唯一 H に基づく正確な25次元下界と、
全 prior129／既受理 \(Z_{nc}\)／三者相互の分離である。
**三像が全空間の成分であるとは、この純入口だけでは結論しない。**
既知より高次元の未発見成分の真の部分に入る可能性は、ambient 局所上界前には残る。

各全源内の構造的に指定した一つの有理点で全百五十係数／三百五十式の微分 rank百を得れば、
Euler 二十四と scalar 一を除く接上界25は既実像下界25と一致し、regular／唯一25次元成分となる。
この条件付き最弱矢印は正しいが、本監査ではその点も operator も新しく選択・実行していない。
rank百未満なら別点列へ自動継続せず、全値を保存して構造を診断すべきである。
g零の nonscalar／Jordan／零固有方向等の全残 incidence、short 源の対応する閉境界、全分類は別に保持する。
