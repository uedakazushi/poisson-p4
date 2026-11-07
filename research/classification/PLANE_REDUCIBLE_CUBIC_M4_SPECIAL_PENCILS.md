# 被約可約 cubic plane・m4 squarefree special j：全 leading pencil

2026-11-07。主義務 G、producer `plane_collision_review`。
\(S=\mathbb P^2\) 上の restriction rank は二、zero cubic \(R=(f=0)\) は被約可約とする。
plane 沿いの実退化因子の最初の normal order は四、
\(G\in H^0(O_S(1))\otimes\operatorname{Sym}^4\mathbb C^2\) は generic squarefree とする。
本ノートは **\(I(G)\equiv0\) 又は \(J(G)\equiv0\) の全 special j 枝** を有限化する幾何証明である。
concurrent、m2/m3、重根、非被約 R は別残余として保持する。
full Poisson source の次元上界は別の有限行列証明書で検査する。

## 1. 既受理の全接続入口と invariant

`PLANE_REDUCIBLE_CUBIC_CONNECTION_ENTRANCE.md` の全 Hamiltonian associated-graded 論証により
\(\gamma=G/f\in\operatorname{Sym}^4N^*\otimes\det N\) は conormal flat connection に parallel。
これは全 simple-pole 接続について成り立ち、log 性は仮定しない。
generic squarefree なので \(\Delta(G)\) の零 divisor は R の成分に支えられる。
support が R 全体に等しいとは主張せず、unused components を全て許す。

normal 二変数を \(u,v\) とし、binomial convention を

\[
 g=a u^4+4b u^3v+6c u^2v^2+4d uv^3+e v^4
\]

とする。古典 invariant は

\[
 I=ae-4bd+3c^2,\qquad
 J=ace+2bcd-ad^2-eb^2-c^3,
 \quad\Delta\sim I^3-27J^2.
\]

係数五つは base 一次式で、全線形写像の像を
\(W\subset\operatorname{Sym}^4\mathbb C^2\) とする。
generic squarefree は \(\Delta|_W\not\equiv0\) と同値。
従って \(I,J\) が両方恒等零の枝はここには無い。

## 2. \(I\equiv0\)：rank 二以下、全 stable pencil の二型

I は五変数上 nondegenerate quadratic form である。
W は totally isotropic だから Witt index より \(\dim W\le2\)。
\(\dim W=1\) は後の scalar 枝へ送る。
二次元なら \(J|_W\) は零でない binary cubic、従って少なくとも一つの root を持つ。
その root の quartic は \(I=J=0\)、すなわち triple root を持つ。
この最後の古典的 nullcone 判定も直接確認出来る。
squarefree quartic の三 root を \(0,1,\infty\) に置いた Legendre 形では I と J は同時に零でなく、
double+double 又は double+simple+simple 型は直接代入で I,J 同時零にならず、
残る triple+simple 又は fourth-power の二型だけが同時零となる。

root quartic を normal GL2 で次の二つに正規化する。

**Fourth power \(A=u^4\) を含む場合。**
もう一つの元 \(B=(a,b,c,d,e)\) は isotropy より
\(e=0\)、\(-4bd+3c^2=0\)。
\(d=0\) なら \(c=0\) で pencil 全体が \(u^3\) を因子に持ち、generic squarefree に反する。
\(d\ne0\) なら \(v\mapsto v+t u\) で c を零に出来、isotropy が b も零にする。
A の倍数を引いて、全 pencil は

\[
 W=\langle u^4,uv^3\rangle,\qquad
 J(G)\sim -L_1L_2^2. \tag{1}
\]

となる。\(L_1,L_2\) は independent base linear forms。

**Triple+simple \(A=4u^3v\) を含む場合。**
isotropy は \(d=0\)、\(ae+3c^2=0\)。
B から A の倍数を引いて b=0 とすれば

\[
 J(sA+tB)=-e s^2t-4c^3t^3. \tag{2}
\]

\(e=0\) なら \(c=0\) で全 pencil が unstable。
\(e\ne0,c=0\) は (1) と変数交換で同型。
\(e,c\ne0\) なら (2) の三 root は相異なる。
代表は \(A=4u^3v\)、\(B=-3u^4+6u^2v^2+v^4\)、
\(J=-s^2t-4t^3\) である。
非零 c,e の残る尺度は diagonal GL2 と B の尺度で除けるので、これがもう一つの stable orbit。
以上で **stable I0 pencil に perfect-cube J は無く**、全 root multiplicity は
三相異又は二重+単純の二場合である。

I0 では \(\Delta\sim J^2\)。三相異の場合、その support は三つの concurrent lines である。
従って nonconcurrent triangle にも conic+line にも入らない。
(1) は二つの異なる line components を要するので conic+line には入らず、
triangle では二つの成分を使い残りの一つは unused となる。
重複 root を捨てる議論ではなく、全 multiplicity を分類している。

## 3. \(J\equiv0\)：singular symmetric spaces と二つの stable pencil

J は Hankel determinant

\[
 \det\begin{pmatrix}a&b&c\\b&c&d\\c&d&e\end{pmatrix}.
\]

まず、三次 symmetric matrix の全 singular linear space Z に関する elementary lemma を使う。
rank 二の元を congruence で \(A=\operatorname{diag}(1,1,0)\) とする。
\(B=\left(\begin{smallmatrix}p&q&r\\q&s&t\\r&t&h\end{smallmatrix}\right)\in Z\) に対し
\(\det(A+zB)\equiv0\) の一次・二次係数は \(h=0\)、\(r^2+t^2=0\)。
linear image \(B\mapsto(r,t)\) は isotropic line 以下である。
image 零なら全 matrix は共通 kernel を持つ。
image が一次元なら符号を替えて \(t=i r\) と出来、
\(\det B=r^2(p+2iq-s)=0\) を Z 上の polynomial identity として用いると
\(p+2iq-s=0\) が全 Z 上で成り立つ。
従って全 bilinear form は共通二平面
\(K=\langle(1,i,0),(0,0,1)\rangle\) 上で零となる compression 型。
rank 二の元が無い space は rank 一以下で、stable quartic を含まない。
これで common-kernel/compression の二場合を漏れなく尽くす。

Hankel 五空間との intersection は両場合とも二次元以下である。
common kernel \(q\ne0\) の三条件の transpose は quadratic q による
\(\operatorname{Sym}^2\to\operatorname{Sym}^4\) の乗法で、単射なので三条件は独立。
compression の二平面を independent binary quadratics \(q_1,q_2\) で書けば、三条件は
\(q_1^2,q_1q_2,q_2^2\) への pairing である。
これらも独立：依存があれば非定数 rational function \(q_1/q_2\) が C 上の二次方程式を満たし定数になってしまう。
従ってここでも \(\dim W\le2\)。

二次元の全 stable pencil は次の二つ。

* common kernel の quadratic が rank 二なら \(q=uv\)、\(b=c=d=0\) となり
  \(W=\langle u^4,v^4\rangle\)（secant）。\(I(G)=L_1L_2\)。
  q が rank 一なら \(W=\langle u^4,u^3v\rangle\) となり全体が unstable。
* compression 二平面 K の annihilator quadratic が rank 二なら
  \(K=\langle u^2,v^2\rangle\)、\(a=c=e=0\) となり
  \(W=\langle u^3v,uv^3\rangle\)（compression）。\(I(G)\sim-4L_1L_2\)。
  annihilator が rank 一なら \(K=\langle u^2,uv\rangle\)、\(a=b=c=0\) で全体が unstable。

両 stable 型は distinct である。secant の二 endpoint は fourth powers、
compression の二 endpoint は triple+simple であり、定数 GL2 は root multiplicity を変えない。
従って **secant だけとする分類は誤り** であり、ここでは全 compression 型も保持する。
両方とも \(\Delta\sim L_1^3L_2^3\) なので二つの異なる line components を要し、
conic+line には入らない。triangle では任意の二成分を使える。

## 4. rank 一と nonconcurrent 全幾何結論

\(\dim W=1\) なら \(G=L g_0\)、\(g_0\) は constant squarefree quartic。
\(\Delta\sim L^6\) より L は R の line component である。
parallel 性と squarefree quartic の projective Lie stabilizer 零から、projective 接続は零。
scalar 接続は \(\beta-\beta_{\rm can}=-d\log(L^3/f)/6\) に一意に固定される。
これは I,J の非零性を使わない rank 一の直接論証なので、
旧 `PLANE_REDUCIBLE_CUBIC_M4_SCALAR_REVIEW.md` の同じ固定 P0・全30/9の上界をそのまま適用出来る。
conic+line 二型の上界24,23、triangle の上界23である。
未使用の conic や第三 line があってよく、support=R を仮定していない。

rank 二の全 nonconcurrent 残余は triangle 上の次の **三つだけ** である。
triangle の finite permutation により、使う二成分を x,y、未使用を z、\(f=xyz\) と出来る。

| 型 | 全 leading space | 非零 invariant |
|---|---|---|
| J0 secant | \(G=A x u^4+B y v^4\) | \(I=ABxy\) |
| J0 compression | \(G=4A x u^3v+4B y uv^3\) | \(I=-4ABxy\) |
| I0 mixed | \(G=A x u^4+4B y uv^3\) | \(J=-AB^2xy^2\) |

A,B は任意の非零定数で、次の上界計算でも両方を自由に残す。
I0 三相異型だけを concurrent に送ったのであって、special j 全部の空性は主張しない。

## 5. 三つの一意な固定接続・正しい raw P0

normal covector の対角接続を \(\beta_u,\beta_v\) とする。
\(u^kv^{4-k}\) に対する \(\gamma\) の weight は
\((k-1)\beta_u+(3-k)\beta_v\) であり、scalar convention は
\(\beta_{\rm can}=d\log f/3\)。
parallel 性は base coefficient の log derivative が
\(d\log f-(k-1)\beta_u-(3-k)\beta_v\) に一致することである。
全三型に対し、次の tuple は \((dx/x,dy/y,dz/z)\) の係数を表す。

| 型 | \(\beta_u\) | \(\beta_v\) |
|---|---|---|
| J0 secant | (1/8,3/8,1/2) | (3/8,1/8,1/2) |
| J0 compression | (0,1/2,1/2) | (1/2,0,1/2) |
| I0 mixed | (1/6,1/3,1/2) | (1/2,0,1/2) |

全 Euler sum は一なので \(\beta_u-\beta_{\rm can}\)、\(\beta_v-\beta_{\rm can}\) は horizontal。
全五 monomial \(k=0,\ldots,4\) を代入すると、base 一次 coefficient を持てるものは
上の二 monomial だけである。
各 displayed connection は flat で全 leading tensor に parallel。
squarefree quartic の twisted GL2 Lie stabilizer も零なので、それが実接続の全体であり、
未処理の off-diagonal/nonlog 接続 parameter は無い。

raw homogeneous tensor の mixed columns は、この二つの閉 scalar 一形式から一意に復元する。
\(\tau=(r_x,r_y,r_z)\)、\(\sum r_i=0\) に対応する全 linear base field は

\[
 V_\tau=\operatorname{diag}(r_z-r_y,\ r_x-r_z,\ r_y-r_x)(x,y,z).
\]

これは \(i_Ei_{V_\tau}\mathrm{vol}=3f\tau\) による、旧 scalar P0 と同じ factor 三の convention。
\(P_f\) に mixed brackets \(\{x_i,u\}=-u(V_{\tau_u})_i\)、
\(\{x_i,v\}=-v(V_{\tau_v})_i\) を加える。
二 columns が異なるため、**raw \(\{u,v\}\) を初めから零としてはいけない**。
normal degree 二の actual Pfaffian が零となる唯一の \(k uv\) を加える必要がある。
係数 k への actual Pfaffian の係数は \(3fuv\) なので、一意に決まり追加 parameter は無い。
この P0 は全座標 logcanonical であり、全350 Jacobi と actual Pfaffian 零を literal に確認する。

## 6. 有限行列で閉じるべき残る矢印

各全接続は P0 一つに固定された。
既受理の global filtration の全 weight 一30・weight 二9の source を残し、
全350 Jacobi weight 一と D3=0 の24行、その完全核を取る。
weight 二では全350 Jacobi に加え、全15 leading coefficient のうち
上の二 monomial 以外を殺す13行を課す。
P2 の linear operator は各任意 P1 について同じで、inhomogeneous fiber は空又は同一 kernel の translate。
P1 の必要 quadratic compatibilities は original363行の左 functional を保存して構成出来る。
真正 quadratic ideal member の有限 rowspace の leading monomials を使う場合も、GB の完了を仮定しない。

固定接続を保つ connected subgroup は、base triangle diagonal torus 二、normal diagonal 二、
normal-to-base shear 六を合わせて次元10。
normal GL2 全四を引いてはいけない（固定した非 scalar 接続を一般には動かす）。
PGL5×source の fiber にはこの実10次元群が自由に入るので、
固定 P0 の全 source dimension 上界 d が得られれば、全像上界は \(24+d-10\)。
この有限検査を `code/classification_plane_reducible_cubic_m4_special.py` が実行する。
結果が必要な上界を与えない場合は、正確な全式を保存し、GB/jet/点列へ自動継続しない。
