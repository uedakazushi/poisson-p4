# 重複超平面のモジュラー分岐と全支持

2026-11-07。課題 E_multiple_hyperplane。全分類の完了は主張しない。
入口と全低次 Jacobi 条件は `REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` §3–5。
このノートはその全 incidence の補集合を保ち、選んだ例だけを被覆する議論を避ける。

## 1. 最高二次部の Koszul 表示が成り立つ範囲

重複超平面を無限遠へ取り、全構造を

\[
 P=C+L+Q,\qquad Z=\operatorname{div}P=c+Ax,
 \qquad \operatorname{tr}A=0,\quad \operatorname{Pf}Q=0
\]

と書く。最初の symplectic 基点では C=J であるが、平行移動後には C の階数を仮定しない。
全条件 `[P,P]=0`、`[Z,P]=0` と `Q∧Ax=0` を保持する。

**補題1。** `rank A≥3` なら、全最高二次部は `Q=Ax∧Bx` と書ける。
A が冪零でなければ `AB=BA`、`tr B=−1`。

**証明。** rank4 は正則列の Koszul 完全性。rank3 では定数ベクトル空間を
`im A=W` とその一次元補空間に分け、`X=Ax=∑_{i=1}^3 f_i e_i` と書く。
三つの一次式 f_i は独立で正則列をなす。`Q=Q_W+V∧e_4` と分けると
`X∧Q=0` は `X∧V=0` と `X∧Q_W=0`。
前者の Koszul 完全性から `V=gX`、後者の三項 syzygy の完全性から
`Q_W=X∧Y_W`。次数を取れば g と Y_W の係数は一次である。
従って `Q=X∧(Y_W+g e_4)`。

div Q=X と tr A=0 より

\[
 [X,Bx]=(1+\operatorname{tr}B)X. \tag{1}
\]

右辺の係数を s とする。行列関係に A^{m-1} を掛け trace を取ると
`s tr(A^m)=0`、m≥1。A が冪零でなければ、Newton 恒等式からいずれかの
`tr(A^m)` は非零であり s=0。従って可換性と trace 条件を得る。□

rank A≤2 にこの補題を適用してはいけない。rank3 冪零でも Koszul 表示は有効だが、
`[A,B]=0` を自動的に課してはいけない。

## 2. 可逆・正則半単純・対を持つ共鳴の全点被覆

非共鳴枝は入口ノート §6 で処理済み。対和零があり一次共鳴を持つ可逆正則半単純
スペクトルは、比例と置換を除き `(1,−1,2,−2)`。
A の唯一のモジュラー零点へ平行移動して

\[
 A=a\operatorname{diag}(1,-1,2,-2),\qquad a\ne0
\]

とする。補題1の可換性により B も対角なので Q は全て log-canonical。
A 不変な定数部は `c∂_{12}+d∂_{34}`、一次部は `p x_2∂_{14}+q x_1∂_{23}` に尽くされる。
従って全 P は、順序 `(12,13,14,23,24,34)` で

\[
 (c+k_{12}x_1x_2,\ k_{13}x_1x_3,\ p x_2+k_{14}x_1x_4,
 q x_1+k_{23}x_2x_3,\ k_{24}x_2x_4,\ d+k_{34}x_3x_4).
\]

モジュラー固有値の四式と `Pf Q=0` は、全ての k を

\[
 (k_{12},k_{13},k_{14},k_{23},k_{24},k_{34})
 =(b,y,-a-b-y,-2b-y,a+3b+y,-2(a+b)) \tag{2}
\]

へ制限する。分母として使ったのは a のみである。
全四つの affine Jacobi 多項式はそれぞれ

\[
\begin{split}
 &-2bc x_3-q(b+y)x_1^2,\\
 &2bc x_4+p(a+2b+y)x_2^2,\\
 &(d(a+b)+pq)x_1-2p(a+2b+y)x_2x_3,\\
 &-(d(a+b)+pq)x_2+2q(b+y)x_1x_4.
\end{split} \tag{3}
\]

従って全 Jacobi はちょうど

\[
 bc=0,\quad q(b+y)=0,\quad p(a+2b+y)=0,
 \quad pq+d(a+b)=0. \tag{4}
\]

pq≠0 なら b=−a,y=a となるが最後の式に反する。従って p または q は零。
p≠0,q=0 の場合、y=−a−2b、`bc=d(a+b)=0`。a≠0 より c,d は同時には非零でない。
p=q=0 でも c,d 同時非零は不可能。さらに c=d=0 なら P=Q は rank≤2 なので、
現在の一般 rank4 支持から除かれる。q-only は座標交換 `(1 2)(3 4)` で p-only と同じである。

次の表は残る全場合を覆う実 Poisson 族である。定数と一次部は表示された source の値を保ち、
六つの k のみを表の通り変更する。t=0 は元の任意の点であり、p,c,d の表示上の非零条件以外に
元のパラメータを特殊化していない。

| 全場合 | 変形後の六係数 | 最高 Pfaffian の係数 κ(t) |
|---|---|---|
| p≠0,c=d=0 | `(b,−a−2b,b+t,a,b,−2(a+b))` | `at` |
| p≠0,c≠0,d=0（b=0） | `(t,−a,−t,a,t,−2a)` | `−2at` |
| p≠0,d≠0,c=0（b=−a） | `(−a,a+t,−a−t,a,−a,t)` | `−at` |
| p=q=0,c≠0,d=0（b=0） | `(t,y,−a−y,−y,a+y,−2a)` | `−2at` |
| p=q=0,d≠0,c=0（b=−a） | `(−a,y,−y,2a−y,y−2a,t)` | `−at` |

各族の全 Jacobi は全パラメータで零。最高部の Pfaffian は
`κ(t)x_1x_2x_3x_4` である。t≠0 で h の重みを1、四 affine 座標の重みを0とする
PGL 極限は非退化 log-canonical bivector で、その実五次因子は
`κ(t)h x_1x_2x_3x_4`、五独立超平面の SNC である。

PGL は連結なので P_t を含む任意の全空間成分はこの SNC 極限も含む。
MPS Theorem 1.4 の SNC locus と交わる成分の尽くしを適用すると P_t は MPS40 の和に入る。
その有限和は閉じているため t=0 も属する。これでこの**全共鳴層**の rank4 支持を被覆した。

正確検証は `code/classification_multiple_hyperplane.py` と
`reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RESONANT_20261107.json`。
全係数式 (2),(3) と五族の Jacobi/Pfaffian を QQ 上の恒等式として確認した。
ジョブは `jobs/multiple_hyperplane_resonant.json`。標本点や有限体推論は使わない。

## 3. 特異・正則半単純 rank3 の全点被覆

**命題2。** A が rank3 の正則半単純行列なら、全 repeated-H rank4 支持は MPS40 の和に入る。
モジュラー場 Z が零点を持たない場合も含む。

**証明。** `A=diag(λ_1,λ_2,λ_3,0)`、λ_i は相異なる非零数、`∑λ_i=0`。
平行移動により c の im A 成分を消し、`Z=X+γ∂_4` とする。
補題1より B は対角で、

\[
 Q=X\wedge(Y+b x_4\partial_4),
 \quad X=\sum_{i=1}^3λ_i x_i\partial_i,
 \quad Y=\sum_{i=1}^3μ_i x_i\partial_i.
\]

`[Z,P]=0` の次数一次式は `[X,L]+γb X∧∂_4=0`。
前項の A-weight 零部分は零なので γb=0、次いで `[X,L]=0`。
三つの非零固有値の間に一次共鳴はない：λ_i+λ_j=λ_k は k=i,j なら一つの λ が零、
第三添字なら `2λ_k=0` を意味する。対和も非零である。
従って全一次部は `L=D∧∂_4`、`D=∑d_i x_i∂_i`。
これは x_4 と独立なので、次数零の式も `[X,C]=0` となる。
対和零がないため C=0。

全低次 Jacobi のうち `[L,Q]=0` は `bD∧X∧∂_4=0` を与える。
b≠0 なら D は X の定数倍（一次式の Koszul 完全性）となり、P 全体が X で割れて rank≤2。
rank4 の仮定より b=0。従って全 P は

\[
 P=X\wedge Y+D\wedge\partial_4,
 \qquad \operatorname{Pf}P=Kx_1x_2x_3,\quad K\ne0. \tag{5}
\]

X,Y,D は互いに可換な対角場なので

\[
 P_t=X\wedge Y+(1+t x_4)D\wedge\partial_4 \tag{6}
\]

は全 t で真の Poisson bivector。最高二次部の Pfaffian は
`tKx_1x_2x_3x_4`。§2 と同じ SNC 極限と有限閉和の議論で全点を被覆する。□

ここでは γ=0 と γ≠0 を一括して処理した。C=0 は初めから仮定したのではなく、
全不変性式と全固有値条件から従った。重複固有値がある rank3、冪零 rank3 はこの命題に含まれない。

## 4. 可逆な非半単純 A の全点被覆

**命題3。** A が可逆で非半単純の場合も、全 repeated-H rank4 支持は MPS40 の和に入る。

以下は全源を保つ実変形による証明である。可逆正則半単純の全支持は §2 と入口ノート
§6–7、可逆な非正則半単純の全支持は同 §8 を用いる。

まずモジュラー零点へ移動する。A の Jordan 分解を S+N とすると `[A,C]=[A,L]=0`
から C,L は S,N のそれぞれに不変である。補題1により A,B は可換。
`[L,Q]=0` は `Ax∧[L,Bx]=0` を与え、Koszul の次数一次部分から
`[L,Bx]=Ax∧v`（v 定数）。div L=0 を使って div を取ると Av=0、従って v=0。
よって **B も L を保存する**。

S の固有値重複型を全て分ける。

* `(a,a,a,−3a)` では S-weight から C=L=0。P=Q は rank≤2 なので除かれる。
* `(a,a,b,c)`、b,c は a および互いと異なる場合、C=0。非零 L を持つのは
  比例・置換を除き `(1,1,2,−4)` または `(2,2,−3,−1)` のみ。
  前者では L は `x_3∂_{12}` の定数倍、後者では
  `x_4(v_1∂_{13}+v_2∂_{23})`。L=0 は rank≤2。
* `(a,a,−a,−a)` では L=0。Jacobi から B は C を保存する。
  rank C=4 は tr B=−1 と矛盾、C=0 は rank≤2 なので rank C=2。

第二項の共鳴一覧は、二つの a の和が残りの固有値となる場合と、a+b=c の場合を
trace 条件へ代入するだけで得られる。非零 L の係数は座標の定数倍で1にできる。
以下の全行列は t=0 でその全正規形へ戻り、t の有限個の例外を除き A_t は可逆正則半単純。

### 重複2+1+1、L=x_3∂_{12}

重複二次元空間で `N=E_{12}` とし、

\[
 N_t=\begin{pmatrix}t&1\\0&-t\end{pmatrix},\quad
 A_t=\operatorname{diag}(aI+N_t,2a,-4a),\quad
 B_t=\operatorname{diag}(\beta I+\gamma N_t,2\beta,-1-4\beta).
\]

これは N の全可換子と L の保存条件を尽くした元の B に戻る。
A_t,B_t は可換、trace はそれぞれ0,−1。L の保存条件は第三固有値が最初の
二次元 block の trace に等しいことであり、全 t で保たれる。

### 重複2+1+1、L=x_4∂_{13}

N 不変性は N(v_1,v_2)^t=0 を課す。従って N=E_{12} と L=x_4∂_{13} を同時に選べる。
同じ N_t を使って

\[
\begin{split}
 A_t&=\operatorname{diag}\left(aI+N_t,\frac{-3a-t}{2},\frac{-a+t}{2}\right),\\
 B_t&=\operatorname{diag}\left(\beta I+\gamma N_t,
 \frac{-1-3\beta-\gamma t}{2},\frac{-1-\beta+\gamma t}{2}\right).
\end{split}
\]

最後の固有値は第一と第三の和であり、L を保存する。可換性と二つの trace 条件も恒等的。
両場合とも `P_t=L+(A_tx)∧(B_tx)` は全 t で真の Poisson bivector である。

### 重複2+2、L=0、rank C=2

S の固有空間を E_+,E_- と書く。C は両者間の rank1 pairing なので
`C=∂_1∧∂_3` とできる。N 不変性により
`N_+=nE_{12}, N_-=mE_{34}`。
B は C を保存するので各 block は上三角、`B_{11}+B_{33}=0`。
nm≠0 なら B と N の可換性から各 block の対角が等しく、tr B=0 となって矛盾する。
従って、必要なら E_+,E_- を交換し、n≠0,m=0 としてよい。

B_+=βI+γN_+、B_- の対角は `−β,−1−β`。
差は1なので、C の image 線を保つ上三角基底変換で B_- の上三角成分を消せる。

\[
\begin{split}
 N_t&=\begin{pmatrix}t&n\\0&-t\end{pmatrix},\\
 A_t&=\operatorname{diag}(aI+N_t,-a-t,-a+t),\\
 B_t&=\operatorname{diag}(\beta I+\gamma N_t,-\beta-\gamma t,-1-\beta+\gamma t).
\end{split}
\]

全 t で A_t,B_t は可換かつ C を保存し、trace は0,−1。
従って `P_t=C+(A_tx)∧(B_tx)` は全 t で Poisson。

以上の各変形は t=0 で元の rank4 点を持つため、rank4 条件も t の非空開で保たれる。
その開で可逆正則半単純の全点被覆を適用し、MPS40 の有限閉和を取れば元の点も入る。□

この節の全場合分け、全行列、trace・保存・可換条件および generic 固有値分離は、
別担当が同じ全正規形から独立に手計算し一致を確認した。追加の有限体実験や標本点を使わない。

## 5. 零固有値が単純な rank3 A の残る全枝

§3 の残りは `spec A=(a,a,−2a,0)`、a≠0。その a-eigenspace に一つの
非零 Jordan block がある場合も含める。零固有値は代数的重複度1と仮定する。
モジュラー定数を `γ∂_4` にし、三次元側を U とする。
S 不変性と補題1により

\[
 C=0,\quad L=D(x)\wedge\partial_4,\quad
 Q=X\wedge(Y+b x_4\partial_4),\quad X=A_Ux,\quad Y=B_Ux.
\]

D の係数行列も S を保存する。次数一次の不変性を A_U^{-1} 倍して trace を取ると
`γb=0`、次いで `[A_U,D]=0`。
全低次 Jacobi は

\[
 [D,B_U]-bD=sA_U \tag{7}
\]

となる（ここで角括弧は行列交換子 DB_U−B_UD）。
これは `X∧([Y,Dx]−bDx)=0` と一次式の Koszul 完全性から得られる。
rank4 条件は `X∧Y∧Dx` が非零であること。

### b=0

A_U が半単純なら、二次元重複 block の trace を式 (7) で取ると s=0。
A_U に非零 Jordan block があるなら B_U,D はその二次元可換子に入り、元から可換である。
どちらも B_U,D の二次元 block は可換。rank4 より両方とも scalar という場合はない。
従ってある非 scalar traceless M を選び、両 block を
`B_U|_2=βI+ηM, D|_2=δI+εM` と表せる。

半単純 A_U では M を semisimple へ変形 M_t し、`A_t|_2=aI+tM_t` とする。
Jordan A_U では M をその非零 nilpotent 部に選び、`A_t|_2=aI+M_t` とする。
同時に B,D の block も同じ係数で M_t を代入する。第三対角成分はそのまま、
A_t の第三成分は−2a とする。M_t の trace は常に零。

全可換性・全 trace 条件を保つため
`P_t=(A_tx)∧(B_tx)+(D_tx)∧∂_4` は真の Poisson 族。
generic t では §3 の全点被覆に入る。rank4 開条件と有限閉和で元の点も被覆する。

### b≠0

非零 Jordan block があれば B_U,D は可換なので式 (7) は D を A_U の倍数にし、
rank4 に反する。従って A_U は半単純。
二次元 block を M,N、第三成分を β,d と置く。式 (7) は

\[
 s=bd/(2a),\qquad [M,N]=-b\left(N+\frac d2I\right).
\]

`N_0=N+(d/2)I` は非零：零なら D が A_U の倍数となる。
trace と二乗 trace を取ると N_0 は非零 nilpotent なので、基底を選び N_0=E_{12}。
M は上三角で、その対角 m_1,m_2 は `m_2−m_1=b`。
b≠0 より N_0 を保つ shear で M の上三角成分を消せる。
D の A_U 比例部分は x_4 の平行移動で除ける。

trace B=−1 は `2m_1+β+2b=−1`。従って全正規形は

\[
 L=x_2\partial_{14},\qquad
 (k_{12},k_{13},k_{14},k_{23},k_{24},k_{34})
 =(ab,-a-2ab,ab,-a,ab,-2ab). \tag{8}
\]

`P_t=P+t x_1x_4∂_{14}` は全 t で Poisson：追加項は log-canonical Q と可換で、
L と同じ二つの vector directions を持つため L とも可換。
最高 Pfaffian は `−at x_1x_2x_3x_4`。§2 の SNC 議論でこの全枝も被覆する。

以上で、**零固有値が代数的に単純な rank3 A** の全点を処理した。
`(a,−a,0,0)` の零 Jordan block によって rank3 となる別の枝は、この節に含めない。

本節は統合担当が全 weight、`A_U^{-1}` trace、式 (7) の符号、b≠0 正規形と
最高 Pfaffian `−at` を独立に再導出し確認した。

## 6. rank3 の零 Jordan block も全点被覆される

残る非冪零 rank3 は、座標 `(x,y,z,w)` で

\[
 Z=a(x\partial_x-y\partial_y)+w\partial_z+\gamma\partial_w,
 \qquad a\ne0. \tag{9}
\]

平行移動で消せる三成分を消し、γ=0 と γ≠0 を保持した表示である。
補題1の可換子は、x,y の対角と、(z,w) 上の `bI+dE_{zw}`。
`B↦B-dA` の自由度も使える。
S-weight 零の全 C,L を並べ、div P=Z と `[Z,P]=0` を課すと、全 P は
順序 `(xy,xz,xw,yz,yw,zw)` で

\[
\begin{split}
 (&-a(1+2b)xy+f_0+f_1w,\quad
 x(abz+\beta w+e),\quad x(abw+\gamma\beta),\\
 &y(-abz+(1+2b-\beta)w+\eta),\quad
 y(-abw+\gamma(1+2b-\beta)),\quad
 bw^2+sw+g),\qquad s=e+\eta,                \tag{10}
\end{split}
\]

ただし `γb=γs=γf1=0`。
例えば Pxy の一次 z 係数は不変性で零、Pxw の一次 x 係数は γβ、
Pyw の一次 y 係数は γ(1+2b−β)、Pzw の一次 z 係数は−2γb=0となる。
従ってこの表示では定数や一次係数を選択的に落としていない。

`f=f0+f1w,F=bw²+sw+g` とすると、上の不変性条件を使った全四 Jacobi は

\[
\begin{split}
 &f((1+2b)w+s)+Ff_1-ab sxy,\qquad \gamma f,\\
 &-abx(sw+2g),\qquad aby(sw+2g).
\end{split}
\]

従って全条件は正確に

\[
 bs=bg=0,\quad \gamma f_0=0,\quad
 (1+3b)f_1=0,\quad (1+2b)f_0+2sf_1=0,\quad sf_0+gf_1=0. \tag{11}
\]

この条件の下で

\[
 \operatorname{Pf}P=fF+xy\{-a(sw+g)-\gamma e\}. \tag{12}
\]

### b≠0 の二枝

この場合 s=g=γ=0。一般 rank4 により f≠0 なので、全解は
`b=−1/2,f0≠0,f1=0` または `b=−1/3,f0=0,f1≠0`。
z の平行移動で e,η=−e を消す。
以下では A の (z,w) block を表示し、その他の成分は対角とする。

第一の場合は C=f0∂xy、L=0。

\[
 A_t=\operatorname{diag}\left(a,-a,
 \begin{pmatrix}t&1\\0&-t\end{pmatrix}\right),\qquad
 B=\operatorname{diag}(-\beta,\beta,-1/2,-1/2).
\]

第二の場合は C=0、L=f1w∂xy。

\[
 A_t=\operatorname{diag}\left(a+t,-a,
 \begin{pmatrix}-2t&1\\0&t\end{pmatrix}\right),\qquad
 B=\operatorname{diag}(-\beta,\beta-1/3,-1/3,-1/3).
\]

どちらも全 t で AB=BA、trace は0,−1。第一は C を保存し、第二は L を保存する。
従って `C+L+(A_tx)∧(Bx)` は実 Poisson 族で、generic t は可逆正則半単純の被覆済み枝。
有限閉和で元の全点を覆う。

### b=0 の全枝と射影的な可換場変形

(11) は f0=f1=0 を強制する。H を斉次座標とし、`D_x=x∂x,D_y=y∂y` と置く。
γ=0 の全 source は

\[
\begin{split}
 P={}&-aD_x\wedge D_y
   +D_x\wedge\{(\beta w+eH)\partial_z\}
   +D_y\wedge\{((1-\beta)w+\eta H)\partial_z\}\\
   &+H(sw+gH)\partial_z\wedge\partial_w,\qquad s=e+\eta. \tag{13}
\end{split}
\]

これは γ=0 の全点であり、s≠0 と s=0 をともに処理する。

**s≠0。** `W=w+(g/s)H` とする。

\[
 A_t=H\partial_z+t^2H\partial_H,\qquad
 B_t=W\partial_z+tW\partial_W.
\]

この二場は可換で、t≠0 では同時対角化可能。
`A_t∧B_t` は t で割り切れ、`(A_t∧B_t)/t→HW∂z∧∂W`。
`u=e−βg/s,v=η−(1−β)g/s` として

\[
 P_t=-aD_x\wedge D_y+D_x\wedge(uA_t+\beta B_t)
 +D_y\wedge(vA_t+(1-\beta)B_t)+(s/t)A_t\wedge B_t. \tag{14}
\]

係数の t=0 の見かけの極は消え、P0 は (13) に等しい。
全四場が可換なので全 t で Poisson。
その四場に関する skew coefficient matrix の Pfaffian は
`−as/t−u(1−β)+βv` であり、generic t≠0 で非零。
同時固有座標で五つの座標超平面の SNC になる。

**s=0。** (12) より g≠0。今度は

\[
\begin{split}
 A_t&=H\partial_z+t^2(H\partial_H+w\partial_w),\\
 B_t&=w\partial_z+tH\partial_w+t^2w\partial_H.
\end{split}
\]

行列として `B_t²=tA_t`、従って可換。B_t の固有値は `0,±√(t³)` なので
t≠0 で同時対角化可能であり、`(A_t∧B_t)/t→H²∂z∧∂w`。
(14) の u,v,s をそれぞれ e,−e,g に置き換えた族を使う。
Pfaffian は `−ag/t−e` で generic 非零。これも P0 が元の全 source に一致する。

**γ≠0。** この場合 s=0、η=−e であり、

\[
 B_0=w\partial_z+\gamma H\partial_w,\qquad
 A_0=H\partial_z=B_0^2/\gamma
\]

を使うと全 source は

\[
 P=-aD_x\wedge D_y+D_x\wedge(eA_0+\beta B_0)
 +D_y\wedge(-eA_0+(1-\beta)B_0)+(g/\gamma)A_0\wedge B_0. \tag{15}
\]

係数行列の Pfaffian は `−(ag+γe)/γ` で、(12)の一般 rank4 条件そのもの。
`B_t=B_0+t z∂z+2t w∂w,A_t=B_t²/γ` と置けば、t≠0 で B_t は固有値0,t,2tを持ち、
A_t,B_t は同時対角化可能。A_t,B_t の像は scalar 行列を法として独立である
（三つの固有値における 1,λ,λ² の Vandermonde 行列式）。
(15)へ代入した族は全 t で Poisson、generic t で SNC。

これらの族は全て五変数上の実二次 Poisson bivector である。
最初の二族の t=0 では分子が t で割り切れるので、射影族も正しく延長する。
各全点は SNC を持つ実既約族の像に入り、MPS の成分定理により MPS40 の和へ被覆される。
従って **非冪零 rank3 A の全支持** が被覆された。

三つの可換場変形は、別担当が全行列を独立に掛け、可換性、wedge の t 割切りと極限、
固有値、四場の独立性を確認した。§2 の certificate はこの節の計算証明書ではなく、
ここは表示した行列と独立手計算による証明である。
統合担当も S-weight で許される全 C,L の列挙から (10),(11) を独立に再導出し、
b≠0 の行列変形、b=0 の rank4 条件と三つの Cartan 極限を確認した。

## 7. 残る有限の分岐

| rank A と固有値型 | 状態・保持する条件 |
|---|---|
| 4・正則半単純・一次非共鳴 | 入口ノート §6 の全点 MPS 被覆 |
| 4・正則半単純・対和零のある共鳴 | §2 の全点 MPS 被覆 |
| 4・正則半単純・対和零なしの共鳴 | 入口ノート §7 の全点 MPS 被覆 |
| 4・固有値重複・半単純 | 入口ノート §8 の全点 MPS 被覆 |
| 4・非半単純 | §4 の全点 MPS 被覆 |
| 3・零固有値が単純 | §§3,5 の全点 MPS 被覆。Z 零点の有無を含む |
| 3・`(a,−a,0,0)`、零 Jordan block | §6 で全支持を MPS 被覆 |
| 3・冪零 | Koszul 表示は有効。式 (1) の非可換枝と c の余核成分を保持して未処理 |
| 0,1,2 | Koszul の追加類を除かず全支持を別に解析中 |

これは場合分けの表であり、未処理行を閉じたとする主張ではない。
既存 Q45 の適用にはさらに `[C,Q]=0`、`[L,L]=0` 等の実際の仮定が必要である。
rank≤2 の最高部だけから、低次項を任意に付けた全点の所属を推論していない。
