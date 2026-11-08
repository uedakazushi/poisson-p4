# K4a の smooth24 親は既知 T37：完成 tensor の実 PGL・Euler 同一性

2026-11-08、PRIMARY E。**全150係数・全実Pfaffianの producer exact PASS。**
別作者による独立 QQ[θ]/(θ²+6) 再生は別 receipt で管理する。
新しい成分は計数しない。

## 1. 比較する二つの実 tensor

source は既受理の **完成した T37** 族
`reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T37.json`
を `p0=t0=t1=t2=t3=1` としたもの。座標順序を X=(h,x,y,z,w) とし、
その raw quadratic skew matrix を P_T(X) と書く。
七つの高次 Maurer–Cartan 補正を全て含む。
未完成の naive T37 tensor は比較にも証明にも使わない。
同完成族が既知成分 T37 に属すことは
`MPS_S61_CATALOGUE_20261103_v62.md` と同独立 review の受理済み入力である。

target は `K4A_SMOOTH_24_PARENT_ACTUAL_CURVE.md` の整数 actual family の e=1。
座標順序 Y=(c,a,b,u,v)、その matrix を P_K(Y) とする。
全十 bracket は
`CLASSIFICATION_K4A_ORDER_TWO_ZERO_LOCUS_20261108.json` の
`generic_actual_brackets` と同一である。
この点は先行 complete torus-slice の別作者証明により、唯一の24次元親の稠密24軌道を代表する。
新しい点・normal・minorを探して代表を変更したのではない。

## 2. 可逆線形写像と非零 tensor scale

θ²=−6 とし、どちらか一つの complex root を選ぶ。
X から Y への線形写像 Y=MX を次で定義する。

\[
M=\begin{pmatrix}
 2&-6&6&3&-3\\
 3&-18&18&6&-6\\
 0&-3\theta&-3\theta&2\theta&2\theta\\
 0&-3\theta&-3\theta&\theta&\theta\\
 11&-54&54&18&-18
\end{pmatrix},\qquad
\det M=2592,\qquad \lambda=-\frac{\theta}{36}\ne0.
\]

従って Φ:[X]↦[MX] は PGL5 の実際の元である。
QQ(θ) は複素数体の有限拡大表示であり、座標の選択・invertibility の未処理補集合は無い。
M の逆行列と両側逆積は同名 exact report に保存した。

target Y の一次 vector field W を

\[
\frac{W}{\theta}=
\left(u-\frac79b\right)\partial_c
 +(2u-3b)\partial_a
 +\left(3c-\frac{13}4a-\frac74v\right)\partial_b
 +\left(3c+\frac54a-\frac34v\right)\partial_u
 +(7u-5b)\partial_v
\]

と置く。次が **全 homogeneous quadratic coefficient の恒等式**である。

\[
 M P_T(X) M^t
   =\lambda P_K(MX)+(MX)\wedge W(MX).                 \tag{1}
\]

右末項は E_Y∧W の matrix である。
すなわち Y 座標では

\[
 \Phi_*P_T(Y)=M P_T(M^{-1}Y)M^t
             =\lambda P_K(Y)+E_Y\wedge W(Y).
\]

Euler 項は projective bivector では零なので、
**Φ_*[P_T]=[P_K]**。tensor の λ 倍も射影 Poisson 空間で同じ点を表す。
source と target を逆にした scalar の読み替えをしていない。

## 3. 全係数と実五次式の検証

一つの有限 E job を gate してから、元二つの public tensor だけを使った。
(1) の skew十entry、各全十五 monomial、合計 **150個の二次係数が全て零**。
さらに E∧差の全十trivector・全三次monomial、合計350係数も零となる。
M と M^{-1} の両側積は厳密な単位行列。
有限体、数値近似、一点での評価又は選んだ submatrix による一致ではない。

実五次 F_T、F_K は、各 matrix P の signed Euler Pfaffian

\[
 F_P=\sum_{i=0}^4(-1)^i X_i
  (P_{ab}P_{cd}-P_{ac}P_{bd}+P_{ad}P_{bc}),
 \quad\{a<b<c<d\}=\{0,1,2,3,4\}\setminus\{i\}
\]

から元括弧のまま構成した。
全126個の degree五係数について

\[
 \lambda^2F_K(MX)=\det(M)F_T(X),
 \quad\text{すなわち}\quad
 F_K(MX)=-559872F_T(X)                                \tag{2}
\]

を確認した。これは normal symbol、cofactor 又は divisor の因子分割だけの一致ではない。
(1) 自体から Poisson structure の実 PGL 同値を得、(2) はその独立した全 Pfaffian readback である。
全 literal identity の最終検証は **1.896581秒、RSS56,920KiB**。
実行時の inline source と full inputs の SHA を report に保存した。

## 4. 全 K4a source と親 component への含意

先行別作者 `K4A_COMPLETE_TORUS_SLICE_SEPARATE_REVIEW.md` は、
P(e) が全 e の actual Poisson family であり、e≠0 が全て一つの24次元 PGL5 orbit、
その閉包 C が e=0 の canonical K4a point の唯一の24次元 component と証明している。
(1) は P(1) を既知 T37 の点へ移す。
T37 は PGL5 で保たれる既知閉 component なので、この24軌道全体とその閉包 C を含む。
dim T37=dim C=24、両者は既約 component だから

\[
                         C=C_{T37}.
\]

従って P(0)∈C_T37。
`PLANE_CONCURRENT_M3_K4A_ORE_NONDOMINANCE.md` と nonzero-quartic 別作者 review により、
元の **全 K4a source** は P(0) の一つの actual23軌道に入る。
PGL5 不変性から、その全点は T37 に入る。
これは「既知 component の境界へ退化したから元の点も既知」という逆推論ではない。
今回は e≠0 の点を既知 component へ実同一視し、その **閉性**で e=0 へ戻す。

全 K4a parent の E認識が完了し、新しい成分の増分は **零**。
既受理の K22 認識、二つの corrected-triple weak壁、他の E未処理源を本同一性に混ぜない。
全 component 総数や全 E の完了は主張しない。

## 5. 写像の構成経路と証明の依存境界

M は arbitrary point／profile／moving basis の列から探していない。
source T37 の全零点は、h≠0で
`zw=−1/4,x=−z²/3,y=w²/3` の rational normal quartic、
h=0で `z=w=0` の secant line となる。
target の既知零点 [0:0:0:0:1] の Heisenberg first jet と、
`a=3c,u=b,v=9c` の零 secant・その二つの distinguished points を使った。
後二つの pure geometry は別担当の `K4A_SMOOTH_PARENT_ZERO_SECANT.md` に保存される。
この marked frame の五parameter全係数一致を一回だけ組み、
α=1、ρ=ω=0、k=−3/4、R=−32を得た。
R=(−2)^5 なので完成 T37 の t2=R,t3=R^{-1} は rational diagonal
diag(1,4,1/4,−2,−1/2) による実変換である。
これを合成し、共通scalar二を掛けると上の M になる。

全係数 frame 構成は94.734007秒で終了した。GB、jet次数の引き上げ、extra point は使わない。
構成時の全290非零形式行と正規化226式、inline source を discovery appendix として保存した。
ただし **認識定理そのものは(1)の一つの可逆写像だけに依存する**。
first-jet正常形や全frame分類の網羅性、zero-locus分類、parameter探索の一意性を
成功証明の追加仮定にする必要は無い。
別作者 reviewer は元二つの公開括弧と M、λ、W だけから全150/126を再構成出来る。

実行契約は `jobs/k4a_t37_intrinsic_frame_identification.json`。
report、discovery appendix、exact dependency manifest は同名 packet に置く。
旧 runner、旧 currentjob、producer report、主稿、global台帳及び Git はこの作業では編集していない。
