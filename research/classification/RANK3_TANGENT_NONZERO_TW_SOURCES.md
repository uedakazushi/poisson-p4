# 階数三接配置の T(w)≠0 枝：二つの全パラメータ源

2026-11-07。主課題 E_quadric3。
Q=w+u² の逆形式の全25次元空間に Pfaffian 条件を課す。
全六余剰クラスを残した表示から、T(w)≠0、r≠0 の全点を二つの既約源へ帰着する。
混合型の全源は以下で新しい22次元成分へ、平方型の全源は独立した実変形の証明を通じて
MPS成分へ接続する。源の構成だけから成分への帰属を推測していない。

## 1. 全 Pfaffian 条件の低次数表示

Y=(v,z,w)、vol3=dv∧dz∧dw とし、

\[
 \alpha=P_Q(\gamma),\quad \gamma=\beta+u\,dr-r\,du,
 \quad \beta=a\,du+\eta,\quad r\in\operatorname{Sym}^2(Y)
\]

と書く。ここで a と η の係数は affine-linear。
\(a=a_u u+a_Y(Y)+a_0\)、\(\eta=u b+\eta_0(Y)\) と置き、
\(\sigma=b-a_Y\)、\(\tau=d_Y\eta_0=\iota_T\mathrm{vol}_3\) とする。
σ、τ、T は定数である。

全六係数の r を残して外積を取ると、Pf(α) は Q で割り切れ、商 F は二次式となる。
\(\lambda\wedge dw\wedge\mu\) は vol3 に対する係数を表すこととして

\[
\begin{aligned}
 F={}&Q\sigma(T)+2wT(r)
      -\sigma\wedge dw\wedge\eta-2dr\wedge dw\wedge\eta
      -u\sigma\wedge dw\wedge dr\\
    &+(a-r)T(w)-2u\eta(T).
\end{aligned} \tag{1}
\]

従って F は r に関して一次である。F=kQ の u² 係数は

\[
 k=-(b+a_Y)(T)\ne0. \tag{2}
\]

特に τ=0 は不可能。T(w)=0 の全枝は
`RANK3_TANGENT_ZERO_TW_RESTRICTION.md` で別に保持する。
以下 T(w)≠0 とし、Y の線形変換と非零尺度で
τ=dv∧dz、T=∂w と正規化する。この変換は Q を保つ。

\[
 X=\eta_{0,z}\partial_v-\eta_{0,v}\partial_z,\qquad
 m=b+a_Y,\qquad \sigma_h=\sigma_vv+\sigma_zz
\]

と書けば div(X)=1。F=kQ の u⁰ と u¹ の条件はそれぞれ

\[
\begin{aligned}
 (2X+2w\partial_w-1)r+w\sigma_w+X(\sigma_h)+a_Y(Y)+a_0&=kw, \tag{3}\\
 2\eta_{0,w}&=a_u-\sigma\wedge dw\wedge b
 +(m_z\partial_v-m_v\partial_z)r. \tag{4}
\end{aligned}
\]

(2)–(4) がこの枝の全 Pfaffian 条件であり、いずれかの係数を省いていない。

## 2. 非零 r は平方型または w と一次式の積

Z0=X+w∂w の線形部分を L0 とする。
L0 の二次元 (v,z) block の trace は1、w の固有値は1である。
(3) の二次係数から \(L_0(r)=r/2\)。
二次同次多項式上の固有値は、L0 の三固有値の二つの和である。
従って r≠0 なら (v,z) block の固有値は

\[
 (1/4,3/4)\quad\hbox{または}\quad(-1/2,3/2).
\]

どちらも w の固有値1を含めて相異なるので L0 は半単純である。
重み1/2の二次固有空間は一次元で、前者では r=L²、後者では r=wL となる。
この時点では L は Y の一次式であり、w 項を持つ可能性をまだ残す。

r=L² の場合、(4) は η0,w=eL+定数を与える。
dYη0=dv∧dz より、X の w 係数は (eL_z,−eL_v) である。
従って X(L) の w 係数は零。L0(L)=L/4 の w 係数を比較すると L_w=0となる。
r=wL の場合、(4) は η0,w=ew+定数を与え、X は w に依らない。
L0(L)=−L/2 の w 係数を比較して再び L_w=0。

τ=dv∧dz を保つ二変数線形変換で、従って非零 r はちょうど

\[
 r=v^2 \quad\hbox{または}\quad r=wv \tag{5}
\]

に正規化できる。r=0 は別の primitive 枝であり、非零 r の議論から削除しない。

## 3. 平方型 r=v² の全源

自由な六パラメータを (b_v,b_z,b_w,A,B,C) とする。

\[
\begin{gathered}
 a_v=(16Ab_z-16C-b_v)/3,\quad a_z=-3b_z,
 \quad a_w=-2b_w-8b_z^2,\\
 \sigma_v=b_v-a_v,\quad\sigma_z=4b_z,\quad
 a_0=-\sigma_vC+4b_zB,\quad a_u=\sigma_zb_v-\sigma_vb_z,\\
 \eta_0=(Av-3z/4-2b_zw+B)\,dv+(v/4+C)\,dz-2b_zv\,dw,\\
 a=a_uu+a_vv+a_zz+a_ww+a_0,\qquad
 \eta=u(b_vdv+b_zdz+b_wdw)+\eta_0.
\end{gathered} \tag{6}
\]

すると \(\operatorname{Pf}\alpha=(b_w+8b_z^2)Q^2\)。
非退化開は k=bw+8bz²≠0 であり、これは既約六次元源である。
逆に (2)–(4) の全係数を比較すると、任意の解は (6) に、
η0,w の定数 D と a_u の加数2Dだけを加えた形になる。
これは β に DdQ を加えただけなので α は変わらない。従って (6) は全解を覆う。

## 4. 混合型 r=wv の全源

同じ六パラメータに対し、今度は

\[
\begin{gathered}
 a_v=(b_v-4Ab_z)/3,\quad a_z=3b_z,\quad a_w=-2b_w-2C,\\
 \sigma_v=b_v-a_v,\quad\sigma_z=-2b_z,\quad
 a_0=-\sigma_vC-2b_zB,\quad a_u=\sigma_zb_v-\sigma_vb_z,\\
 \eta_0=(Av-3z/2+B)\,dv+(-v/2+C)\,dz+2b_zw\,dw,\\
 a=a_uu+a_vv+a_zz+a_ww+a_0,\qquad
 \eta=u(b_vdv+b_zdz+b_wdw)+\eta_0.
\end{gathered} \tag{7}
\]

と置く。\(\operatorname{Pf}\alpha=(b_w+2C)Q^2\) であり、k=bw+2C≠0 が既約源を与える。
やはり (2)–(4) の全係数比較によって全解が得られ、DdQ の核方向だけを除いている。

## 5. 全パラメータ検証と支持への接続

`code/classification_double_quadric.py rank3_extra` は (6),(7) の全パラメータを保持し、
α の次数≤2、Pfaffian、対応する二次斉次 lift の全 Jacobi 恒等式を QQ 上で確認した。
二つの源にそれぞれ全体尺度4,2を掛け、bw=1、残りのパラメータ零の点で
全350×150接行列と真の六方向の像を検査した。

| r | 全 Jacobi 階数下界 | 軌道・鉛直・尺度の階数 | 六源方向を含む階数 | 射影接次元上界 |
|---|---:|---:|---:|---:|
| v² | 101 | 47 | 47 | 24 |
| wv | 102 | 47 | 47 | 23 |

全 D F=0 は厳密な恒等式、有限体 minor は標数零の階数下界だけに使った。
証明書は `reports/CLASSIFICATION_RANK3_TANGENT_EXTRA_SOURCES_20261107.json`、実測約1.24秒。
契約は `jobs/rank_three_tangent_sources.json`、計算前 gate は E_quadric3 active 下で PASS。

この表だけから成分性を結論しない。以下の追加証拠が全支持への橋を与える。

## 6. 混合型の全源は新しい22次元成分に入る

(7) の bw=1、他の五パラメータ零で得る代表は、上表と同じ全体尺度2により

\[
 (P_{uv},P_{uz},P_{uw},P_{vz},P_{vw},P_{zw})
 =(-hv,3hz-2uw,2hw+2u^2,2hw+2u^2+2wv,-2uv,4w^2+6uz).
\]

`RANK3_MIXED_COMPONENT_REVIEW.md` は全350方程式・全群作用を別実装で再構成し、
全微分 D の QQ 階数102、実軌道・鉛直・尺度の階数47、これに一法線 R を加えた階数48を確認した。
さらに全方程式上の十項の整数線形汎関数 λ が
\(\lambda D=0\)、\(\lambda F(R)=-140\ne0\) を満たす。
従って、実軌道方向を切った真の一変数形式切片は非零二次式を持ち、零次元である。
全局所次元は射影化後22以下、実 PGL 軌道は22次元なので、その閉包 Cmix は全空間の22次元成分である。
同じノートが従前58成分との相違も証明している。

この点を通る任意の全空間成分 C は、PGL5 の連結性によって当該22次元軌道の閉包を含む。
一方、この点での全局所次元上界22から dim C≤22。従って C=Cmix となり、
この点を通る被約成分は一つだけである。非被約接空間を無視した滑らかさは主張しない。

(7) の非退化六パラメータ源の像閉包は既約でこの点を含む。
その像閉包を含む全空間成分もこの点を通るため、Cmix に一致する。
従って **(7) の全点**が Cmix に属する。
第2節の全点正規化を戻せば、T(w)≠0、r=wL≠0 の全支持は Cmix で被覆された。
一般パラメータ点での微分階数だけから全像を決めた議論ではない。

## 7. 平方型の全源への接続

`RANK_THREE_QUADRIC_SUPPORT.md` 第5節は、(6) の指定点における実二次 Poisson 曲線と
SNC 軌道極限、全周囲の二法線切片の初期式を用いて、その点を通る唯一の23次元成分を同定する。
同じ既約源の論理により、(6) の全像はこの MPS 成分に入る。
有限の接空間方向や図式の一致だけで所属を主張していない。

第2–4節の正規化と全係数比較は、Tw0担当者と平方型担当者が独立に読解し、
余剰 L_w の排除および DdQ 以外の原始形式自由度を落としていないことを確認した。
