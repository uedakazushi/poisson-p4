# 階数三接配置の T(w)=0 枝：六余剰クラスへの Pfaffian 制約

2026-11-07。主課題 E_quadric3 の支持に必要な補題。
`RANK_THREE_QUADRIC_SUPPORT.md` の全逆形式表示を入力とする。
以下は多項式恒等式の手計算であり、新しい計算実行はない。
この補題自体は全源の構成までを扱う。得られた全源のPGL飽和閉包が新しい24次元成分であり、
この枝の全支持を覆うことは `RANK3_TWZERO_COMPONENT.md` の独立検証付き証明に接続する。

## 命題

`Q=w+u²`、体は複素数とする。全逆形式を

\[
 \alpha=P_Q(\beta+u\,dr-r\,du),\qquad
 r\in\operatorname{Sym}^2\langle v,z,w\rangle
\]

と書く。`Pf(α)=kQ²`、`k≠0` を仮定する。
`Y=(v,z,w)`、`vol3=dv∧dz∧dw` として

\[
\begin{gathered}
 \beta=\theta\,du+\eta,\quad
 \theta=t_0+t u+b(Y),\quad \eta=u c+\Gamma(Y),\\
 s=c-b,\quad d=b+c,\quad
 d\beta=du\wedge s+\tau,\quad \tau=d_Y\Gamma=\iota_T\mathrm{vol}_3.
\end{gathered}
\]

b,c は定数一形式、Γ は affine 一形式である。T(w)=0 の枝では、
Q を保つ Y の線形座標変換と非零定数の正規化により

\[
 \tau=dw\wedge dv,\qquad T=\partial_z
\]

としてよい。この規約で、次が必要である。

1. `k=−d_z≠0`。
2. `r=w(Av+Bz+Cw)`。特に `r_vv=r_vz=r_zz=0`、二次形式 r の階数は高々二。
3. `Γ_v=mw+n`、`Γ_z=ℓw+o` は w のみに依り、
   `Γ_w=(m−1)v+ℓz+gw+h`。
4. `p=s_v/2,q=s_z/2` と置けば、ある affine 一変数式 W に対し

\[
 (Aw+p)^2+w=(Bw+q)W(w).                         \tag{1}
\]

と書ける。`Bw+q` は零多項式ではない。
B=0 の場合は A=0、q≠0。B≠0 の場合は

\[
 (pB-Aq)^2=qB.                                 \tag{2}
\]

さらにこの枝の双対 Poisson 二ベクトルは `{u,w}=0` を満たし、
base 座標 (u,w) と fiber 座標 (v,z) の間の四つの混合括弧は (u,w) のみに依る。
残る `{v,z}` は (v,z) に関して affine である。

## 証明

余剰六係数を一つも落とさずに展開すると

\[
\begin{aligned}
 \alpha&=du\wedge A_1+B_2,\\
 A_1&=Q s-2u\eta+(\theta-r)dw+2w\,dr,\\
 B_2&=Q\tau-dw\wedge\eta-u\,dw\wedge dr.
\end{aligned}                                                   \tag{3}
\]

従って `Pf(α)` は Q で割り切れる。次の Pfaffian 商 F は次数高々二である。
`K(λ,μ)=λ_v μ_z−λ_z μ_v` と置けば

\[
\begin{aligned}
 F={\operatorname{Pf}(\alpha)\over Q}
 ={}&Q s_z-2u\eta_z+2w r_z\\
    &-K(\eta,s)-2K(\eta,dr)-uK(dr,s).
\end{aligned}                                                   \tag{4}
\]

一般の T の式では `(θ−r)T(w)` も現れるが、この枝では零である。
F は r の六係数に関して一次であり、r² の項はない。
u² 係数の比較から `k=s_z−2c_z=−d_z`。

`R=r+s(Y)/2` と置く。`F=kQ` の u 係数と u⁰ 係数はそれぞれ

\[
\begin{aligned}
 2\Gamma_z+d_vR_z-d_zR_v&=0,                    \tag{5}\\
 \Gamma_zR_v+(w-\Gamma_v)R_z+{d_z\over2}w&=0.   \tag{6}
\end{aligned}
\]

ここでは `K(d,s)=2K(c,s)` を使って定数項を相殺した。
(5)を(6)へ代入すると

\[
 d_z(R_v^2+w)=R_z(d_vR_v+2\Gamma_v-2w).          \tag{7}
\]

R_v,R_z,Γ_v は affine 一次式であり、`d_z≠0` である。
R_v が v または z に非自明に依れば、R_v と w は独立な affine 座標にできるので
`R_v²+w` は `C[v,z,w]` で既約な二次多項式である。
一方(7)の右辺は二つの affine 一次式の積である。
一方が非零定数の場合も積の次数は高々一だから等式は不可能である。
従って `R_v∈C[w]` となり、`r_vv=r_vz=0`。

`R_v²+w` は零多項式ではない。右辺の両因子の (v,z) に関する次数は、
整域 `C[w]` 上で積に関して加法的である。従って両因子は `C[w]` に属し、
特に `R_z∈C[w]` となる。よって `r_zz=0`、すなわち主張2が得られる。
同じ理由で Γ_v は w のみに依り、(5)から Γ_z もそうなる。
`d_YΓ=dw∧dv` の三係数を積分すると主張3の Γ_w が得られる。

`R_v=Aw+p,R_z=Bw+q` を(7)へ代入し、d_z で割れば(1)となる。
B=0 なら q≠0 であり、右辺の次数は高々一なので A=0。
B≠0 なら w=−q/B を代入することで(2)を得る。

最後に(3)では `α_vz=0` である。四次元の双対係数を取ると `{u,w}=0`。
主張2–3を(3)へ代入すると、`α_uv,α_uz,α_vw,α_zw` は全て (u,w) のみに依り、
`α_uw` は (v,z) に関して affine である。これが括弧についての主張を与える。□

## 全支持を覆う一つの既約パラメータ源

(1)の affine 一次式を `U=Aw+p,V=Bw+q,W=Ew+f` と書く。
`VW−U²=w` と同値な全係数条件は

\[
 A^2=BE,\qquad p^2=qf,\qquad qE+Bf-2Ap=1.      \tag{8}
\]

対称行列 pencil

\[
 M(w)=\begin{pmatrix}q&p\\p&f\end{pmatrix}
       +w\begin{pmatrix}B&A\\A&E\end{pmatrix}
\]

の二つの係数行列は階数一であり、どちらも零行列ではない。
零なら(8)の最後の式に反するためである。
複素数体上の非零対称階数一行列は `xxᵗ` と書ける。
二つの列 x0,x1 を並べた行列の行列式の二乗は一なので、一列の符号を選び直せば
次の **SL2 全源** が得られる。

\[
 \begin{pmatrix}a&c\\b&d_0\end{pmatrix}\in SL_2,
 \quad
 V=a^2+c^2w,\quad U=ab+cd_0w,\quad W=b^2+d_0^2w. \tag{9}
\]

逆に(9)は恒等的に `VW−U²=w` を満たす。
従ってこれは非零開だけのパラメータ表示ではなく、(8)の全複素点を覆う。
B≠0 は c≠0、B=0 は c=0 である。後者では `ad0=1` なので q≠0、A=0となる。
二つは同じ既約 SL2 源の開部分と閉境界であり、別個の支持成分と決めつけない。

全 β と r は、(9)と自由パラメータ

\[
 \delta\ne0,\quad D,D_0,S,C,g,h,t_0,t\in\mathbb C
\]

から次のように復元できる。

\[
\begin{gathered}
 s=2ab\,dv+2a^2\,dz+S\,dw,\qquad
 d=D\,dv+\delta\,dz+D_0\,dw,\\
 b_Y=(d-s)/2,\qquad c_\eta=(d+s)/2,\qquad
 r=w(cd_0v+c^2z+Cw),\\
 \Gamma_v=(\delta W-DU+2w)/2,\qquad
 \Gamma_z=(\delta U-DV)/2,\\
 \Gamma_w=(\partial_w\Gamma_v-1)v+(\partial_w\Gamma_z)z+gw+h,\\
 \theta=t_0+t u+b_Y(Y),\qquad
 \eta=u c_\eta+\Gamma,\qquad \beta=\theta\,du+\eta.
\end{gathered}                                                   \tag{10}
\]

(10)は `d_YΓ=dw∧dv`、(5)–(6)、従って全 `Pf(α)=−δQ²` を満たす。
逆に前節の必要条件から任意の元を(9)–(10)で復元できる。
`β↦β+κdQ` は α を変えず、h を h+κ、t を t+2κ に変える。
従って h=0 としてよく、固定した `τ=dw∧dv` の全逆形式支持は

\[
 SL_2\times\mathbb G_m\times\mathbb A^7
\]

という既約な11次元源の正則像で尽くされる。
ここで像の次元が11に等しいことや、全周囲の成分性はまだ主張しない。

全自由パラメータを特殊化した簡潔な実点として、次の二つを使える。

| 点 | r | β | Pf(α) |
|---|---|---|---|
| B=0 の点 | `Cw²` | `2u dz+2w dv+v dw` | `−2Q²` |
| B≠0 の点 | `w(v+z+w)` | `2u dz+2w dv+w dz+(v+z)dw` | `−2Q²` |

前者の C は任意であり、C=1 は余剰クラスを実際に持ち、C=0 は primitive 部分である。
両者は(9)の行列がそれぞれ恒等行列と `[[1,1],[0,1]]`、δ=2、D=0、
他の自由パラメータが零（Cのみ指定値）の場合に得られる。

primitive Lagrangian 族の λ=3 の点もこの全源に含まれる。
具体的に `SL2=I,δ=4,D=D0=S=C=g=h=t0=t=0` とすると

\[
 \beta=z\,du+3u\,dz+3w\,dv+2v\,dw,\quad r=0,
 \quad d\beta=2du\wedge dz+dw\wedge dv,
 \quad\operatorname{Pf}(\alpha)=-4Q^2.
\]

`ι_Z(dβ)=β` を満たす場は
`Z=(3u/2)∂u−2v∂v−(z/2)∂z+3w∂w` であり、`Z(Q)=3Q` である。
この点に対する全周囲の成分・滑らかさ検査と、本節の全源の既約性を接続できる。

もしこの源の像が含む一点において、全 Poisson 空間が既知一成分だけを持つこと
（例えばその既知成分の全周囲での滑らかさ）が別途証明できれば、
既約像の閉包はその既知成分に含まれ、全 T(w)=0 枝の支持被覆が従う。
既知成分への一点の所属だけでは、この一意性の代用にはならない。

## 制限と次の使用範囲

この結果は全六余剰クラスを保ってから三係数の消失を証明する。
残る `wv,wz,w²` の三クラスを零とする根拠にはならない。
また r の階数二以下という条件だけで成分を同定することもできない。
(5)–(7)はこの枝の全 Pfaffian 恒等式を保持しているので、
今後の族の構成・閉包への変形にはその全式を使うことができる。
二つの base 座標が Poisson 可換である事実は、合理写像
`[h:u:v:z:w]↦[h:u:w]` の一般二次元線形 fiber が Lagrangian であることを表す。
この構造だけから既知成分閉包への所属は主張しない。
