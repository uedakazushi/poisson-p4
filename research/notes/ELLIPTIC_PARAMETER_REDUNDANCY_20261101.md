# 楕円二族の完全な一般表示重複：ρ₁=ρ₂=3000

日付：2026-11-01。対象は複素数体上、滑らかな楕円曲線に対応する一般パラメータ。
コード：`code/elliptic_parameter_redundancy.py`。
証明書：`reports/ELLIPTIC_PARAMETER_REDUNDANCY_20261101.json`。

## 1. 結論と数えているもの

v29 の正規化表示

\[
 \Psi_k:(\lambda,p,c,S)\longmapsto S\cdot\operatorname{Norm}_p(c\Pi_{k,\lambda}),
 \qquad p\in\mathbf P^4,\ c\ne0,\ S\in\operatorname{Sp}_4
\]

を、その有理な標準 chart 上で考える。ここで像は **全 Poisson 双ベクトル**の
S61 正規化であり、一次ジェットだけではない。両族について一般の表示重複は

\[
 \boxed{\rho_1=\rho_2=3000}.
\]

従来の 50 は λ を固定する一般安定化群の位数だった。残りの因子 60 は λ 上の
icosahedral 群の軌道長である。一般の一次ジェット上で、五変数模型の長さを N_k、
対応する楕円成分から W16 への次数を d_k とすれば

\[
 N_k=3000d_k.
\]

d₁,d₂ 自体はこの計算では決まらない。全 S61 ファイバーの有限性・尽くしにも
結論を加えない。以下では正規化群内の同値を厳密な係数計算で構成し、群外の同値を
零点集合の幾何で排除する。

## 2. 一次資料と採用する族

Ortenzi–Rubtsov–Tagne Pelap,
[arXiv:1001.4422v2](https://arxiv.org/pdf/1001.4422), 式 (16),(18),(19)
の係数を使う。式 (18) の B₃=+λ を採用する。式 (21) の −λ は既存の Jacobi
検証と両立しない。この符号問題は v27 の証明書で確認済みである。

添字を mod 5 とし、

\[
 \{x_i,x_{i+1}\}=a_1x_ix_{i+1}+a_2x_{i+2}x_{i+4}+a_3x_{i+3}^2,
\quad
 \{x_i,x_{i+2}\}=b_1x_ix_{i+2}+b_2x_{i+3}x_{i+4}+b_3x_{i+1}^2.
\]

係数ベクトルは

\[
 f_1(\lambda)=\left(\frac{-3\lambda^2+\lambda^{-3}}5,
 -2\lambda^{-1},\lambda^{-2},\frac{-\lambda^2-3\lambda^{-3}}5,2,\lambda\right),
\]
\[
 f_2(\lambda)=\left(\frac{2\lambda^2+\lambda^{-3}}5,
 \lambda,-\lambda^{-1},\frac{-\lambda^2+2\lambda^{-3}}5,-\lambda^{-2},1\right).
\]

5λ³f_k(λ) の六成分は次数 ≤5 の全多項式空間の基底をなす。従って各射影係数曲線は
有理正規五次曲線である。特に λ は q(5,1) では −a₂/(2a₃)、q(5,2) では a₂/b₃
から復元され、一般パラメータは一意である。六つの斉次双ベクトルの射影化が独立で
あることも、アフィン chart の係数行列の非零 6×6 小行列式で検証した。

## 3. 正規化群の全生成元と有理変換

ζ を原始五乗根、φ=1+ζ+ζ⁴=−ζ²−ζ³ とする。V=ℂ⁵ の基底に関して

\[
 C e_j=e_{j+1},\quad D e_j=\zeta^j e_j,\quad
 F_{ij}=\zeta^{ij},\quad Q=\operatorname{diag}(\zeta^{i^2})
\]

と置く。射影 Heisenberg 群 H=⟨C,D⟩ は位数 25。
N=N_{PGL(V)}(H) に対して

\[
 1\to H\to N\to\operatorname{SL}_2(\mathbf F_5)\to1,
 \qquad |N|=25\cdot120=3000.
\]

この正規化群の記述は Fisher,
[The Hessian of a genus one curve](https://arxiv.org/html/math/0610403v2),
Lemma 7.9、および
[Invariant theory for the elliptic normal quintic I](https://arxiv.org/html/1110.3520v1),
Theorem 3.6 と一致する。本計算ではさらに生成と位数を直接点検する。

具体的には、共役作用を (D,C) の順で記述すると

\[
 F\mapsto\begin{pmatrix}0&1\\-1&0\end{pmatrix},\qquad
 Q\mapsto\begin{pmatrix}1&2\\0&1\end{pmatrix}.
\]

これらが生成する SL₂(𝔽₅) の全120行列を列挙した。核が H であることは次の通り。
H を射影的に中心化する行列は D,C をそれぞれ五乗根倍する。その二つの倍数を H の
元で打ち消すと D,C と厳密に可換になるので、Schrödinger 表現の既約性からスカラー
である。共役は中心交換子を保つので、商は SL₂(𝔽₅) に入る。これで上界と生成が揃う。

pushforward の規約は A_*B(x)=A B(A⁻¹x)Aᵀ である。
A_*Π_{k,λ}=s_{A,k}(λ)Π_{k,m_{A,k}(λ)} の λ 変換と係数は次の通り。

| A | q(5,1): m | q(5,1): s | q(5,2): m | q(5,2): s |
|---|---|---|---|---|
| C,D | λ | 1 | λ | 1 |
| Q | ζ³λ | ζ⁴ | ζ⁴λ | ζ² |
| F | (1−φλ)/(λ+φ) | γ₁(1−φλ)³(λ+φ)²/λ³ | ((φ−1)λ+1)/(λ−φ+1) | γ₂((φ−1)λ+1)³(λ−φ+1)²/λ³ |

ここで

\[
 \gamma_1=\frac{-7\zeta^3+\zeta^2-6\zeta-3}{25},\qquad
 \gamma_2=\frac{\zeta^3+7\zeta^2+8\zeta+4}{25}.
\]

F² は射影的に反射 H₀:e_j↦e_{−j} であり、H₀ は λ を固定して Π を −Π に送る。
係数計算は六成分の抽出だけで終わらず、各生成元が六次元空間を保つことを全斉次
係数で点検した。その後、上表の族の恒等式を ℚ(ζ)(λ) で検証した。

F,Q の Möbius 作用は各族とも位数60であり、証明書に全60行列と生成語を収録した。
任意の元の有理変換は生成語で完全に指定される。順序 AB の意味を A∘B とすれば

\[
 m_{AB}=m_A\circ m_B,\qquad
 s_{AB}(\lambda)=s_B(\lambda)s_A(m_B(\lambda)).
\]

全3000個の射影変換は、SL₂(𝔽₅) の120個の生成語の lift に CᵃDᵇ（0≤a,b<5）を
掛ければ得られる。λ 作用の核は位数50で、既知の H⋊⟨H₀⟩ と一致する。
非恒等 Möbius 変換の固定点を除けば λ の安定化群もこの核に等しい。

## 4. 正規化群の外から同値が来ないこと

この節は計算証明書とは別の幾何学的証明である。Poisson 双ベクトルのスカラー倍を
許す射影同値 g が二つの滑らかな族の元を結ぶとする。g はその被約零点集合を結ぶ。

### q(5,1)

被約零点集合は埋め込まれた楕円正規五次曲線 E そのものである。
Markarian–Polishchuk,
[Compatible Feigin–Odesskii Poisson brackets](https://arxiv.org/html/2207.07770v2),
Proposition 2.3.1 の n=4 の証明では、拡張の中項の型 (1,4) が階数0、(2,3) が
階数2を与え、前者が E であることを直接確認できる。

射影自己同型の E への制限は忠実である。実際、点ごとに E を固定する非スカラー
射影変換の固定点集合は真の線形部分空間の有限和だが、E は既約かつ非退化である。
埋込みを保存する E の自己同型は L=O_E(1) を保存する。
そのうち平行移動は E[5] に正確に対応する。

H の位数5の元は E 上の位数5自己同型になる。楕円曲線の原点固定自己同型には
位数5の元がないため、すべて平行移動であり、位数比較から H は E[5] 全体である。
この E[5] は曲線間の同型で必ず共役される。従って gHg⁻¹=H、すなわち g∈N。

### q(5,2)

零点集合は曲線ではなく滑らかな楕円 scroll

\[
 S_E=\mathbf P_E(V^\vee)\hookrightarrow\mathbf P H^0(E,V)^*,
 \qquad \operatorname{rk}V=2,\quad\deg V=5,
\]

であり、V は安定である。これは Markarian–Polishchuk,
[Compatible Poisson brackets associated with elliptic curves in G(2,5)](https://arxiv.org/html/2310.18759v2),
Lemma 3.2.1(i)、その証明後の零点集合の記述、および Lemma 3.2.3 を使う。

この scroll の ruling π:S_E→E は内在的である。任意の有理曲線から E への写像は
定数なので、S_E の同型は ruling の P¹ を ruling の P¹ に送る。従って基底の同型 f
を誘導する。射影変換は O_{S_E}(1) を保存するから、π_*O_{S_E}(1)=V の同一視を
通して f^*V≅V（異なる二つの scroll の場合は対応する二つの V の同型）を得る。

基底上恒等な polarized scroll 自己同型は V の自己同型から来るが、V は安定なので
End(V)=ℂ。従ってこの核は射影的に自明である。これで polarized scroll の自己同型は
Aut(E,det V) に忠実に入る。H の作用も基底上忠実であり、位数5の各元は平行移動。
det V の次数は5なので、ここでも H の基底上の像は E[5] 全体となる。
二つの scroll を結ぶ g は E[5] を共役する。基底作用への単射性により、その持上げも
一意であるから gHg⁻¹=H、従って g∈N。

この証明は j≠0,1728 という追加仮定を必要としない。原点固定自己同型の可能な位数
が 2,4,6 であることだけを使う。一般 λ の安定化群が50であることは §3 の Möbius
作用から従う。特殊 λ の安定化群について一般値50を主張するものではない。

以上により、両族で滑らかな一般元の全射影同値は N の3000個から尽くされる。
λ の復元式と六次元空間の射影的独立性により、各 g に対応する λ′ と s_g は一意。

## 5. S61 表示への移行と五変数模型

g∈N に対応して

\[
 \lambda'=m_g(\lambda),\qquad p'=g(p),\qquad c'=s_g(\lambda)c
\]

とする。両辺を v29 の同じ有理 Darboux frame と primitive perspective gauge で
正規化すると、残る座標変換は一意の Sp₄ 変換である。定数項 J を保つので微分は
symplectic、primitive な一次項同士を結ぶ perspective 部分は零になるためである。

v29 の共変 frame R が可逆な開集合では、source primitive jet を b,b′ と書くと

\[
 S'=S\,R(b)R(b')^{-1}
\]

で残りの変換を有理に復元できる。よって N は表示空間に有理に作用し、その軌道は
同一の全正規化 Poisson 双ベクトルを与える。逆に二つの全正規化表示が等しければ、
正規化射影変換を合成して元の族の間の射影同値を得る。§4 によりそれは N の元であり、
上記の λ′,p′,c′,S′ に一致する。

有限個の N 元の全 chart 条件、滑らかさ、det R≠0、および必要な scalar invariant の
非零条件を同時に満たす稠密開集合に縮める。非恒等な有限射影変換の固定点集合は真の
閉集合なので、一般 p を取れば3000表示は相異なる。これで ρ_k=3000 を得る。
Sp₄ の中心 −I による追加の2倍はない。S は frame の式ですでに一意に決まる。

標数0なので一般表示ファイバーは被約であり、関数体次数も3000。
適切な不変開集合ではこの表示被覆は N を群とする有限 étale Galois 被覆になる。
これは「表示空間→楕円成分」の被覆についてであり、「楕円成分→W16」の Galois 群や
次数についての主張ではない。v29 の五変数模型は c と S を一意に消去しているので、
一般長さについて N_k=3000d_k がそのまま成立する。

## 6. 再現と独立点検範囲

```sh
.venv/bin/python code/elliptic_parameter_redundancy.py \
  --output /tmp/elliptic-parameter-redundancy-review.json
```

出力先は未存在ファイルを指定する。主な exact check は、五つの射影行列の逆行列、
六次元係数空間への作用、両族の有理恒等式、H の共役作用、SL₂(𝔽₅) の120要素、
各 λ 上の60要素、六成分の射影的独立性と有理正規五次曲線の係数独立性である。
幾何学的上界は §4 に明示した一次資料と証明に依存する。

elliptic_jets 担当が、S61 正規化の逆向き、primitive perspective の一意性、
frame/scalar に追加重複が生じないことを独立に点検した。成分次数 d₁,d₂ を決める
次の有限問題では、50だけで割る従来の暫定処理を3000に置き換える。
