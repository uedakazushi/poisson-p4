# AB の proper incidence、全ファイバー、特異点集合

2026年11月1日。\(\mathbb C\) 上の被約多様体として扱う。

[零一次部分の尽くし](TORIC_ZERO_JET_EXHAUSTION_20261101.md)により、実際の T01 の
零一次部分の被約台は十次元既約閉集合 \(\overline{AB}\) である。
以下、これを \(B\) と略記する。
このノートは B の明示的な proper birational incidence と、その全幾何学的ファイバー、
被約な特異点集合を決定する。**incidence 射は正規化射ではない。**
正規化は有限射である一方、ここでは \(\mathbb P^1\) と \(LG(2,4)\) のファイバーが現れる。
B 自身の全域正規性や、T01 のスキーム論的零ファイバーの厚みは主張しない。

## 1. 固定 Lagrangian 平行移動平面での全線形表示

四次元シンプレクティック空間で \(L=\langle\partial_y,\partial_w\rangle\) を固定し、
位置を \(q=(x,z)^t\)、運動量を \((y,w)^t\) とする。
\(L\) の全平行移動で不変な斉次二次 Q は、各係数が x,z の二次式となる。
その十八係数に \([J,Q]=0\) を課すと rank 七の線形条件となり、全解は

\[
 Q=\sum_{i,j=1}^2 A_{ij}(q)\partial_{q_i}\wedge\partial_{p_j}
      +g(q)\partial_y\wedge\partial_w,\qquad A=DF,
\]

である。ここで F は二成分の斉次三次式、g は任意の二次式。
F の八係数と g の三係数が独立である。特に q-q ブロックは零で、
A の各行は三次ポテンシャルの勾配になる。

\(s=\operatorname{div}F\) とおくと、一意に

\[
 F=\frac14q s+J_2\nabla H,\qquad
 J_2=\begin{pmatrix}0&1\\-1&0\end{pmatrix}
\]

と書ける。s は二元二次式、H は二元四次式である。

\[
s=s_0x^2+s_1xz+s_2z^2,\qquad H=\sum_{i=0}^4 h_i x^{4-i}z^i.
\]

次の線形写像を定義する。

\[
L_H=\begin{pmatrix}
h_2/6&-h_1/4&h_0\\
h_3/2&-h_2/3&h_1/2\\
h_4&-h_3/4&h_2/6
\end{pmatrix},\qquad
C_s=\begin{pmatrix}
2s_1&-2s_0&0\\4s_2&0&-4s_0\\0&2s_2&-2s_1
\end{pmatrix},\qquad N=L_H+C_s/96.
\]

\((s,H)\mapsto N\) は八次元空間から \(\mathfrak{sl}_3\) への線形同型である。

## 2. 共変性と固定 L の AB 閉包

位置を \(q\mapsto Tq\)、運動量を余接変換で変えると、

\[
 F\mapsto T^{-1}F(Tq),\qquad s\mapsto s(Tq),\qquad
 H\mapsto (\det T)^{-1}H(Tq).
\]

\(\rho(T)\) を二次式の基底 \((x^2,xz,z^2)\) における代入作用とすれば、

\[
 N\mapsto (\det T)\rho(T)N\rho(T)^{-1}.
\]

証明書は八係数をすべて保ち、\(\mathfrak{gl}_2\) の四生成元について微分恒等式
\(\delta N=(\operatorname{tr}T)N+[d\rho(T),N]\) を確認する。
\(GL_2\) は連結なので、この有理表現の Lie 同変性は群の同変性を与える。
補 Lagrangian の変更 \(p\mapsto p+Cq\)、\(C=C^t\) は A を変えず、g だけを変える。
従って \(\operatorname{rank}N\le1\) という条件は座標や補空間の選択によらない。

AB の対角代表 \(F=(x^3,z^3)^t/3\) では

\[
 s=x^2+z^2,\quad H=(x^3z-xz^3)/12,\qquad
 N=\begin{pmatrix}0&-1/24&0\\0&0&0\\0&1/24&0\end{pmatrix}.
\]

これは rank 一で \(N^2=0\)。その \(GL_2\) 軌道の微分 rank は四である。
一方

\[
 Y=\{N\in\mathfrak{sl}_3:\operatorname{rank}N\le1\}
\]

は既約四次元の最小冪零錐である。実際、\(N=u\otimes v\)、\(v(u)=0\) の表示により、
非零部分は \(\mathbb C^3\times(\mathbb C^3)^*\) の既約な一つの二次条件を
一つの逆数スケールで割った四次元族である。従って上の \(GL_2\) 軌道は Y に稠密。
固定 L における AB の閉包は、正確に

\[
 Y\times\mathbb A^3_g
\]

となる。rank 一の九小行列式からなる ideal に関して、Q の全 Jacobi 係数が
零に簡約されることも証明書が確認している。

Y は正規で、特異点は原点だけである。正規性は次の環表示から直接従う。

\[
\mathbb C[Y]=
\left(\mathbb C[u_1,u_2,u_3,v_1,v_2,v_3]/\left(\sum u_iv_i\right)\right)^{\mathbb G_m},
\]

作用の重みは u に一、v に負一とする。重み零単項式は \(u_iv_j\) で生成され、
その関係は rank 一の二次小行列式と trace 零である。
括弧内の五次元超曲面の特異点は原点だけであり、超曲面の Cohen–Macaulay 性と
余次元一の正則性から正規。正規環の不変部分環も正規である。
非零 N の射影化は滑らかな旗多様体
\(\{([u],[v]):v(u)=0\}\subset\mathbb P^2\times(\mathbb P^2)^*\) なので非零部分は滑らか。
原点では全定義式が二次で、接空間八次元、局所次元四のため特異である。

## 3. 大域 proper incidence

\(K(Q)=\{v:\partial_vQ=0\}\) を定数平行移動核とする。
次の incidence を考える。

\[
\widetilde B=\{(L,Q):L\in LG(2,4),\ L\subset K(Q),\
                    \operatorname{rank}N_L(Q)\le1\}.
\]

上の共変性により条件は大域的に定義される。
\(LG\) の座標近傍で \(\widetilde B\) は
\(Y\times\mathbb A^3_g\times\mathbb A^3\) と同型なので、既約・正規・十次元。
その特異集合は \(N_L(Q)=0\) の部分である。

射影 \(f:\widetilde B\to Q_{45}\) は projective、従って proper である。
各固定 L の稠密な対角軌道と任意 g は既知 AB 族であり、
\(\mathrm{Sp}_4\) は \(LG\) に推移的に作用する。
従って像は正確に閉集合 B となる。一般点では \(K(Q)=L\) なので f は birational。
\(\dim K=2\) の開部分では、線形写像の核を取る操作が代数的な逆写像を与え、
f は同型である。

## 4. 全幾何学的ファイバー

非零斉次二次 Q では \(\dim K\le3\)。B の点では上の incidence により
\(\dim K\ge2\)。全場合は次の通り。

| 点 Q | \(f^{-1}(Q)\) |
|---|---|
| \(\dim K(Q)=2\) | 被約な一点 |
| \(\dim K(Q)=3\) | 被約な \(\mathbb P^1\) |
| \(Q=0\) | 被約な \(LG(2,4)\) |

第二の場合を詳しく示す。\(K=\ker\ell\) とすると、斉次性から
\(Q=\ell^2R\) で R は定数二重ベクトル。
compatibility は \(J^\sharp\ell\wedge R=0\) を与えるので

\[
 Q=\ell^2(J^\sharp\ell\wedge w).
\]

\(v=J^\sharp\ell\) は三次元超平面 K の symplectic radical。
K に含まれる Lagrangian 二平面は、v を含む二平面全体であり
\(\mathbb P(K/\mathbb Cv)\simeq\mathbb P^1\) となる。
これらすべてが incidence の rank 条件も満たす。
実際、任意のそのような L を標準位置に取ると、符号を w に吸収して

\[
Q=x^2\partial_y\wedge(a\partial_x+b\partial_y+c\partial_z+d\partial_w).
\]

混合ポテンシャルは \(F=(-ax^3/3,-cx^3/3)^t\)、
\(s=-ax^2\)、\(H=cx^4/12-ax^3z/12\) で、対応する N の全二次小行列式は恒等的に零。
証明書はこの恒等式も確認する。

ファイバーの被約性も従う。L を K に含める条件は Grassmannian の通常の
\(Gr(2,K)\) を定め、\(\dim K=3\) では isotropic 条件がその \(\mathbb P^2\) 内の
被約な直線を切る。N の条件はこの直線上で恒等的に零なので追加の厚みを与えない。
\(\dim K=2\) の場合は被約な一点、Q=0 では N 条件も消えて LG 自身である。

## 5. 被約な特異点集合

次の二つの閉集合を定める。

\[
\begin{aligned}
\Sigma_{\mathrm{plane}}
 &=\{Q=g\,u\wedge v:\langle u,v\rangle=L\text{ が Lagrangian},\
                           g\in\operatorname{Sym}^2(V/L)^*\},\\
\Sigma_{\mathrm{one}}
 &=\{Q=\ell^2(J^\sharp\ell\wedge w)\}.
\end{aligned}
\]

零も含める。第一は LG 上の rank 三ベクトル束
\(\operatorname{Sym}^2(V/L)^*\otimes\wedge^2L\) の proper な像、第二は
\(\mathbb P(V^*)\) 上の rank 三ベクトル束の proper な像であり、
それぞれ既約六次元である。一般点でパラメータ平面または \([\ell]\) が一意なので
この次元は像の正確な次元である。
第一は \(\widetilde B\) の特異集合の像、第二は \(\dim K\ge3\) の集合。

\(\dim K=2\) の開部分で f は同型なので、そこでの特異点は正確に
\(N=0\)、すなわち \(\Sigma_{\mathrm{plane}}\) の点である。
その一般点の局所模型は \(Y\times\mathbb A^6\)。
従って \(\Sigma_{\mathrm{plane}}\) 全体が特異集合に入る。

\(\Sigma_{\mathrm{one}}\) のトレース非零開部分も特異である。
これを直接示す。\(\Sigma_{\mathrm{plane}}\) はトレース零なので、この開部分を除いた
B の領域では、f の全逆像は滑らかである。
正次元ファイバーを持つ集合は \(\Sigma_{\mathrm{one}}\) の逆像だけで、
その次元は \(6+1=7\)、滑らかな十次元 \(\widetilde B\) の中で余次元三。

もしその点で B が滑らかなら、小さい滑らかな Zariski 開近傍に制限した
projective birational 射は、滑らかな十次元多様体間の射となる。
\(\det(df)\) は直線束の generically nonzero な切断であり、
\(\mathbb P^1\) ファイバーの接方向を殺すため、そのファイバーで零になる。
従ってその零集合には余次元一の部分が必要である。
しかし f は \(\dim K=2\) の部分で同型で、零集合は上記余次元三の例外集合にしか
存在できない。矛盾である。
この議論は B の正規性を仮定せず、仮に滑らかだとした近傍だけを使う。

\(\Sigma_{\mathrm{one}}\) のトレース非零部分は稠密なので、その閉包全体も特異。
残る点は \(\dim K=2,N\ne0\) であり、f の同型と Y の滑らかさから滑らか。
従って正確に

\[
 \boxed{\operatorname{Sing}(B)_{\mathrm{red}}
       =\Sigma_{\mathrm{plane}}\cup\Sigma_{\mathrm{one}}\ }.
\]

二つは相異なる六次元既約成分である。その交叉は
\(Q=\ell^2(J^\sharp\ell\wedge w)\)、\(\ell(w)=0\) の部分、すなわち
\(\mathbb P(V^*)\) 上の rank 二族の像で、既約五次元。
第一にだけ属する一般点は g が二変数に依存する場合、第二にだけ属する一般点は
トレースが非零の場合である。

## 再現と未決事項

```
.venv/bin/python code/ab_incidence_20261101.py --output /tmp/ab-incidence-review.json
```

[証明書](../../reports/AB_INCIDENCE_20261101.json)は固定 L の全十八係数からの還元、
八次元線形同型、普遍 GL2 共変性、四次元稠密軌道、rank 一 ideal による Jacobi の簡約、
\(\dim K=3\) の全係数での rank 一恒等式を保存する。
properness、全ファイバーと特異集合の議論は上記の幾何学的証明である。
B の正規化、\(\Sigma_{\mathrm{one}}\) に沿う正確な横断局所環、
実際の T01 の零ファイバーのスキーム構造は、次に切り分けて調べる課題として残る。
