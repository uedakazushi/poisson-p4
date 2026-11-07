# nil22 image double_C：全六磁気源の従属 HHHQ 成分への独立監査

2026-11-07。PRIMARY E_multiple_hyperplane。**PASS**。
全六磁気源とその PGL 飽和は既存の閉成分
`R4_DEPENDENT_THREE_HYPERPLANES_QUADRIC_24` に含まれる。新成分は加えない。

独立実装 `code/classification_nil22_double_c_dependent_review.py` は producer や CAS を import せず、
既存の独立 Fraction/Laurent 疎多項式 helper のみを用いる。
有限 job と結果はそれぞれ
`reports/CLASSIFICATION_NIL22_DOUBLE_C_DEPENDENT_REVIEW_JOB_20261107.json`、
`reports/CLASSIFICATION_NIL22_DOUBLE_C_DEPENDENT_REVIEW_20261107.json`。
入力の確定 SHA256 と実装全文を結果に保存した。
gate PASS 後、0.0251秒、最大 RSS 14,080KiB、全検査 PASS。

## 実族と端点

座標順は (q1,q2,x,y)。任意の六係数を残して
\(R=m_0+m_1x+m_2y+m_3x^2+m_4xy+m_5y^2\)、
\(Q=x(y+1)-1,L=x+1,K=1+\varepsilon L\) とする。
producer の族の全括弧は

\[
\begin{aligned}
\{q_1,q_2\}&=R+[-(1+\varepsilon)y-2-3\varepsilon-\varepsilon x]q_1,\\
\{q_1,x\}&=0,&\{q_1,y\}&=-Q,\\
\{q_2,x\}&=LK,&\{q_2,y\}&=-(1+\varepsilon)(y+1)^2-(1+2\varepsilon)(y+1)-\varepsilon,\\
\{x,y\}&=0.
\end{aligned}
\]

全係数は QQ[m0,…,m5,ε] の多項式である。独立実装はこの式から全四 strict Jacobi と
全350 projective Jacobi を再構成し、全て零であることを確認した。
ε=0 の全六括弧を、元の公開源表
`CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json` の case2
`nil22_image_double_C` と一致させた。m の特殊化・反転は一つもない。
符号付き Euler-Pfaffian も直接構成し

\[
D_\varepsilon=-h(x+h)\{(1+\varepsilon)h+\varepsilon x\}(xy+hx-h^2)
\]

と一致した。ただし成分帰属をこの因子の形だけからは推論しない。

## 全八源そのものへの同定

ε(ε+1)≠0 上の actual affine 座標変更

\[
X=-\varepsilon(x+1),\quad
Y=y+2+\frac{\varepsilon(x+1)}{\varepsilon+1},\quad
Z=q_2,\quad W=-\varepsilon(\varepsilon+1)q_1
\]

の行列式は −ε²(ε+1)。逆は

\[
q_1=-\frac{W}{\varepsilon(\varepsilon+1)},\quad q_2=Z,\quad
x=-\frac X\varepsilon-1,\quad y=Y-2+\frac X{\varepsilon+1}.
\]

両向きの合成、Jacobian と行列式、全六括弧の push を独立に展開して検査した。
変換後は f=X(X−1) として厳密に

\[
\{X,Z\}=f,\quad
\{Y,Z\}=Y(a_1Y+2X-1),\quad
\{Y,W\}=(a_0+a_1X)Y+f,\quad
\{Z,W\}=a_1YW-R_{\rm new},
\]

他二括弧は零であり、

\[
a_0=\varepsilon(\varepsilon+1),\quad a_1=\varepsilon+1,\qquad
R_{\rm new}=-\varepsilon(\varepsilon+1)
 R\left(-X/\varepsilon-1,Y-2+X/(\varepsilon+1)\right).
\]

これは既済 `THREE_HYPERPLANES_QUADRIC_DEPENDENT.md` §6,(12) の全八源である。
符号はその源の \(B=-a'YW+R_{\rm new}\) が \(\{Z,W\}=-B\) に現れる規約と一致する。
源の immutable snapshot
`CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_SOURCE_20261107.json` に含まれる
`projective_source()` の関数 AST を、手で再構成した全八源の式と照合した。
その SymPy コードを実行・import してはいない。
全八 parameter の literal binding と全六括弧の一致を独立に検査した。

六磁気係数の変換は行列式 ε²(ε+1)⁶ を持ち、全磁気を保持する。
対象源の既約 quadric chart も
\(a_0(a_0+a_1)=\varepsilon(\varepsilon+1)^3\ne0\) で満たす。
分母の比較には v=ε+1 を導入し、ε と v の負冪を全て払ってから v=ε+1 を代入した。
従って有限標本や数値比較ではなく、その開集合上の正確な恒等式である。

## 閉性による全支持の結論

既済の独立証明
`CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_REVIEW_20261107.json` と上記ノート §8 は、
**全八源の PGL 飽和像閉包そのもの**が24次元全成分であることを証明している。
その成分性は今回の再証明対象でなく、明示された既済依存である。
従って任意の固定 (m0,…,m5) に対して、上の実族は ε≠0,−1 で同成分に入る。
族は全 ε で正則かつ非零であり、同成分は閉なので ε=0 も入る。
m は任意なので全六源を覆う。成分が PGL 飽和閉包なので全 PGL 像も覆う。

ε=−1 は座標変更の chart の外だが実族自体は正則である。
結論には稠密開 ε≠0,−1 と閉性だけで足り、ε=−1 の別成分帰属を追加仮定しない。
源の境界や磁気の零枝を除かず、局所環、jet、平坦性、divisor の類似性を用いていない。
