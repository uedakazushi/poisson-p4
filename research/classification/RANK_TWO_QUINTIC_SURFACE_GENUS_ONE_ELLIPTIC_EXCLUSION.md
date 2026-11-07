# 次数五・rank 二支持：正規化種数一の elliptic 全枝は空

2026-11-07。主義務 G、producer `plane_collision_review`。

**命題。** \(\pi\) を \(\mathbb P^4_\mathbb C\) 上の rank 四 Poisson tensor、
実退化因子 \(D\) を既約被約 quintic とする。
非正規支持の integral surface 成分 \(S\) は次数五で、\(\pi|_S\) の一般 rank が二とする。
正規化 \((X,L)\)、\(L=\nu^*O_D(1)\) の sectional genus が一で、
polarized birational model が elliptic \(\mathbb P^2\)-bundle 型となる全場合は存在しない。
nef だが非 ample な型、非自明な \(\operatorname{Pic}^0\) の twist、single/non-split/split conductor を含む。
同じ genus 一の canonical index 二型は本命題の明示補集合である。

これは純粋な幾何証明であり、新しい CAS 計算・有限標本は無い。
既受理入力の scope と新しい矢印を以下で区別する。
構造化記録は `reports/CLASSIFICATION_RANK_TWO_QUINTIC_GENUS_ONE_ELLIPTIC_20261107.json`、
有限監査契約は `jobs/rank_two_quintic_genus_one_elliptic.json`。
本ノートの独立監査は未実施であり、producer の証明完了として保存する。

## 1. rank 非依存の全 elliptic model

length 公式 \(\sum\delta_i\deg S_i=6-g\) より、\(g=1\)、\(\deg S=5\) のとき
\(\delta_S=1\)、他に二次元非正規支持は無い。
`RANK_ZERO_QUINTIC_SURFACE_GENUS_ONE_RECOGNITION.md` §1 と、その別作者監査
`RANK_ZERO_QUINTIC_SURFACE_GENUS_ONE_SEPARATE_REVIEW.md` §2 の全 model 還元は
restriction rank を使用しない。
Fukuma Proposition 3(iii) を smooth resolution に適用し、共通 resolution 上の polarization の一致から
全 section ring を同定する、その既受理部分だけを使う。
elliptic 型では

\[
 Y=\mathbb P_B(E),\quad p:Y\to B,\quad H=O_Y(1),\quad
 \operatorname{rank}E=3,\quad\deg E=H^3=5
\]

で、\(B\) は smooth elliptic curve、\(H\) は nef・big・globally generated。
\(E=p_*H\) も globally generated であり、section ring の ample model への birational morphism

\[
 f:Y\longrightarrow X,\qquad H=f^*L
\]

を持つ。\(X\) 自身が最初から smooth であるとも、degree 一 section ring が生成されるとも仮定しない。
以下では \(H\) の nef 境界を全部残す。

## 2. 全 zero summand と収縮 locus

非零 \(E\to O_B\) の image は globally generated な \(O_B(-Z)\) なので \(Z=0\)。
global generation は 1 へ写る global section も与えるため、この quotient は split する。
繰り返して

\[
 E=O_B^{\oplus r}\oplus F,\qquad H^0(B,F^*)=0
 \tag{1}
\]

と出来る。\(F\) は globally generated、degree 五である。
\(r=3\) は degree に反し、従って \(r=0,1,2\) の三場合が全てである。

ここで \(F\) が ample であることと、(1) の零 block 以外に収縮曲線が無いことを同時に確認する。
\(H\)-degree 零の integral curve の normalization を \(T\) とする。
fiber 上の \(O(1)\) は ample なので \(T\to B\) は非定数、従って surjective。
tautological quotient は \(g^*E\to M\)、\(\deg M=0\) で、\(M\) は globally generated、従って \(M\simeq O_T\)。
全 evaluation \(H^0(F)\otimes O_T\to g^*F\to O_T\) の係数は \(T\) 上の定数である。
もしこれが非零なら、その constant functional は \(T\) 上で evaluation kernel を消す。
\(T\to B\) の surjectivity により \(B\) の全点でも kernel を消し、非零 \(F\to O_B\) を与える。
これは \(H^0(F^*)=0\) に反する。
従って \(g^*F\to O_T\) は零で、商は \(O_B^r\) だけから来る。
その係数も定数なので、\(T\) は \(B\times\{\text{one point of }\mathbb P^{r-1}\}\) に含まれる。

同じ議論を \(\mathbb P(F)\) に行うと、globally generated な \(O_{\mathbb P(F)}(1)\) に零次数の曲線は無い。
その tautological morphism は finite、従って \(O(1)\) は ample、すなわち \(F\) は ample である。
この証明は有限 cover 後に新たな degree 零 quotient が出る可能性も評価 functional によって除いている。

以上から、\(r=0\) では \(f\) は finite birational、normal \(X\) 上同型。
\(r=1\) では唯一の zero-summand section \(B\subset Y\) が一点へ収縮し、他の収縮曲線は無い。
正次元の projective fiber は曲線を含むので、\(f\) はこの section の外で quasi-finite、
また proper birational と normality から \(f_*O_Y=O_X\) で fiber は connected。
従って section を含む fiber に別の孤立点は無く、その fiber は section そのもの、
それ以外では proper quasi-finite birational と normality から同型である。
従って \(r=1\) の \(f\) は small contraction で、exceptional locus は codimension 二。
section ring の degree 一生成、\(F\) の tautological embedding、又は特異 model の分類は使っていない。

## 3. \(r=2\)：line vertex から次数五支持は不可能

この場合 \(F=M\) は degree 五 line bundle である。
実 \(D\) への五 section の空間 \(W\subset H^0(E)=\mathbb C^2\oplus H^0(M)\) は
\(Y\) で basepoint-free である。
その \(\mathbb C^2\) への restriction は rank 二でなければ、\(B\times\mathbb P^1\) に basepoint を持つ。
従って \(W\cap H^0(M)\) は三次元で、各 \(b\in B\) でも同時に零とはならない。
もし同時に零なら、その \(\mathbb P^2\) fiber で残る二 section の係数を零にして basepoint を作れるからである。

各 \(\mathbb P^2\) fiber の image は、共通 vertex line と
三 section による plane curve の点を結ぶ plane である。
従って実像 \(D\) はその line を vertex とする cone。
これは全 section ring の degree 一生成を仮定せず、実五 section から直接従う。
vertex の外の局所環は plane base curve の局所環と affine space の積である。
従って二次元非正規支持は、base curve の singular point と vertex line の join、すなわち plane の有限和。
vertex line 自体は一次元であり、他の二次元支持を生まない。
次数五 integral \(S\) は存在しない。
ここで rank 零用 high-degree-cone 排除は使っていない。

## 4. \(r=0,1\)：全導手と閉留数一形式

以後 \(f\) は同型又は small で、divisor を一つも縮めない。
有限双対性から \(\omega_D\simeq O_D\) の generator の \(X\) 上の rational pullback は
divisor \(-C_X\) を持つ。
さらに \(Y\) への pullback の divisor は、extra exceptional Weil divisor を持たず
\(-\overline C\) である。projective bundle の canonical formula から

\[
 \overline C\sim-K_Y=3H-p^*\det E. \tag{2}
\]

elliptic \(B\) では \(K_B\simeq O_B\)、\(\deg\det E=5\) である。
\(\operatorname{Pic}Y=\mathbb ZH\oplus p^*\operatorname{Pic}B\) であり、
\(\operatorname{Pic}^0B\) の twist を数字二つの Pic 群と取り違えない。

\(\delta_S=1\) の全 generic transverse type は node 又は cusp。
cusp なら \(\overline C=2\Gamma\) となり、(2) の \(H\) 係数三に反する。
従って node、\(\overline C\) は reduced であり、一つの prime（connected non-split cover を含む）又は
二つの prime の全場合が残る。

node 完成環の各 minimal prime は Hamiltonian derivation に保存され、各 branch は環の surjective 商。
従って normalization branch の Poisson bracket は regular で、導手 ideal に tangent である。
係数体が Poisson central とは仮定しない。
有限分離 residue extension・completion からの descent を含む全局所論証は、既受理の
`RANK_TWO_SEXTIC_SURFACE_SEPARATE_REVIEW.md` §3 の node 補題と同じである。

ambient inverse symplectic 二形式の closed Poincaré residue を \(\alpha\) とする。
\(\alpha=c\,\iota_{\pi_X}\Omega_X\)、\(c\ne0\) の inverse-Pfaffian 恒等式より、
\(Y\) 上の rational 一形式 \(\alpha\) は閉じており、全 divisorial pole は \(\overline C\) にある。
導手一般点で \(\Omega_Y=(du/u)\wedge ds\wedge dt\) と取れば、
tangency により mixed normal bracket は \(u\) で割れ、

\[
 \alpha=\lambda\,du/u+\text{regular}.
\]

閉性から各 prime 上の residue \(\lambda\) は定数であり、
\(\pi|_S\) の一般 rank 二と characteristic zero の generic separability から \(\lambda\ne0\)。
small exceptional section は codimension 二で、追加の divisorial pole を作らない。
locally free \(\Omega_Y^1\) の Hartogs 延長も使用出来るが、
以下の矛盾は非零 residue だけで得るので \(H^0(\Omega_Y^1)=0\) を仮定しない。
elliptic bundle では実際にこの holomorphic 一形式空間は非零である。

## 5. 全 ample 切断の residue theorem で single と split を排除

数値 fiber class を \(F\) と書く。
\(H^3=5\)、\(H^2F=1\)、\(F^2=0\)。
\(B\) 上の degree 一 line bundle \(N\) を選び、
\(A_n=H+p^*N^{\otimes n}\) とする。
\(n\gg0\) では \(A_n\) は very ample であり、数値 class は \(H+nF\)。
一般二 section の intersection curve は smooth で、small exceptional section、
\(\overline C\) の singular locus、全 local-residue bad locus という codimension 二集合を避けて取れる。
従って \(\alpha\) の curve 上の全 pole は導手との transverse intersection だけである。
compact curve の residue theorem をこの全 curve に適用する。

single prime なら \(\overline C=\Gamma\)、(2) より

\[
 \Gamma A_n^2=(3H-5F)(H+nF)^2=10+6n.
\]

従って \((10+6n)\lambda=0\)、\(\lambda\ne0\) に矛盾する。
non-split node で二つの局所 branch が一つの大域 prime を作る場合も、同じ一つの residue を数える。

split なら \(\overline C=\Gamma_1+\Gamma_2\)、各 \(\Gamma_i\to S\) は degree 一なので
\(H^2\Gamma_i=5\)。Pic 群を正確に保ち

\[
 \Gamma_i\sim a_iH+p^*M_i,\qquad b_i=\deg M_i
\]

と書く。一般 \(\mathbb P^2\) fiber への restriction より \(a_i\ge0\)。
\(a_i=0\) の prime は fiber 自身だが、その \(H^2\)-degree は一なので不可能。
従って \(a_i\ge1\)、\(a_1+a_2=3\)、\(5a_i+b_i=5\) より、交換を除き

\[
 (a_1,b_1)=(1,0),\qquad(a_2,b_2)=(2,-5).
\]

ここでは \(M_1\) を trivial とせず、すべての \(\operatorname{Pic}^0\) twist を保持している。
全交点数はそれらの twist によらず

\[
 \Gamma_1A_n^2=5+2n,\qquad\Gamma_2A_n^2=5+4n.
\]

従って、十分大きい相異なる二つの連続整数 \(n,n+1\) について residue theorem は

\[
 \begin{pmatrix}5+2n&5+4n\\7+2n&9+4n\end{pmatrix}
 \begin{pmatrix}\lambda_1\\\lambda_2\end{pmatrix}=0.
\]

この行列の determinant は正確に 10、従って
\(\lambda_1=\lambda_2=0\)。両 residue の非零性に矛盾する。
一つの固定切断の residue 和だけでは得られない矛盾を、二つの ample class が与える。

以上で \(r=0,1\) の cusp/node、single/non-split/split を全て排除し、
§3 の \(r=2\) と合わせて elliptic normalization 全枝の不存在が証明された。

## 6. 保存する補集合と証明の状態

本命題は genus 一の elliptic polarized model 全体の空性であり、一般 bundle だけの上界ではない。
unstable/stable の区別や target projection の標本、normality・smoothness の追加仮定、
未読の elliptic scroll 分類、GKKP、rank 零 cone 排除を使っていない。
canonical index 二 model、genus 零の別処理、その他の次数や全 E は本ノートから閉じない。
既受理入力の rank 非依存部分、新しい全収縮 locus、二 ample class の留数矛盾を区別し、
独立監査前の producer 証明として保存する。global 台帳・main・Git は変更しない。
