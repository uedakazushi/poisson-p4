# Full V56 の actual m(0) rank24 と正則比較 operator

2026-11-04、v63。新しい source solve・jet を用いず、保存済み actual Artin20 値から
full 比較行列 (V:3\times56\) の canonical cofactor (m\) の中心を求めた。

\[
 \operatorname{rank}_{\mathbb Q}m(0)=24, \qquad
 \det M(0)=-\frac{397267578125}{50096498540544}\ne0.
\]

固定した24列、その24×24行列と逆行列、全56列の数値因子分解
\(M(0)\Lambda=m(0)\) を保存して全成分で検査した。
有限値の計算は0.112秒・約29 MB、証拠は150 KB、typed operator は2502 nodesである。

## 正確な source と係数抽出

(V\) の最初48列は existing (Q48\) の \((\tau_A,\tau_B,\tau_A\tau_B)\) の三行であり、
四つの tensor 乗数の順をそのまま保つ。残り八列は actual honest C0 の odd 像の負号、
すなわち \(-z_{A,0},\ldots,-z_{A,6},-g_A\) が最初の target に入る。
定数 target の unit1 を消した正確な source 順序である。

(Q48\) は既存105解の十二組の \(\delta,b_A,b_B\)、両 \(\nu\) から四成分代数の
十六乗法則で評価した。元の105方程式は再実行していない。
actual source scalar 関数は、先行 typed chainlift の same- (R_0\) 関数参照で束縛する。

実 (H_{20}\) から \(\mathbb Q[r]/r^5\) への二つの map

\[
 z_i\mapsto0,\qquad z_i\mapsto-c_i r^2
\]

を、全20基底・全800積恒等式で検査した。
各 sector の関数 (f_s\) について、最初の七 cofactor の中心は
actual (H_{20}\) 表示の一次 (z_i\) 係数である。
三次の basis alias から一次係数を推測していない。
最後の cofactor は対応する curve map 上の (r^{n_s}\)、\(n_s=3,3,2\) の係数である。
全56列について、その次数未満の係数が全て零であることも検査した。

最後の actual 関数 (g_A\) には未指定の四次 tail があるが、
\(g_A=(125/24)r^3+O(r^4)\) であり、その tail は全ての (m(0)\) slot に寄与しない。
証拠では最後の有限値を「中心抽出に十分な cubic representative」と明記した。
実 (H_{20}\) 関数の四次係数が零だとは主張していない。
後述の typed 関数参照は actual (g_A\) 自体である。

## 全次数の typed m,n

各 actual (V\) 成分を、実 sector 座標の七つの ordered-prefix Hadamard と、
exact な (r^{n_s}\) 除法により分解する有限 operator を保存した。
それぞれの関数は先行 chainlift の実 source operator を \(\ell=0\) へ正則に制限したものである。
閉点で非単位となる元の (r/g\) を新しく割り直さない。

全次数の割り切れは、既受理の exact Cbar augmentation と実三 cyclic sector 表示から従う。
有限値の零だけから高次の割り切れを推測していない。
この operator で (d_1m=V\) が成り立つ。
選んだ24列 (M\) は上の中心行列式により formal unit matrix であり、その有限逆行列 operator を定義する。
source の24列への埋込を (E_S\) として

\[
 n=E_S M^{-1},\quad mn=I_{24},\quad Vn=d_1, \quad h=0
\]

が全次数で成り立つ。scalar の全 Taylor 展開や348環への係数評価を行ったという主張ではない。

full source の (X\) は、chainlift の \((B48,P8)\) の640列と、
complete normalizer \((\operatorname{Syz}Q48,0)\) の128列である。
除いた (P0\) は actual (UP=QB\) の定数行によって正確に整合する。
各 cycle (mX_j\) を literal Koszul (d_2\) に持ち上げれば、full段階で116生成・936関係の
正確な (T\) 表示が得られる。proper adapter の16 free-zero unit の消去後は100/920になる。

今回その768 cycle の係数持上げは行っていない。
actual (T\) の長さ、(K\) の最小生成数・torsion、元の深さは未計算である。

主要証拠：
`reports/AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63/actual-fullV56-m0-and-typed-comparison.json`。

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_fullV56_comparison_unit_20261104_v63.py
```
