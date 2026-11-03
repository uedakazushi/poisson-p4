# 真の ordinary 初期イデアルの η≠0 chart

2026-11-03。\(Q_\star=Q_{2/9}\) の35変数有限モデルにおいて
\(u=a-2/9,z_0,\ldots,z_{33}\) を全て ordinary degree 1 とする。
\(I\) は実際の full350 残余関係からなる完成局所イデアル、
\(I^*=\operatorname{in}_{\mathfrak m}(I)\) はその全 ordinary 初期イデアルとする。

既に証明済みの57個の二次初期式と17個の新しい三次初期式で生成される homogeneous
ideal を \(J\) と書く。従って \(J\subset I^*\) である。\(J\) が \(I^*\) 全体である
とは仮定しない。ここでは \(\eta=z_{33}\neq0\) の chart を厳密に決定する。

## 1. 証明した結果

下記の26本の線形関係が定める素イデアルを \(P\) とすると、

\[
 J[\eta^{-1}]=I^*[\eta^{-1}]=P[\eta^{-1}],
\]

かつ

\[
 (\mathbb Q[u,z_0,\ldots,z_{33}]/P)[\eta^{-1}]
 \simeq\mathbb Q[\eta,\eta^{-1},z_{22},z_{24},z_{27},z_{28},z_{29},z_{30},z_{31},z_{32}].
\]

従って **実際の ordinary 接錐の η≠0 の部分は被約かつ滑らかな9次元**であり、
\(\mathbb G_m\times\mathbb A^8\) である。これは既知の滑らかな T39 成分の接空間と一致する。
真の接錐の他の既約成分は全て \(\eta=0\) に含まれる。

これは根基だけの一致ではない。元の \(J\) の localization 自体が線形素イデアルとなる。
一方、**実際の完成イデアル \(I[\eta^{-1}]\) が線形イデアルであるとは主張しない**。
ordinary 初期イデアルと実際のイデアルを区別し、接錐の全体や全 Poisson 枝の分解も
今回の結論に含めない。

## 2. 三つの線形関係は根基を取らずに得られる

\[
 b=z_{11},\quad c=z_{12},\quad d=z_{23},\quad\eta=z_{33}
\]

と書く。57個の二次式の **定数 \(\mathbb Q\) 行空間**と \(\mathbb Q[b,c,d,\eta]_2\) の
共通部分は厳密に9次元であり、次を基底に持つ。

\[
\begin{gathered}
 b^2-16\eta^2,\ bc-24\eta^2,\ bd+16\eta^2,\ b\eta-4\eta^2,\quad
 c^2-36\eta^2,\\
 cd+24\eta^2,\ c\eta-6\eta^2,\ d^2-16\eta^2,\ d\eta+4\eta^2.
\end{gathered}
\]

これは直線 \((b,c,d)=(4\eta,6\eta,-4\eta)\) 上で消える全ての二次式の空間である。
その生成する4変数イデアルは

\[
 (b-4\eta,c-6\eta,d+4\eta)(b,c,d,\eta)
\]

に等しい。特に三つの線形式の平方だけでなく

\[
 \eta(b-4\eta),\qquad\eta(c-6\eta),\qquad\eta(d+4\eta)
\]

が既に \(J\) に入る。よって \(\eta\neq0\) では三つの線形式そのものが \(J\) に入る。
nilpotent direction を人為的に除いた quotient に置き換える必要はない。

証拠ファイルには基底9本と平方3本の **全ての定数 cofactor** を保存した。例えば
\(Q_k\) を full350 残余式の0-based 番号 \(k\) の二次初期式とすると、

\[
 \eta(b-4\eta)=
 -\frac{41}{912}Q_{126}+\frac12Q_{224}-\frac{137}{456}Q_{310}
 -\frac{16}{1045}Q_{204}-\frac{29}{1045}Q_{10}
 -\frac{12}{95}Q_{47}-\frac{12}{1045}Q_{92}.
\]

全ての恒等式は630個の二次単項式について \(\mathbb Q\) 上で直接確認した。

## 3. η=1 chart の正確な polynomial triangular elimination

homogeneity により \(S[\eta^{-1}]\) は、\(\eta\) と比座標 \(x_i/\eta\) の Laurent
polynomial ring と同一視できる。従って \(\eta=1\) でのイデアル計算を行い、最後に
homogenize すればよい。まず既に \(J\) にある三つの線形関係で
\(b=4,c=6,d=-4\) とする。

残る74本の式に対して、選んだ式が

\[
 f_j=c_i(x_i-r_i),\qquad c_i\in\mathbb Q^*,\quad r_i\text{ は }x_i\text{ を含まない}
\]

となる23段の消去を行う。各段で全ての式について

\[
 f_k-f_k|_{x_i=r_i}=(x_i-r_i)D_{ki}
\]

を polynomial dictionary の積と差で確認した。従ってこれは単なる点での代入検査ではなく、
各段での **生成イデアルの厳密な等式**を与える。全ての pivot は非零有理数であり、
変数に依存する式を invert していない。GB は用いていない。

23段の後には全74式がゼロになり、自由変数は

\[
 z_{22},z_{24},z_{27},z_{28},z_{29},z_{30},z_{31},z_{32}
\]

の8個だけとなる。全ての中間 RHS を後ろから代入し直すと、最終 graph は次の
**線形**関係で与えられる。右辺に現れる8変数と \(\eta\) は自由である。

| 変数 | 最終 RHS | 変数 | 最終 RHS |
|---|---|---|---|
| \(u\) | \(0\) | \(z_0\) | \(-(z_{22}+z_{32})/3\) |
| \(z_1\) | \(-5z_{31}/6\) | \(z_2\) | \(z_{22}\) |
| \(z_3\) | \((4z_{27}+z_{29})/3\) | \(z_4\) | \(0\) |
| \(z_5\) | \(5z_{24}/2+2z_{28}/3\) | \(z_6\) | \(z_{27}+z_{29}/2\) |
| \(z_7\) | \(-2z_{30}+4z_{31}/3\) | \(z_8\) | \((z_{22}+10z_{32})/3\) |
| \(z_9\) | \(z_{30}/4\) | \(z_{10}\) | \(z_{32}\) |
| \(z_{11}\) | \(4\eta\) | \(z_{12}\) | \(6\eta\) |
| \(z_{13}\) | \(-8z_{28}/9\) | \(z_{14}\) | \(-2z_{27}/3-z_{29}\) |
| \(z_{15}\) | \(0\) | \(z_{16}\) | \(4z_{28}/3\) |
| \(z_{17}\) | \(5z_{29}\) | \(z_{18}\) | \(3z_{30}-2z_{31}\) |
| \(z_{19}\) | \(-z_{29}/4\) | \(z_{20}\) | \(-2z_{31}/3\) |
| \(z_{21}\) | \((z_{22}+z_{32})/3\) | \(z_{23}\) | \(-4\eta\) |
| \(z_{25}\) | \(0\) | \(z_{26}\) | \(4z_{28}/3\) |

この26本の線形式は異なる26変数について monic な graph をなすので、生成イデアル
\(P\) は素イデアルである。生成式の逆向きの確認として、全74式に最終 graph を代入し、
全係数ゼロを確認した。さらに、選ばれなかった式も含め **全350本の二次初期式**でも
同じゼロを確認した。以上より \(J[\eta^{-1}]=P[\eta^{-1}]\) が証明された。

## 4. 既知の実際の T39 quotient との一致

ここでは次の既存の独立検証済み入力だけを再利用する。

- `reports/Q_LOCAL_Q51_QUOTIENT_20261103_v63/T39-nine-parameter-model.json.gz`
  は exact polynomial \(\mathbb Q^9\) モデルを持ち、実際の full350式と Q51 の45式を
  全て満たす。自由な親座標9個への projection が厳密な left inverse であり、従って
  completed full35 ring から正則9次元環への全射を与える。
- 自由座標の順は
  \((z_{22},z_{24},z_{27},z_{28},z_{29},z_{30},z_{31},z_{32},z_{33})\) である。
- 同じ packet の `certificate.json` の `Q51_jacobian_45x35` は26階数の線形関係行列である。
  ここで列は新しい Smith frame でなく、今回と同じ元の
  \((a-2/9,z_0,\ldots,z_{33})\) の順である。

今回の producer は、その7527-node DAG を \(\mathbb Q\) の ordinary 1-jet algebra で
独立に評価した。得られた35行9列の線形 map の **全315成分**が、上表の graph の係数と
一致する。9個の自由座標の行は単位行列であり、rank は9。45行35列の Q51 行列の rank
は26で、この map を消すことも直接確認した。従って T39 の接平面の素イデアルは正確に
\(P\) である。

実際の T39 quotient への代入では、\(I\) の全ての式がゼロになる。式の最初の ordinary
homogeneous term に注目すれば、\(I^*\subset P\) である。よって

\[
 J\subset I^*\subset P
\]

であり、既に証明した \(J[\eta^{-1}]=P[\eta^{-1}]\) から
\(I^*[\eta^{-1}]=P[\eta^{-1}]\) が従う。従って \(V(P)\) は真の接錐の1成分であり、
\(\eta\neq0\) と交わる唯一の成分である。その open には nilpotent や embedded
structure もない。残る接錐の成分・構造は \(\eta=0\) に限定される。

ここでは過去の full350 nonlinear T39 identity を理由なく再計算していない。今回新規に
確認したものは全35行9列の first-order map、線形行列の階数と消去、および今回の
\(J\) との正確な chart ideal equality である。既存 nonlinear identity の証拠と
independent review は今回の input hash で固定した。

## 5. 再現と保存範囲

新規 producer は `code/q_true_initial_eta_open_20261103_v63.py`。約3秒、通常の
Python/SymPy のみで実行でき、GF による推測や巨大な GB を使用しない。

```sh
python code/q_true_initial_eta_open_20261103_v63.py \
  --output work/q-true-initial-eta-v63/replay
```

`reports/Q_TRUE_INITIAL_ETA_OPEN_20261103_v63/` に保存する。
`certificate.json.gz` は定数 cofactor 12本、準備後の74生成式、23段全ての選択生成式・
非零 pivot・全74行の division cofactor、完全に back-substitute した graph、26線形式、
T39 の fresh 1-jet map を持つ。`report.json` は入力 hash と主張の範囲を固定する。
全ての式番号と変数番号は0-based である。

\(\eta=0\) 上の \(J\) は、より大きい次元の false initial supports を持ち得る。
今回の9次元 chart を global dimension upper bound 9や12に置き換えてはならない。
高次初期式によるその境界の縮小と、実際の全完成局所環の枝分解は別の未解決問題である。
