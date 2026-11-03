# Qstar の真の四次初期式と、残存する二つの15平面

2026-11-03。\(Q_\star=Q_{2/9}\) の実際の35変数有限モデルにおいて、
\(u=a-2/9,z_0,\ldots,z_{33}\) を全て ordinary degree 1 とする。
\(I\subset\mathbb Q[[u,z_0,\ldots,z_{33}]]\) は full350 関係からなる実際のイデアル、
\(I^*=\operatorname{in}_{\mathfrak m}(I)\) はその ordinary 初期イデアルとする。
補助的な古い重みと ordinary degree は区別する。

## 1. 結論と限界

実際の full350 式を ordinary degree 4 まで厳密に計算し、独立な物理的 Jacobi 計算でも
全350本の係数を確認した。次に低次項の全ての cancellation を含めて真の四次初期式を
作った。その結果、先に構成された **二つの15次元線形平面は、真の初期イデアルの四次以下の
全ての homogeneous component 上でも消える**。

従って

\[
 K_{\leq4}=((I^*)_2,(I^*)_3,(I^*)_4)
\]

について \(\dim\mathbb Q[u,z]/K_{\leq4}\geq15\)。どの12本の同次線形式を加えても商の
次元は3以上であり、四次以下だけを用いた12-linear-section Artinian 証明は不可能である。

これは **実際の完成局所環や全初期イデアルの次元が15以上という主張ではない**。
15次元の実際の Poisson 枝を発見したという主張でもない。五次以上の初期式、別の構造的
議論、または実際の枝の構成が必要である。今回、五次 jet や大きい GB は開始していない。

## 2. 入力と二つの線形平面

使用する全350関係は、独立検証済みの正確な有限35変数 DAG
`reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/model.json.gz` である。
これは低次近似だけのモデルではない。既存の二次・三次 packet は57個の独立な二次初期式と
17個の新しい三次 class を与え、その生成するイデアルを \(J\) と書く。

二つの平面 \(L_0,L_1\simeq\mathbb A^{15}\) は
`reports/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63/certificate.json.gz` の
`branches[*].all35_graph` を homogenize したものである。両者の自由変数は

\[
 u,z_5,z_{13},z_{14},z_{15},z_{16},z_{17},z_{22},z_{24},z_{25},z_{26},z_{27},z_{28},z_{29},z_{31}.
\]

残る座標はゼロとし、次のものだけを追加する。

\[
 z_{19}=-(z_{27}+z_{29})/2,\qquad z_{20}=-2z_{31},
\]

\[
 \begin{array}{c|cc}
 &z_6&z_3\\\hline
 L_0&0&-z_{14}+2z_{27}\\
 L_1&5z_{27}/7&-z_{14}+4z_{27}/7.
 \end{array}
\]

入力の affine chart は \(z_{22}=1\) である。そこに現れる唯一の定数1を \(z_{22}\) に
置き換えると、この線形 graph が得られる。今回の coverage producer は全350二次式と
下記の全293三次 tail が両平面上で消えることを改めて直接確認した。

\(\eta=z_{33}\neq0\) では真の接錐は正確に T39 の滑らかな9次元接平面であることが
別 packet で証明済みである。今回の15平面は \(\eta=0\) 上にあり、その結果と矛盾しない。

## 3. 全350式の正確な四次 jet

各関係を

\[
 F_i=Q_i+C_i+D_i+O(\mathfrak m^5)
\]

と書く。定数・一次項は全てゼロである。
`q_ordinary_quartic_initial_20261103_v63.py` は有理数係数の truncated polynomial algebra
を ordinary degree 4 まで実装する。\(a\) の単元逆数は非零定数の周りの幾何級数を四次まで
取り、各積が1となることを確認する。DAG の一時値は最終使用後に解放する。

全350本の \(Q_i,C_i\) は frozen packet と完全一致した。復元された物理的 full150
係数から、5変数の exponent dictionary の微分・積だけで divergence と Jacobi を再計算し、
全350本の四次 jet が一致した。DAG 計算75.4秒、cofactor/checkpoint 込み112.8秒、物理的
再検証まで201.6秒、peak RSS 約1.334 GB、AS cap 3 GB で完了した。全weighted高次数の
1800万項を再展開したものではない。

## 4. 定数二次消去を四次まで運ぶ

既存の正確な定数行操作で350本を、57本の二次独立な実際の関係 \(A_i\) と、二次項がゼロの
293本の実際の関係 \(G_j\) に取り替える。この取替えは可逆な定数行操作であり、実際の
生成イデアルを変えない。

\[
 A_i=Q_i+C_i+D_i+\cdots,\qquad
 G_j=\widetilde C_j+\widetilde D_j+\cdots.
\]

同じ frozen 操作を四次 tail と actual350 に対する定数 cofactor に正確に適用した。
293本の \(\widetilde C_j\) のうち158本が恒等的ゼロである。その158本の
\(\widetilde D_j\) のうち非零なのは **raw \(F_{34}\) 自身の2790項の四次式だけ**であり、
これは両15平面上で恒等的ゼロであった。番号は全て0-basedである。

これら158行について actual350 の cofactor を直接 replay し、二次・三次項がゼロ、四次項が
保存した tail に一致することを確認した。raw quartic tail を無条件に初期式と呼んでいない。

次に非零の135本の \(\widetilde C_j\) だけの定数行空間を調べた。古い重みごとの厳密 QQ
rank の合計は80、定数 kernel は55次元である。全55個の実際の cofactor の四次項は
**全変数で恒等的ゼロ**であった。これは四次までの identity であり、全次数での syzygy と
いう主張はしない。

## 5. 全ての三次 cancellation と四次 carry

三次式の行列

\[
 T:S_1^{\oplus57}\oplus\mathbb Q^{293}\longrightarrow S_3,
 \quad(\ell_i,b_j)\longmapsto\sum_i\ell_iQ_i+\sum_jb_j\widetilde C_j
\]

を用いる。行数は \(35\cdot57+293=2288\)、既存の正確な QQ rank は1911であり、
kernel の次元は377である。

まず \(35\cdot57=1995\) 本の線形倍 \(x_vQ_i\) だけの kernel を計算した。その次元は101。
実際の cofactor は \(\sum h_{iv}x_vA_i\) であり、三次項がゼロとなることを全35変数で
確認した。従ってその四次 carry は

\[
 \sum h_{iv}x_v C_i.
\]

101本全ての carry が両15平面上でゼロとなった。

さらに1995本の線形倍と135本の非零 \(\widetilde C_j\) を同時に入れた行列の kernel を
全て計算した。その次元は219であり、上の101方向、定数三次 kernel 55方向、および混合する
63方向を含む。古い重みごとに rational row の分母を払い、FLINT の exact integer
nullspace を用いた。native code は kernel 次元を既存の QQ rank と照合し、basis 行列の
rank と元行列との積がゼロであることも確認する。計算は約0.6秒であった。

219本の全てについて、actual350 に対する **定数＋一次の cofactor polynomial** を保存した。
この cofactor を original \(Q_i,C_i\) に適用し、二次項・三次項が全て恒等的ゼロとなる
ことを別途直接確認した。よって四次 carry は真の四次初期式またはゼロであり、その式は

\[
 \sum h_{iv}x_vC_i+\sum k_j\widetilde D_j
\]

で与えられる。全219本が両15平面上で厳密にゼロとなった。特殊な数値点だけの確認ではない。
残る158個の zero-cubic row を合わせ、377次元 kernel の全体を網羅している。

## 6. 四次 component を漏らしていない理由

任意の \(H\in I\) を \(A_i,G_j\) の formal coefficient を用いて書く。
\(\operatorname{ord}H\geq4\) なら、二次項に寄与する \(A_i\) の定数係数は、独立な57本の
\(Q_i\) により全てゼロとなる。三次項の cancellation は正確に \(\ker T\) である。
四次項には、これに加えて \(A_i\) の二次係数と \(G_j\) の一次係数が寄与するだけなので、

\[
 (I^*)_4=operatorname{span}_{\mathbb Q}
 \{S_2Q_i,\ S_1\widetilde C_j,\ \text{quartic carry of }\ker T\}.
\]

全 \(Q_i\) と全 \(\widetilde C_j\) は両平面上でゼロであり、\(\ker T\) 全377方向の carry
もゼロである。従って四次 component 全体が両平面上で消える。四次 component の rank や
Hilbert function は今回求めていないが、この包含は全体についての証明である。

\(L_0,L_1\subset V(K_{\leq4})\) よりその次元は15以上。homogeneous な12線形形式による
切断では次元は高々12しか減らないため、切断の次元は3以上となる。

## 7. 証拠、再現、次の有限問題

保存先は `reports/Q_ORDINARY_QUARTIC_INITIAL_20261103_v63/`。

- `report.json`, `quartic-tails.json.gz`: full350 の四次 tail、全293 cubic/quartic carry と
  actual350 定数 cofactor、full physical four-jet replay の結果。
- `plane-restrictions/`: 158 zero-cubic directions の restriction。
- `constant-kernel/`: 非零135 cubic tails の55定数関係と zero quartic identities。
- `linear-syzygies/`: 101 linear-quadratic syzygies と全 cofactor recipe。
- `mixed-syzygies/`: 全219 kernel directions、actual350 cofactor polynomial、両平面への
  全 restriction。native の整数入力・出力は圧縮保存し、実行 binary は保存しない。
- `coverage.json`: 入力 hash の整合、全350二次式・全293三次式の平面上の直接ゼロ確認、
  377方向の被覆と結論の範囲。

全 producer は `code/q_ordinary_quartic_*_20261103_v63.py` と同名 native `.cpp`。
新しい出力先に、まず full4 jet を作り、その出力を後続の `--quartic` に指定すれば再現できる。
native nullspace には既存 Cloud 設定を用いる。

```sh
python code/q_ordinary_quartic_initial_20261103_v63.py \
  --output work/quartic-replay/jets --cpu-limit 240
source /workspace/cloud-setup/activate.sh
python code/q_ordinary_quartic_mixed_syzygies_20261103_v63.py \
  --quartic work/quartic-replay/jets --output work/quartic-replay/mixed
```

最初の四次 jet を理由なく繰り返す必要はない。保存された exact tail と cofactor から、
後続の kernel・restriction 証明だけを約15秒で再現できる。追加の constant/linear-only
producer はこの完全な219方向の前に試した bounded 方法の receipt でもある。

次に必要なのは、五次以上の **真の cancellation を含む**関係、または別の構造的 argument
である。raw quintic tail だけを加えることや、同じ四次 ideal に大きい GB/12-linear-section
を繰り返すことは、この障害を解決しない。今回の結果を実際の枝の次元下界に取り違えない。
