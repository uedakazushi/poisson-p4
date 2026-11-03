# 実際の blowup による二つの15平面の一般方向の除去

2026-11-03。対象は \(Q_\star=Q_{2/9}\) の actual full35 完成局所 Poisson モデルである。
ordinary 変数を \(x=(u,z_0,\ldots,z_{33})\)、\(u=a-2/9\) と書く。
\(F_1,\ldots,F_{350}\) はその実際の関係、\(I=(F_i)\)、\(I^*=\operatorname{in}_{\mathfrak m}I\)
とする。全 \(F_i\) は二次以上であり、\(F_i=Q_i+C_i+D_i+\cdots\) と展開する。

## 1. 証明したことと範囲

`Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63` の二つの15次元線形平面 \(L_0,L_1\) は、
真の ordinary 初期式の四次以下を全て満たすことが別 packet で確認された。しかし今回、
**各平面全体が実際の接錐の台に含まれることはない**と証明した。従って

\[
 \dim\bigl(L_b\cap\operatorname{Supp}(\operatorname{Spec}S/I^*)\bigr)\leq14,
 \qquad b=0,1.
\]

両平面の一般方向は実際の接錐から除かれる。これは五次以上の真の初期条件が両平面の一般点を
除くことの存在証明でもあるが、具体的な高次 homogeneous polynomial は今回展開していない。

全ての平面上の点が除かれるという主張ではない。後述の有理点 \(p_0\) 自体が除かれるとも
主張しない。また **full35 局所環全体の次元上界14または12ではない**。二つの特定の15平面に
ついての、実際の blowup と saturation を用いる除去である。

ordinary な \(\varepsilon^2\) correction が存在しないことだけでは除去できない。
\(\varepsilon^{3/2}\) 等の Puiseux normal directions もあり得る。今回の証明はそのような
方向を禁止する仮定を置かず、actual strict transform の有限性と次元を用いる。

## 2. 実際の blowup chart

二つの平面の同次 graph は
`reports/Q_TRUE_INITIAL_BOUNDARY_CONES_20261103_v63/certificate.json.gz` に厳密に保存され、
`research/notes/Q_ORDINARY_QUARTIC_INITIAL_20261103_v63.md` §2にも記載されている。
\(z_{22}\neq0\) の projective chart における自由な14変数を

\[
 p=(u,z_5,z_{13},z_{14},z_{15},z_{16},z_{17},z_{24},z_{25},z_{26},z_{27},z_{28},z_{29},z_{31})
\]

とし、その方向 graph を \(v_b(p)\in\mathbb A^{35}\)、\((v_b)_{z_{22}}=1\) とする。
graph で決まる残る20座標に独立な normal 変数 \(n\) を加える。全35変数空間の原点 blowup の
この chart は正確に

\[
 z_{22}=r,\qquad x=r\bigl(v_b(p)+n\bigr)
\]

である。\(n\) は自由な14座標と \(z_{22}\) の位置には成分を持たない。従って変数数は
\(14+20+1=35\) であり、この変換は全 ambient blowup chart の座標を与える。

全ての \(F_i\) が二次以上なので、weak transform

\[
 G_i(p,n,r)=r^{-2}F_i\bigl(r(v_b(p)+n)\bigr)
\]

は正則な formal power series である。元モデルの \(a\)-only rational denominators は
\(a=2/9+r p_u\) で依然単元である。展開は

\[
 G_i=Q_i(v_b+n)+r C_i(v_b+n)+r^2D_i(v_b+n)+\cdots
\]

となる。これは選んだ低次式だけの人工的なモデルでなく、実際の full350 weak transform
である。原点以外の strict transform は \((G_i):r^\infty\) で定まる。

## 3. 両平面上の同じ有理点での exact certificate

両方について

\[
 p_0=(1,2,2,2,2,5,5,7,2,1,5,3,5,6)
\]

を、上に記した14変数の順で取る。\(Q_i(v_b(p_0))=0\) であり、\(n=r=0\) は weak
transform の点である。

\(p=p_0\) を固定した21変数 \((n_{20},r)\) における全350式の線形行列は

\[
 \mathcal L_b=[dQ_{v_b(p_0)}|_{n_{20}}\ \ C(v_b(p_0))]
\]

である。厳密な有理数計算により両方とも階数19であり、選択した19個の unknown は \(r\) を
含む。正規化した2次元 kernel の \(r\) 成分は厳密にゼロである。

従ってこの kernel 上の \(G_i\) の二次項は、**正確に \(Q_i\) 自身を kernel に代入したもの**
である。\(rC_i\) の二次項や \(r^2D_i\) は kernel の \(r=0\) により消える。
五次以上の式、または full quartic tail を無視した近似に依存していない。

19個の pivot 式を formal implicit function theorem で消去する。\(K\) を正規化した
21行2列 kernel matrix、\(B\) を \(Q_i(K_1),Q_i(K_1+K_2)-Q_i(K_1)-Q_i(K_2),Q_i(K_2)\)
の350行3列行列とする。選択した pivot 行を \(R\)、pivot 列を \(C\) と書けば、残る式の
二次係数行列は厳密に

\[
 B'=B-\mathcal L_{:,C}(\mathcal L_{R,C})^{-1}B_R.
\]

この行列の階数は両方とも3である。従って残る2変数の二次式はその maximal ideal の平方全体を
張る。transverse な閉 fiber の formal local ring は長さ3であり、basis は \(1,s,t\)、
\((s,t)^2=0\) となる。高次項を追加してもこの有限性の上界を弱めない。線形階数が19なので
embedding dimension は2であり、実際に長さはちょうど3である。

保存した非零小行列式は次の通りである。

| 平面 | linear 19次 minor | implicit quadratic 3次 minor |
|---|---|---|
| \(L_0\) | \(-956393211332282501953125/100352\) | \(-27221241815728515625/311328730570208324171328\) |
| \(L_1\) | \(-18169764766174127197265625/205641981272\) | \(9471284996723/21652080156000\) |

全350行の行列、選択行・列、正規化した kernel、補正前後の二次行列は証拠ファイルに保存した。
どの値も有限体標本や浮動小数からの推測ではない。

## 4. 相対有限性と saturation の次元論

\(R=\mathbb Q[[p-p_0]]\)、\(\dim R=14\) とする。weak transform のこの点での completion を

\[
 S=R[[n_{20},r]]/(G_1,\ldots,G_{350})
\]

と置く。19個の pivot の消去後、\(S\) は \(R[[s,t]]\) の quotient となる。前節の証明に
より \(S/\mathfrak m_R S\) は長さ3である。formal finiteness theorem
（complete Nakayama / Weierstrass）により、\(S\) は \(R\) 上有限であり、高々3個の module
generator を持つ。ここでは \(R\)-flatness や generic rank 3を仮定していない。

actual strict transform の completion は

\[
 T=S/(0:_S r^\infty)
\]

である。Noetherian 環では principal saturation は有限段で安定する。completion は平坦で
各 \(r^N\) の kernel と可換するので、元 chart での saturation 後に completion しても、
completion 後にこの saturation を取っても同じである。

\(T\) は有限 \(R\)-module の quotient なので依然有限である。\(T=0\) なら除去に問題はない。
\(T\neq0\) の場合、定義により \(r\) は \(T\) 上 nonzerodivisor である。従って

\[
 \dim T\leq14,\qquad \dim T/(r)\leq13.
\]

blowup の exceptional divisor は \(\operatorname{Proj}(\operatorname{gr}_{\mathfrak m}(S_0/I))\)
である（ここで \(S_0\) は元の35変数 regular ring）。従って \(T/(r)\) は実際の ordinary
接錐のこの projective chart の completed local ring である。

もし平面 \(L_b\) の全体が実際の接錐の台に含まれれば、その projective chart は14次元で
\(p_0\) を通るため、この completed local ring の次元は14以上になる。しかし上では13以下
と証明した。これは矛盾である。従って \(L_b\) は実際の接錐に含まれない。

接錐の台は閉集合であり、\(L_b\) は既約な15次元 affine cone であるから、両者の交叉は
proper closed subset、従って次元14以下となる。これが一般方向の除去の意味である。

重要なのは saturation 後の \(r\) の nonzerodivisor 性である。未飽和 \(S\) には
\(r=0,n=0\) の14次元 section が存在してよく、それを実際の接錐の成分とはみなさない。
Puiseux correction の可能性も含めて actual strict transform 全体に次元論を適用している。

## 5. 証拠と残る問題

producer は `code/q_plane_blowup_finite_transverse_20261103_v63.py`、証拠は
`reports/Q_PLANE_BLOWUP_FINITE_TRANSVERSE_20261103_v63/`。

```sh
python code/q_plane_blowup_finite_transverse_20261103_v63.py \
  --output work/q-plane-blowup-v63/replay
```

約2秒、通常の Python/SymPy のみで再現できる。新しい jet 展開、GB、setup を必要としない。
原料の actual finite model と二次・三次係数は既に独立検証されたものを使用し、今回の
19次・3次 minor と全行の implicit correction は新たに厳密計算した。

先行する診断では両平面の複数の有理点で \(\operatorname{rank}dQ=18\)、
\(\operatorname{rank}[dQ,C]=19\) を得たが、それだけでは一般点の除去を主張していない。
normal 20列の polynomial kernel の一次 ansatz（2048式、300未知数）と二次 ansatz
（13606式、2400未知数）も試し、どちらも kernel はゼロであった。二次 ansatz は両平面合計
31.4秒、2 GB/90秒 cap 内で完了し、kernel 次数の試行はそこで止めた。
今回の成功した証明は、一般点での rank18の上界も、そのような ansatz も必要としない。

残るのは両平面の proper closed residual subsets の決定、接錐の他の層、ならびに実際の
full35 局所環の全枝・次元である。今回の方法を full model の global upper bound と
取り違えない。
