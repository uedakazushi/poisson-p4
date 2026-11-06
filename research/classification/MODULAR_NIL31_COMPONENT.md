# モジュラー Jordan 型 (3,1) から得る新しい25次元成分

2026-11-07。課題 E_multiple_hyperplane。ここでの「階数二」はモジュラー線形部を指し、
Poisson 構造自身の一般階数は4である。旧 positive-weight NEW25 とは別の族。
全分類の完了、primitive β 像の優越性はこのノートの結論に含めない。

## 1. 全90係数から得た源と n7≠0 開の尽くし

affine 座標を `(x,y,z,w)` とし、モジュラー場を `Z=y∂x+z∂y` と固定する。
全定数6、全一次24、全二次60係数を残して

\[
 \operatorname{div}P=Z,\qquad [Z,P]=0,\qquad Q\wedge Z=0
\]

を線形に解いた。全解は15変数 n0,…,n14 であり、最高部の Pfaffian と Jacobi は恒等的に零。
残る全 Jacobi 係数は次の六式に尽くされる（非零定数倍を除く）。

\[
\begin{split}
 &(2n_{10}+n_{14}+1)n_7+2n_{11}n_6,\
 &n_{11}n_7,\
 &n_{11}n_5+n_7n_9,\
 &(6n_{10}+2)n_6+3n_{13}n_7,\
 &(2n_{10}+1)n_5+n_{12}n_7+2n_6n_9,\
 &n_5n_9.                                                   \tag{1}
\end{split}
\]

`n7≠0` では全解が正確に

\[
 n_{11}=n_9=0,\quad n_{14}=-1-2n_{10},\quad
 n_{13}=-\frac{2(3n_{10}+1)n_6}{3n_7},\quad
 n_{12}=-\frac{(2n_{10}+1)n_5}{n_7}.                         \tag{2}
\]

自由変数は `n0,n1,n2,n3,n4,n8,n10,n5,n6,n7` の全十個。
これは既約源 A9×Gm で、標本的な subfamily ではない。

読みやすく `r=n10,p=n5,q=n6,t=n7` とし、
`A(z,w)=n0+n1 z²+n2 zw+n3 w²+n4 z+n8 w` と書くと全六括弧は

\[
\begin{split}
 P^{xy}&=A(z,w)+\tfrac r2y^2-\tfrac t2y,\
 P^{xz}&=z(ry-t/2),\
 P^{xw}&=p+qz+tw-\frac{(2r+1)p}{t}y
          -\frac{2(3r+1)q}{3t}yz-(2r+1)yw,\
 P^{yz}&=rz^2,\
 P^{yw}&=-\frac{(2r+1)p}{t}z-\frac{2(3r+1)q}{3t}z^2-(2r+1)zw,\
 P^{zw}&=0.                                                   \tag{3}
\end{split}
\]

全 Jacobi は恒等的に零であり、

\[
 \operatorname{Pf}P=-\tfrac12z^2\left(p+\tfrac23qz+tw\right). \tag{4}
\]

t≠0 により全点で一般 rank4。
h-Casimir 斉次化の実退化因子は

\[
 h^2z^2\left(ph+\tfrac23qz+tw\right),                        \tag{5}
\]

であり、三つの一次因子は独立である。

## 2. 全周囲接空間による成分性

全自由変数を1とし、二ベクトル全体を6倍する。
この整数代表の六括弧は

\[
\begin{split}
 P^{12}&=6h^2-3hy+6hz+6hw+3y^2+6z^2+6zw+6w^2,\
 P^{13}&=-3hz+6yz,\
 P^{14}&=6h^2-18hy+6hz+6hw-16yz-18yw,\
 P^{23}&=6z^2,\
 P^{24}&=-18hz-16z^2-18zw,\
 P^{34}&=0.
\end{split}                                                  \tag{6}
\]

P の係数が x に依らないことは全周囲の変形を制限する理由にはしない。
全150二次 bivector、全350射影 Jacobi 方程式で線形化 D を作った。
V は全 vertical24方向、M は V と全 PGL 方向・scalar、F は M と全十パラメータの実微分。
正確な `DF=0` と非零有理整数 minor から

\[
 \operatorname{rank}D=100,\quad \operatorname{rank}V=24,\quad
 \operatorname{rank}M=47,\quad \operatorname{rank}F=50.          \tag{7}
\]

producer は非零有限体 minor を標数零の階数下界にだけ使用する。
`DF=0` と F の下界50が D の上界100も与えるため、ここでの等号は有理数上の等号である。
同様に F の階数上界も ker D から50となる。

よって全射影接空間の上界と実族の次元下界は

\[
 150-100-24-1=25=50-24-1.                                    \tag{8}
\]

源の PGL 飽和像閉包を C_N31 とすると、実微分の階数から dim C_N31≥25。
全局所次元は接空間次元25以下なので、この点は滑らかで唯一の全空間成分を通る。
既約な C_N31 はその成分と同次元であるから一致する。
従って **C_N31 は全 Poisson 空間の25次元既約成分であり、(3)の全源を覆う**。
精密な局所 nil 厚みや完成環の追加計算は不要。

## 3. 従前成分との相違

成分性を得た後では、次元が25と異なる既知成分とは異なる。
同次元の MPS 成分 T01,T39,T40 は SNC 点を持つ。
一方、C_N31 の一般 rank4 部分では、射影的掛け算写像の閉像により
全点の実退化因子が `H²L²M` の形を持つ。従って SNC 点を持たず、三つの MPS 成分とは異なる。

旧 positive-weight NEW25 の全 rank4 部分は三乗一次因子を持つ閉条件 `H³Q` に入る。
(5)の独立な三因子の重複度は2,2,1なので、その条件に入らない。
従って旧 NEW25 とも相異なる。CA/CB は28次元であり、同じ因子形を部分的に持つことを
理由に25次元全成分を同一視する必要はない。
一般 rank2 の成分とは閉 rank 条件で区別する。

この相違は従前台帳を基準にする。他の新しい nilpotent 成分は別に記録された22次元なので、
それらとの相違も次元から従う。正式な個数更新は統合担当が行う。

## 4. n7=0 の補集合を残す

n7=0 では (1) と一般 rank4 から n11=0 となり、

\[
 \operatorname{Pf}P=z(n_5+n_6z)(n_{10}z+n_9).
\]

全残条件から二枝だけが残る。

| 枝 | 全固定条件 | Pfaffian |
|---|---|---|
| n5≠0 | `n10=−1/2,n6=n9=n11=n7=0` | `−n5 z²/2` |
| n6≠0 | `n10=−1/3,n5=n9=n11=n7=0` | `−n6 z³/3` |

両枝で n12,n13,n14 は自由であり、(2)の n14 条件を押し付けない。
これらの全支持を本成分へ含めたとは主張しない。
他の nilpotent rank2 affine 型も別義務として保持する。

## 5. rank1 の別の全十変数源もこの成分へ入る

`MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §4.3 の rank1 nilpotent、c=0、b≠0 の全源は
実退化因子 `h²x²w` を持つ。そこで用いるパラメータ A,d0,f,d のまま、
もう一つの重複超平面 x=0 を無限遠に取る。
canonical representative の x-cofactor は

\[
 Y_x=-(Az+d_0h+(f+1/5)x+dw/3)\partial_y+h\partial_z.
\]

x-Casimir gauge の発散は−5Y_x。x=1 上で

\[
 [5Az+5d_0h+5dw/3+5f+1]\partial_y-5h\partial_z.                \tag{9}
\]

A≠0 なら線形部は Jordan 型 (3,1)、定数部はその image に属する。
従って平行移動と線形変換で §1 の型へ戻る。
n7≠0 も必要である：もし n7=0 なら §4 の Pfaffian は、鎖末端の固有 covector z のみの
多項式となる。一方 (9) の `im(N*)²` は h の線であり、実 Pfaffian は定数倍の h²w。
w は h と独立な invariant covector なので矛盾する。

従ってこの rank1 十変数源の A≠0 全開部分は C_N31 に入る。
元の全源は既約で、C_N31 は閉集合なので A=0 の全点も入る。
この橋は源担当が独立に cofactor、Jordan 型、n7 開の必要性まで再導出し確認した。

## 6. 再現資料と検証の射程

* 全五つの nilrank2 affine 型と全90係数：
  `code/classification_multiple_hyperplane_rank2_nil.py`、
  `reports/CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json`。
* 本成分の全源と全150接空間：`code/classification_rank2_nil_component.py`、
  `reports/CLASSIFICATION_RANK2_NIL_31_COMPONENT_20261107.json`。
* 有限契約：`jobs/multiple_hyperplane_rank2_nil.json`、`jobs/rank2_nil_31_component.json`。

独立再構成は `reports/CLASSIFICATION_RANK2_NIL_31_COMPONENT_REVIEW_20261107.json`。
全90線形系 rank75/null15、全15源の係数一致と階数15、六残式の必要性、(2)の三角消去、
全十変数の raw Jacobi/Pfaffian、QQ 上の D100/V24/M47/F50、全61実方向の D 消滅を確認した。
rank1 十変数源からの橋も同じ独立担当が確認した。
本 producer の実測は全五源14.65秒、成分確認0.53秒。
全134 S61式の再生成や primitive β 像の優越性を、この二計算が確認したとはしない。
