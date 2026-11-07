# Ga9 の二根型から得られる既約非正規五次因子の23次元成分

2026年11月7日。`E_reduced_nonnormal` と Ga9 の exceptional membership のための結果。
以下は有理数体上の全係数計算と本文の幾何的証明による。別作者の Fraction-only レビューも
全項目 PASS し、統合担当も受理した。全 Poisson 空間の尽くしは主張しない。

## 1. 明示的な一般元

座標を `(a,b,c,d,e)` とする。次の二次斉次二ベクトルを考える。
記載しない `{a,b}` は零で、反対称性を補う。

    {a,c}=48a²,                    {a,d}=48ab,
    {a,e}=32ac+b²,                 {b,c}=−72ab,
    {b,d}=−12b²,                   {b,e}=−40ad−8bc,
    {c,d}=72ad+48bc,               {c,e}=−48ae+6bd+32c²,
    {d,e}=72be−88cd.

これは projective Jacobi の全350係数を満たす。符号付き4行 Pfaffian を `pf_i` とし、
実際の反標準因子を `D=Σx_i pf_i` と取ると、直接積展開で `D=60F` となる。ここで

    F=−80a³d²−240a²b²e+160a²bcd+6ab³d−b⁴c.                 (1)

この元の PGL5 軌道閉包を `C` とする。`C` は射影次元23の全周囲成分であり、
一般階数4、一般 PGL5 軌道次元23を持つ。従って一般元は上の一つの明示的軌道で表される。

## 2. 全周囲成分性

計算の出発点は、§4 の実際の cubic curve の非零元 `P1` である。
全150斉次係数を用い、Euler 鉛直24方向を保持した全350行 Jacobi 微分行列 `D1` と、
実 GL5 方向・Euler 鉛直方向を合わせた行列 `G1` に対し、厳密に

    rank D1=102,  rank G1=48,  D1 G1=0

を得た。非零の有理小行列証明書は
`CLASSIFICATION_GA9_DISTINCT_CURVE_ORBIT_20261107.json` に保存される。
射影接空間の次元は `150−24−1−102=23`、実射影軌道次元は `48−24−1=23` である。
軌道が与える局所次元下界と接次元上界が一致するので、この点は滑らかで、
局所既約成分は一つ、軌道閉包そのものが全 Poisson 空間の23次元既約成分になる。
この証明には Ga9 原点の complete IFT、六法線切片、Gröbner 基底を用いない。

上の簡約一般元への変換も有理可逆である。旧座標 `(h,x,y,z,w)` と新座標の間の行列は

    (h,x,y,z,w)^T = B (a,b,c,d,e)^T,
    B=[[  0,−4, 0,−18, 0],
       [  6, 0, 4,  0, 2],
       [863, 0,87,  0,−9],
       [  0, 3, 0,  1, 0],
       [ 99, 0,11,  0, 3]].

旧 P1 の Euler 鉛直成分を除き、`B^(-1) P1(Bu) B^(-T)` と非零定数倍を取ると§1になる。
行列、全9括弧、旧因子と式(1)の正確な比例定数は
`CLASSIFICATION_GA9_DISTINCT_COMPONENT_AUDIT_20261107.json` に保存した。
旧点の一次元射影安定化群の Lie 代数は半単純で、B の列の固有値は
`(−50,−25,0,25,50)` である。これは成分性の別の仮定ではない。

## 3. 絶対既約性と非正規性の直接証明

式(1)を `C[a,b,c,d][e]` の一次式と見る。e の係数は `−240a²b²` である。
残項は `a=0` で `−b⁴c`、`b=0` で `−80a³d²` なので、係数と残項の最大公約因子は1。
従って原始一次多項式であり、Gauss の補題から C 上で既約である。
QQ 上の因数分解だけから絶対既約性を推測していない。

`∂eF=−240a²b²` より特異点は a=0 又は b=0 上にある。
a=0 では `∂cF=−b⁴` より b=0、b=0 かつ a≠0 では
`∂dF=−160a³d` より d=0 となる。逆包含は全微分への代入で従うので

    Sing(F)_red = {a=b=0} ∪ {b=d=0},

二つの射影平面であり、その交わりは射影直線 `{a=b=d=0}` である。
既約三次元超曲面の特異集合が二次元なので R1 条件が破れ、因子は非正規である。
特異スキームの正確な非被約構造はこの結論に必要なく主張しない。

最初の平面に沿う F の一般重複度は3、二番目では2である。
§1 の projective 二ベクトルは最初の平面上で一般階数2、二番目上で零となる。
後者の零とは斉次代表そのものの零ではなく Euler 鉛直を除いた零である。
この違いを保ち、全ての特異面を一様な ordinary double 型とは扱わない。

## 4. 全 Ga9 二根型源の包含

対象は cofactor `Z=y∂x+∂z`、`c(y)` が相異なる二根を持つ型である。
`NIL211_SHORT_TRANSLATION_REVIEW.md` の全源分類は、可逆正規化によって
canonical 五変数 `(D1,D2,D3,D4,d0)` とその normalizer から全支持を覆う。
元の全28係数、全 affine Jacobi/Pfaffian、補集合はその証明で保持される。
その PGL 飽和射影像閉包 `S` は既約であり、同ノートの全源上界

    dim S ≤ 10+9+4−1=22

が成り立つ。五変数零の canonical 点 P0 は tensor scale3 で

    {x,y}=2h²−3y²,       {x,z}=−3hx+6yz,
    {x,w}=−3hz−3wy+3xy, {y,z}=3hy,
    {y,w}=3y²−3h²,      {z,w}=−3hw+9hx−9yz,
    {h,−}=0

となる。この点の全実軌道行列は rank47 なので射影軌道は22次元で、上界と合わせて
`S=closure(PGL5·[P0])` である。`GA9_DISTINCT_WEIGHT_BLOCKS` の全五変数微分も保持済みである。

ここで `CLASSIFICATION_GA9_DISTINCT_CURVE_20261107.json` は、P0 から P1 への
literal な ε 三次の二ベクトル `Pε` の全10括弧を保存する。
全350 projective Jacobi が `QQ[ε,h,x,y,z,w]` で恒等的に零であり、ε=0 は上の P0、
ε=1 は§2の P1 である。さらに

    g(t)=diag(1,t,1,t,t),  t^(-1) g(t)_*Pε=P_(ε/t)

を全係数で満たす。従って全 ε≠0 は同じ実 PGL/scalar 軌道で、P0 は C に属する。
C は PGL 安定な閉集合なので `S⊂C` となり、canonical 点だけでなく、
正規化を戻した全二根型 Ga9 源とその全パラメータ境界が被覆される。

この計算の発見に使用した complete 六法線 IFT/support は公開のまま保持するが、
成分性・全源包含の最小証明は実曲線、実軌道、全源上界だけで閉じている。
他の三つの `c=0` 源への包含は主張しない。

## 5. 凍結した既存98成分との相違

比較台帳は `CLASSIFICATION_GA9_DISTINCT_PRIOR98_20261107.json` で凍結した。
全98項目と source hash を保持する。同次元以外の成分とは、両者が全周囲成分であることと
次元の違いから相異なる。同次元の既存成分は MPS29個と非MPS27個である。

MPS29個については SNC 点の存在から推測せず、既に全源を覆うと証明された
一般 period と全 smoothing parameter の**完成された族そのもの**の実五次因子を計算した。
`CLASSIFICATION_MPS23_DIVISOR_FACTORS_20261107.json` は次を全変数の恒等式で記録する。

| IDs | 全一般源の座標因子次数 |
|---|---|
| T04–T15 | 1+1+3 |
| T18–T20 | 1+2+2 |
| T22–T31 | 1+4 |
| T32–T35 | 2+3 |

全29件で、Euler 収縮した実 Pfaffian 因子と保存因子の積を `QQ(p,t)[x0,…,x4]` 上で直接比較した。
全源閉包が対応する MPS 成分であることは
`research/notes/MPS_S61_CATALOGUE_20261103_v62.md` §1–2 と、その MPS 局所 Torelli 適用による。
従ってこれは一点の因数分解や smoothing graph の直観ではない。

非MPS27個は、各既証明の密な源で 1+4、5H、1+1+1+2、2+3、1+1+3、1+2+2、2H+F3 のいずれかを持つ。
全項目の元の証明参照と因子次数・重複度は `GA9_DISTINCT_COMPONENT_AUDIT` に残した。
可約五次形式の集合は、各次数分割の射影乗法写像の像の有限和であり閉である。
その affine cone は零形式も含め閉なので、Poisson の Pfaffian 写像で引き戻せば階数落ちの境界も含む。
従って全同次元既存成分はこの閉可約条件に入り、式(1)の絶対既約因子を持つ C とは一致しない。
これで凍結した98個のいずれとも異なる。

## 6. 独立再現と状態

公開の最小入力は、実曲線 report、P1 全周囲/orbit report、五源の有限分類と二根型上界、
凍結 prior98、MPS29全源因子 report である。対応コードと job はそれぞれ

* `classification_ga9_distinct_curve.py` / `jobs/ga9_distinct_curve.json`、
* `classification_ga9_distinct_curve_orbit.py` / `jobs/ga9_distinct_curve_orbit.json`、
* `classification_ga9_distinct_component_audit.py` / `jobs/ga9_distinct_component_audit.json`、
* `classification_mps23_divisor_factors.py` / `jobs/mps23_divisor_factors.json`。

job gate は契約と入力 hash を検査するだけで数学の証明ではない。
本 producer の正確計算は完了し、式(1)の絶対既約・非正規性は上の手証明で閉じている。
別作者 `classification_ga9_distinct_component_review.py` は complete IFT/GB、CAS、producer import に
依存せず、literal curve、全350、D102/G48、P0G47、可逆座標変換と全9括弧・式(1)、
全29一般 MPS source の因子積、全98比較を Fraction 算術だけで再構成し、全項目 PASS した。
証拠は `CLASSIFICATION_GA9_DISTINCT_COMPONENT_REVIEW_20261107.json`（1.176秒、19,328KiB）と
`GA9_DISTINCT_QUINTIC_COMPONENT_REVIEW.md`。統合担当が受理済みであり、台帳の更新は同担当に一元化する。
