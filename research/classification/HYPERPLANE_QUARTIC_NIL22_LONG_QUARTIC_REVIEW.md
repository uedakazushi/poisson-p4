# 型7 nil22 長 translation：全 quartic 七源の独立監査

2026-11-07。PRIMARY E の有限入口を、producer を import しない別の
`Fraction`/疎 Laurent 多項式実装で監査した。**PASS**。
成分性、局所環、既存成分への帰属はこの監査の対象でない。

現行 V2 は公開 producer packet・公開証明ノート・元の nilpotent 源表だけに結合される。
初版 V1 と embedded reviewer は変更せず保存した。移行前に四つの旧 structured math 節が
公開版の対応節と全て一致することも直接比較した。

- 実装：`code/classification_nil22_long_quartic_review.py`
- 有限 job：`reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_REVIEW_JOB_20261107.json`
- 現行証明書：`reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_REVIEW_V2_20261107.json`
- 公開数学入力：`reports/CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_SOURCES_20261107.json`
  の `potential_geometry`, `actual_core_action`, `rank_resonance` および全七 fiber 節と、
  `research/classification/HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_SOURCES.md`。
  全 SHA256 は job と証明書に保持する。再実行には work 内のファイルを使わない。
- 元の全源：`CLASSIFICATION_HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_20261107.json`
  の零始まり row6、20 自由係数。元の全90線形係数、全六括弧、全 Jacobi、Pfaffian
  を既存の独立 Fraction 構築から再生した。

## 1. potential と共鳴の全場合分け

固定 cofactor は (Z=y\partial_x+\partial_y+w\partial_z)。

\[
 (t,u,v,w)=(y,x-y^2/2,z-wy,w)
\]

は Jacobian −1 の多項式自己同型である。これは因子分解の補助であり PGL 同値ではない。
全六括弧を直接 push すると

\[
 P=\partial_t\wedge V+\operatorname{Jac}(\Phi),\qquad
 \operatorname{div}V=-1,
 \quad J(P)=(Q_w,-Q_v,Q_u,0),\quad Q=V\Phi+\Phi.
\]

元の実 Pfaffian は厳密に \(\Phi-Q\)。従って Poisson 支持上は
\(Q=k\) 定数、\(F=\Phi-k\)、\(VF=-F\) である。
実装は20係数を残したままこれらを展開して一致させた。

\(q=n_{14}\)、\(H=\operatorname{Hess}_{u,v}\Phi\)、
\(M=\partial_{u,v}(V_u,V_v)\) とすると

\[
 M^TH+HM+H=0,\quad \operatorname{tr}M=-q-1,\quad
 \operatorname{tr}(\operatorname{adj}H\,(M^TH+HM+H))=-2q\det H.
\]

以下の被覆は複素数体上の全点についての手証明であり、計算による標本推測ではない。

**rank (H=2)。** q=0。定数線形変換で二次項を UV とすると conformal 条件により
M は diag(α,−1−α) となる。\(F=UV+p(w)U+q_1(w)V+r(w)\)、\(Vw=s\) と置く。
\(VF=-F\) の係数比較は \(g=r-pq_1\) に \(sg'+g=0\) を与える。
多項式解は g=0 のみなので F は積に分かれる。この枝全体を排除できる。

**rank (H=1)。** 定数線形変換と active 変数の平方完成で
\(F=U^2+f(w)V+g(w)\)、deg f≤2、deg g≤4。
active 固有値は −1/2、passive 固有値は −q−1/2。
f≠0 では \((n_{13}+qw)f'+(1/2-q)f=0\) なので、最高次数 d により
\(q(d-1)+1/2=0\)。d=1 は不可能、d=0,2 は q=1/2,−1/2。
q=0 でも非零多項式解はない。
f=0 では U 係数が active drift を消し、\((n_{13}+qw)g'+g=0\)。
q=0 は g=0。q≠0 は流れの根からの一冪で、零または偶数冪なら平方差として可約。
既約の場合は d=1,3 のみで q=−1,−1/3。

**rank (H=0)。** \(F=\ell(w)\cdot(u,v)+r(w)\)、deg ℓ≤2、deg r≤4。
q≠0 の流れの根 w0 で ℓ(w0)=0 なら定数項式から r(w0)=0 となり、w−w0 が因子。
従って既約なら ℓ0≠0、M は固有値 −1 を持つ。
根を原点にして \(\ell=\sum\ell_d w^d\) とすると
\((M^T+(1+dq)I)\ell_d=0\)。d>0 の項があれば固有値は −1 と異なり、
trace 式から \(q(d-1)=-1\)。d≤2 より q=−1 のみ。
ℓ が定数だけなら補助座標から元座標へ戻しても ℓ 部分は次数≤2である。
quartic 条件は r の四次項を要求し、4q+1=0、すなわち q=−1/4。
q=0、n13≠0 では M の固有値は −1,0、ℓ の多項式微分方程式は ℓ を定数にし、
r も次数≤2。q=n13=0 では ℓ=f(w)e と r=−f(w)e·b(w) で F が積となる。
この二つの q=0 枝も落としていない。
実 core に H=0 と det(M+I)=0 を入れると q=−1/4 は単位理想、q=−1 では
\(n_{11}=0,n_{19}=-1/2,n_{12}=n_8,n_2=4n_8^2\)。

## 2. 実 PGL 正規化と有限 core

\[
 x'=x+az,\quad y'=y+aw,\quad z'=rz,\quad w'=rw,\quad r\ne0
\]

は固定 Z の actual affine centralizer。det=r² と全20 parameter map、元の全六括弧の
push を直接再生した。q は不変である。

rank1、n17≠0 は a=n18/(2n17), r²=1/n17 により H=diag(−2,0) へ送る。
複素平方根は全非零点を覆い、追加の連続 parameter を生まない。
四つの q の全 core 理想は q=−1/2 だけを残し、n2=−3/16,n8=±1/4。
r=−1 が符号を交換する。

rank1、n17=0 は n18=0,H22≠0。q=1/2,−1,−1/3 は一つの a で
n2=n8=n12=0 へ送る。q=−1/2 の全 core は

\[
 (3n_{19}+1)n_8=0,\qquad
 4n_{12}n_8+(4n_{19}+1)n_2=0.
\]

最高四次は \(w^3(n_8y-n_2w/4)\)。n8=n2=0 の連続 core は全部四次なし。
n8≠0 は n19=−1/3 と n2=12n12n8、a=−3n12,r=n8 で一 core へ送る。
n8=0,n2≠0 は n19=−1/4、n2=ρ² として a=−4n12,r=ρ で一 core へ送る。
rank0 の残りは a=−2n8,r=1 で一 core へ送る。
これらの零枝も含めた七 core は producer 表と一致した。

全 q≠0 なので、さらに actual affine 変換 z'=z+sy,w'=w+sh、s=n13/q で n13=0。
全 core が不変であり、この20 parameter map も全六括弧で検査した。

core の理想一致は producer の GB を実行して再確認したのではない。
元の全 Jacobi から取り出した八変数 core 多項式を使い、各候補について**両方向の理想所属**を
次数≤2の単項式倍による一つの有限 Macaulay 行列で作った。
証明書の各 target は `[generator index, multiplier, rational coefficient]` の和であり、
その和を直接積み戻して一致させた。二次 multiplier は n19 の一次消去後の n19*n8² を
除くために必要で、無制限 GB は使わない。単位理想の枝では 1 の明示所属を保持する。

## 3. 七つの全 fiber と全補集合

各固定 core と n13=0 の元の**全 Jacobi 理想**について、報告された graph 理想と同じ方法で
両包含を証明した。従って単なる源の代入や一方向の解族の検査ではない。

|順|全既約 quartic parameter 空間|逆 core 正規化の parameter 数|中心移動を戻した固定源の上界|
|---|---|---:|---:|
|U、q=−1/2|A4|2|7|
|V、q=1/2|D(n15)⊂A5|1|7|
|V、q=−1|D(n3−4n6n7)⊂A5|1|7|
|V、q=−1/3|D(n5−n1)⊂A5|1|7|
|V、yw³型|A4|2|7|
|V、w4型|A4|2|7|
|rank0、q=−1|D(n6)⊂A5|1|7|

最初と五・六番目は一変数について一次、係数は非零定数倍の w²。
w=0 の余項にはそれぞれ x² または z² の非零定数係数が残り、Gauss により全 parameter で既約。
二・七番目は単位 parameter を係数とする一次式。
三・四番目は平方を除いた余項が非零定数倍の w、w³ なので、C(w) 上で判別式が平方にならない。
各単位の零補集合は平方または w 因子である。
**rank0 の生の compatibility n6*n10=0 は残した。** n6=0 で n10 を勝手に零とせず、
生の全 Pfaffian に n6=0 を入れて w 因子を確認してから n6≠0 上で n10=0 とする。
全20代入、全六括弧、全 Jacobi、実 Pfaffian、四次最高部、canonical C* weights を検算した。

各 fiber は上記の既約開集合。逆正規化と中心移動を戻すと、各固定 cofactor 源は次元≤7。
固定 Z の全 affine orbit 接行列20×20も別途組み、正確な有理数消去で rank12 を再生した。
従って H の選択4、全体 scalar の projectivization1 を含め

\[
 \dim\overline{\operatorname{PGL}_5\cdot\{
 \text{型7の }H+\text{既約四次}
 \}}\ \le 7+12+4-1=22.
\]

よって未知の次元25–27の一般成分をこの型が支配することはない。
成分帰属や anchor の orbit 次元22との同一性は別の証拠を要し、この結論に含めない。

## 4. 実行記録と範囲

有限 gate PASS 後の公開 V2 最終実行：0.7204秒、最大 RSS 16,096KiB、全項目 PASS。
標準ライブラリと既存の独立 sparse/Fraction helper のみを用いた。
初回の reviewer 固有の `shift`/`s` および `^` 表記読替を修正し、
一次 Macaulay 倍では不足する n19*n8² の消去を上記の有限二次行列へ置き換えた。
producer の入力・主張・証拠 bytes は変更していない。
