# native128 syzygy の一次係数：既存 source からの有限 nil 代数 preflight

2026-11-04。対象は actual \(\operatorname{Res}_4(N^tS)\) の128列の
\(\bar R/\mathfrak m^2\) 上の値である。新しい physical source jet や105単位方程式は解かない。
この packet は source・指数・support の検査だけであり、128列の一次係数はまだ計算しない。

## 存在する入力

\(\bar R\) の座標は実 CA F9 の \((z_{A,0},\ldots,z_{A,6},r_A)\)、
\(\ell=5(s_A-r_A)=0\) である。既存 Artin105 の正確な値から、
全105未知数のこの八座標に関する一次係数を読み取れる。
原点と一次 monomial は quotient の高次 alias と混同しない。
両 normalizer48の全 parameter pair、\(s,t\)、physical 共通九底、両 \(\nu\) の一次値を保存した。

元の scalar coefficient DAG を dual9 で評価すれば、全 physical150出力の
normalizer trace-pair \((a_i,b_i)\) の一次値が得られる。
必要閉包は17,926 nodes、最大liveは400ベクトル、係数上界3600である。
source48解の補正は行わず、既存105一次値を代入するだけである。
fold由来の \(p^2,p^3\) の寄与は、元の係数DAGの \(s,t\) に既に正しく含まれる。

旗 frame は元の exact 式
\(H_i=G_{\gamma,i}\operatorname{core\_group}_i(\tau_i+s_i/2)\) であり、
両 core group は保存済み二次行列である。GL18 の位置順序も元の source と同じものを固定した。

## 新たに必要な小さい gauge 微分

元の group30 は group15 を形式的補助変数 \(\eta_A\) の偶奇へ展開したもの。
その中心は even15が零、odd15が
\((-1,0,0,2,0,0,-2,0,-2,0,0,-1,0,-4,0)\) である。
これは physical \(\tau_A\) を任意の独立変数と取り違えた値ではない。

group15 constraint の Jacobian を saved nil center で
\(A_0+\eta_A A_1\) とすると、原点group30 Jacobian は

\[
 J_{30}=\begin{pmatrix}A_0&0\\A_1&A_0\end{pmatrix},\qquad\det J_{30}=4/9.
\]

(A_0\) は保存済み15×15 matrixである。
新しく必要なのは、同じ十五constraint・十五group方向の
\(\eta_A\delta\) 係数225個だけである。
\(\mathbb Q[\eta_A,\delta]/(\eta_A^2,\delta^2)\) の長さ4で、同じ literal flag/fiber callback を微分すればよい。
\(A_1=0\) と仮定しない。

RHS は dual9 に対する十五constraintの偶奇成分30×8＝240係数である。
\(\eta_A^2=\nu_{A,1}\) とする長さ18代数で saved gauge center を評価し、
固定 \(J_{30}^{-1}\) を一回適用する。補正の積は base の二乗に入るため、
この一回で全30 component constraint を厳密に検査できる。

最後の native係数は

\[
 C_{36}=\frac{\mathbb Q[\epsilon_1,\ldots,\epsilon_8,\tau_A,\tau_B]}
 {(\epsilon_i\epsilon_j,\tau_A^2-\nu_{A,1},\tau_B^2-\nu_{B,1})}
\]

で計算する。基底は \((1,\epsilon_1,\ldots,\epsilon_8)\times
(1,\tau_A,\tau_B,\tau_A\tau_B)\) で、maximal ideal の五乗は零。
元の \(\tau_A^3\tau_B\)、\(\tau_A\tau_B^3\) が一次 base 係数へ落ちる寄与を消さない。
単に両 \(\tau^2=0\) を固定したまま係数だけ微分してはいけない。

## 原点だけの shortcut を使わない

原点の長さ4代数では physical (F_{12}=0\) であったが、(C_{36}\) では一般に非零である。
必要な正しい検査は

\[
 F_{17}(\mathrm{actual14})=R_{17}F_{12},\quad
 dH_{12}=R_{12}F_{12},\quad TH_{12}=(G_8,0),\quad
 NF_{12}=d(G_8,0).
\]

actual first5 pivot と rational roots は一次では異なるので、最初の五式も
**完全な monomial telescoping**で行う。原点でだけ有効だった derivative shortcut は使わない。
(d\) は反転しない。原点の (G_8\) も零ではなく、
\((1/2,25/8,0,0,0,5,-125/2,0)\tau_A\tau_B\) を正確に再現する。
全144 native (N\) entry の四成分と (G_8\) が凍結 X0 と一致することが必要である。

## 有限 support と出力 gate

native変数は全て maximal ideal にあるので、直接多項式評価には次数4以下、
telescoping係数には元の次数5以下が必要十分である。
実入力の必要項数は、17式186、H12が365、first roots32、second roots53、
membership cofactors355で、最大一多項式46項である。
全 source係数は既存ファイルに存在し、新たな高次 physical係数は不要。

新しいQQ入力は225混合Jac係数、240 RHS、240解係数。
最大の出力128×48×9は55,296 scalar slotsで、diagnosticを含めた保守的上限を120,000係数、
32 MBとする。一回の有限AD・固定Jac微分 job を60秒・2 GBで実行し、各段階を保存する。
cap・異常時に自動延長しない。現在この derivative job 自体は未実行である。

次担当へ渡す schema は `source_columns128_rows48_sparse_dual9`、
各列は `[row,[[dualIndex,QQ],...]]`。
dualIndex0は定数、1..7はzA0..6、8はrA。full X の下八行は零である。
data hashをreportへ直接固定し、原点全係数のX0一致後にのみμTの組立へ渡す。
