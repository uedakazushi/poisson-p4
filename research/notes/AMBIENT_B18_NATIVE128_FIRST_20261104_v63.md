# actual native128 syzygy の全一次係数

2026-11-04、v63。実 \(\operatorname{Res}_4(N^tS)\) の48×128成分を、
\(\bar R/\mathfrak m^2\) 上で正確に構成した。
八変数は actual CA F9 の \((z_{A,0},\ldots,z_{A,6},r_A)\) であり、
\(\ell=5(s_A-r_A)=0\) を正確に使う。
出力3862個の非零QQ係数、実行2.098秒、約49 MB。
新しい physical source jet、105単位方程式、点やrank searchは用いない。

## 入力と有限代数

既存 Artin105 の正確な解から両48 source の全一次係数を読み取る。
元 scalar coefficient DAG を dual9 に評価し、両 physical150 trace-pair の一次値を得た。
この段階は元の正則係数関数への代入であり、source IFT を再実行しない。
全 physical odd 係数の定数は零、even 定数は既存 (P_0\) と一致することを確認した。

同じ source の flag matrix は
\(H_i=G_{\gamma,i}\operatorname{core\_group}_i(\tau_i+s_i/2)\)。
γ0 は tensor scalar なので旗行列へ含めない。
両 frozen core group の (p^2\) 項を含め、trace shift と monic
\(\tau_i^2=\nu_{i,1}\) を保つ。

本体の代数は、square-zero な八つの  \(\epsilon\) を持つ dual9 と
\((1,\tau_A,\tau_B,\tau_A\tau_B)\) の積で長さ36である。
その maximal ideal は五乗零なので全単位逆元を有限幾何級数で正確に計算する。
\(\epsilon\tau_A\tau_B\) や、\(\tau_A^3\tau_B\) が monic を通して一次baseへ落ちる寄与を捨てない。

## 同じ group30 の一次補正

saved nil center の group15 Jacobian を
\(A_0+\eta_AA_1\) とする。\(A_0\) は既存15×15行列と全成分一致した。
225個の mixed \(\eta_A\delta\) 係数を長さ4代数で計算し、(A_1\) を零とは仮定しない。

\[
 J_{30}=\begin{pmatrix}A_0&0\\A_1&A_0\end{pmatrix},\quad
 \det J_{30}=4/9
\]

の逆行列を固定QQ算術で確認した。
次に長さ18の CA cover で240 RHS係数を評価し、一回の固定Jacobian補正で
30×8個の未知一次係数を求めた。全30 component constraint をその代数で正確に零と確認した。
補助 \(\eta_A\) と physical \(\tau_A\) を混同せず、source と同じ偶奇順序を保持する。

## 交叉外の正確な式

実 normalizer pair は自由四係数環で扱い、まだ (F_{12}\) を課さない。
一次では (F_{12}\) が非零なので、原点だけの「actual pivot=rational root」や
「(dG_8=0\)」という shortcut を用いない。
最初の五式も完全な monomial telescoping で構成した。
次の全式を長さ36で確認した。

\[
 F_{17}(\mathrm{actual14})=R_{17}F_{12},\qquad
 dH_{12}=R_{12}F_{12},\qquad
 TH_{12}=(G_8,0),\qquad NF_{12}=d(G_8,0).
\]

(d\) は反転しない。(R_{12},T,N\) の原点四成分は既存 X0 と全て一致する。
\(G_8\) の原点も、零でなく
\((1/2,25/8,0,0,0,5,-125/2,0)\tau_A\tau_B\) と一致した。

native (S\) の source 符号をそのまま使い、(N^tS\) を形成する。
四つの基底乗数による Res4 を、変化する monic \(\nu_{A,1},\nu_{B,1}\) を保ったまま抽出した。
全6144原点係数が、別に凍結済みの actual X0 の48×128行列と一致した。

## 出力と範囲

主要データ `actual-native128-X-mod-m2.json` のキーは
`source_columns128_rows48_sparse_dual9`。
128列、各 entry は `[row,[[dualIndex,QQ],...]]` である。
row0..47、dualIndex0は定数、1..7はzA0..6、8はrA。
fullX の下八行は零である。data hash を report へ直接固定した。

五つの phase file は source firsts、固定30Jac、240解、actual14/R17/R12、T/N/G8を保存する。
全出力35,825QQ slots、225 KBであり、保存時の累積出力＋指定livebufferの150,000QQ guard と
32 MB guard を通過した。
この packet はµT、Tの長さ、Kのtorsionや原点の深さを決めない。
次の小行列組立は、別担当が既受理の first640 一次データとこの128列を結合して行う。

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_native128_first_20261104_v63.py
```
