# nil22 A=1/4：全源の実Poisson曲線とSNC被覆

2026-11-07。E_multiple_hyperplane。
`MODULAR_NIL22_TOP_SUPPORT.md` の A=1/4 の全 rank四源は MPS の既知成分に入る。
証明は全源の一軌道化、真正な射影 Poisson 曲線、その一つの非零点の SNC 軌道極限を使う。
成分の具体的な T 番号は同定しない。接空間や形式切片の次元はこの包含証明に不要である。

## 1. 全七変数源の正確な一軌道化

u=n10 とし、最初に `W=w+4u,Z=z+4uy` と置く。
残る係数を

\[
 b=n_6,\ c=n_9,\ t=n_{12},\ d=n_2-4un_6,\
 e=n_7-4un_8+16u^2n_9,\ f=n_8-8un_9
\]

とする。さらに

\[
 \delta=4c/3,\ j=-2f,\ g=-4(e-\delta t),\
 e_0=4(b-\delta f)/5,\ c_0=d+16tc^2/9-\delta e
\]

と置き、行列式一の affine 変更

\[
 X=x+\delta Z+e_0W+c_0,\quad Y=y+\delta W,
 \quad T=Z+jW+g
\]

を行う。元の全七変数源の全六括弧は正確に

\[
\begin{aligned}
 P^{XY}&=X-Y^2/2,&P^{XT}&=XW-YT/4,&P^{XW}&=YW/4,\\
 P^{YT}&=-T/4+YW/2,&P^{YW}&=W/4,&P^{TW}&=t+W^2/4
\end{aligned}                                                    \tag{1}
\]

となる。Pf(P)=t(X−Y²/2) だから一般 rank四は t≠0 と同値である。
斉次化すると

\[
 P=N\wedge S+t h^2\partial_T\wedge\partial_W,
 \quad N=Y\partial_X+h\partial_Y+W\partial_T,
 \quad S=-X\partial_X-Y\partial_Y/2-T\partial_T/4+W\partial_W/4.
\]

二つの Jordan block の相対 scale によって t は任意の非零平方倍へ移り、
複素数体上では全 t≠0 が一軌道になる。
公開された全源の凍結 witness は t=7 であり、その値を保って以下を実行した。
全係数の affine 恒等式、全 Jacobi、Pfaffian、同一 witness の輸送は
`CLASSIFICATION_NIL22_QUARTER_CANONICAL_20261107.json` に保存した（0.172秒）。
この変換は別担当の手導出を、本ノート担当が全 coefficient identity として検算した。

## 2. 元点からの真正な曲線

以下では (h,X,Y,T,W) を (h,x,y,z,w) と書く。
4倍した (1) の t=7 を s=0 として、非零係数

\[
\begin{aligned}
 P_s^{hz}&=-2s w(sx+y),&P_s^{xy}&=4hx+2s^2x^2-2y^2,\\
 P_s^{xz}&=4xw-yz,&P_s^{xw}&=yw,\\
 P_s^{yz}&=-hz-s^2xz-2sxw+2yw,&P_s^{yw}&=w(h+s^2x),\\
 P_s^{zw}&=28h^2+3szw+w^2
\end{aligned}                                                    \tag{2}
\]

を持つ二次 lift を取る。(2) は全 s で元の350本の射影 Jacobi 方程式を満たす。
この恒等式は有限な Taylor 解ではなく、全パラメータの多項式恒等式で確認した。
実際の曲線は有限 IFT で発見されたが、その導出手順や切片の完全性をこの被覆の仮定にしない。
曲線そのものと元点の一致だけを用いる。

## 3. s=1 の明示SNC軌道極限

s=1 において `H=h,U=x+y,V=y−x,Z=z,W=w` と線形変更すると、非零括弧は

\[
\begin{aligned}
 P^{HZ}&=-2UW,&P^{UV}&=4H(U-V)-4UV,\\
 P^{UZ}&=2UW-UZ-HZ,&P^{VZ}&=(-2U+4V)W+VZ-HZ,\\
 P^{UW}&=(U+H)W,&P^{VW}&=(H-V)W,\\
 P^{ZW}&=28H^2+3ZW+W^2.
\end{aligned}                                                    \tag{3}
\]

重み `(H,U,V,Z,W)=(1,0,0,−1,1)` を取る。
各非零係数単項式について「係数の重み−微分二本の重み」は非負である。
従って対応する対角 PGL の一径数族（必要ならパラメータを逆に取る）の極限は

\[
 P_0^{UV}=-4UV,\quad P_0^{UZ}=-UZ,\quad P_0^{VZ}=VZ,
 \quad P_0^{UW}=UW,\quad P_0^{VW}=-VW,\quad P_0^{ZW}=3ZW.
\]

その実 Pfaffian 五次式は `−12HUVZW` で、五つの独立な超平面の SNC 因子を持つ。
全変換係数、全非零単項式の重み、極限の全括弧と Pfaffian は
`CLASSIFICATION_NIL22_QUARTER_SNC_20261107.json` に保存した（0.093秒）。

曲線 (2) を含む全 Poisson 空間の既約成分 C を一つ取る。
連結 PGL は C を保つので、C は (3) の全軌道とその SNC 極限を含む。
従って既存 MPS 定理により C は MPS の40成分の一つである。
元点 P0 も曲線 (2) に属し、§1の全源は P0 の一軌道だから、全源が C に入る。
generic-structure 担当が全重み、係数行列の Pfaffian −12、閉成分へのこの論理を独立に手確認した。

## 4. A=−2 の小さい別証拠との共通実装

同じ focused checker の `--minus-two` は、親担当の A=−2 実曲線を別に検証する。
入口は `CLASSIFICATION_NIL22_MINUS_TWO_WEIGHTED_SLICE_20261107.json` の全曲線へ
`tau0=−4s,tau1=3s` を同時代入したもの。
全350恒等式、s=0 の元点一致、全 s≠0 の同一射影軌道を与える重みを確認した。
s=1 の五次式は厳密に Q5·C3 で、Q の行列式は 229582512≠0。
C=wA+g の二次式 A は rank三であり、`A(0,3,1)=0,g(0,3,1)=288`。
よって A は絶対既約で gcd(A,g)=1、Gauss の補題により C は絶対既約である。
既存 `SMOOTH_QUADRIC_CUBIC_SUPPORT.md` の全点定理を適用すれば、非零曲線全体が MPS に入り、
閉性により元点とその全源が入る。ここにも二法線の scheme 全構造の決定は不要である。
この小証拠を `CLASSIFICATION_NIL22_MINUS_TWO_MPS_20261107.json` に保存した（0.491秒）。

共通 checker は `code/classification_nil22_quarter_snc.py`。
minus-two 拡張前の quarter checker の完全 bytes と元 hash は minus-two report 内の
`historical_quarter_producer_bytes/sha256` に保存し、元 quarter report と job は変更していない。
旧入口の正確な再現には、その文字列を一時ファイルへ戻して使える。
現在の checker で quarter を再実行する場合は、既存 output を保護した別 checkout で行う。
有限契約は `jobs/nil22_quarter_snc.json` と `jobs/nil22_minus_two_mps.json`。
