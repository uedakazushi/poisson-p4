# 元 D0 の nilideal 上の実乗法と44出力 cube operator

2026-11-04。既受理の完全比較複体

\[
 R_0^{780}\xrightarrow{X}R_0^{44}\xrightarrow{V}R_0^3,
 \qquad H_1\simeq\operatorname{Nil}(D_0)
\]

に、元の環乗法を literal に付ける。新しい nil 生成元を求めたり、全 cycle を列挙したり、
cube の零性を判定したりするものではない。既存 source 係数の展開、780列の再構築、
新しい jet・点・CAS 計算は行わない。

## 1. 同じ source 上の代表元

actual (C_0=A_{A,0}=R_0+I_A\tau_A\) の九生成元は

\[
 1,\ \beta_1\tau_A,\ldots,\beta_8\tau_A,\qquad
 (\beta_1,\ldots,\beta_8)=(z_{A,1},\ldots,z_{A,7},g_A),\quad\tau_A^2=\nu_A.
\]

source 内の正確な scalar node は
\((7148,7199,7454,7505,7688,7943,8045,505938)\) と \(\nu_A=20086\)。
旧 (g_A\) の関数を近似や新たな非単元除算で置き換えない。
九生成の free module は28 Koszul 関係で (C_0\) へ写る。free module 自体に
associative な環構造があるとは仮定しない。

正規化側の coefficient algebra (C\) は同じ (R_0\) 上で free4、basis は
\((1,\tau_A,\tau_B,\tau_A\tau_B)\)。reduced target は \(E_{\rm red}=C/(dR_9)\)。
\(Q:R_0^{36}\to C\) は9つの (dR_9\) と四 basis 倍であり、その36列を全部保存している。
元の source map の九列は文字通り

\[
 U_0=(1,0,0,0),\qquad U_i=(0,\beta_i,0,0).
\]

これは型だけでなく受理済み実 C0/normalizer source の同一関数である。

cycle を \(z=(a,c)\)、\(a\in R_0^{36},c\in R_0^8\) と書く。単位列消去後の
cycle 条件は \(Vz=0\)。消した第零成分を戻すと

\[
 b_0=Q_0a,\qquad B=\sum_i\beta_i c_i,\qquad Qa=b_0+B\tau_A.
\]

従って (H_1\to\operatorname{Nil}(D_0)\) は (z\) を実 C0 の元
\(b_0+\sum_i c_i\beta_i\tau_A\) の D0-class へ送る。

## 2. literal bilinear operator

二つの cycle (z=(a,c),z'=(a',c')\) に対し、第二の (b_0',B'\) も同様に定める。
9つの四成分 block ごとに

\[
 a''_j=\mu_C(b_0+B\tau_A)a'_j,\qquad
 c_i''=b_0c_i'+b_0'c_i
\]

と置く。ここで

\[
 \mu_C(b_0+B\tau_A)=
 \begin{pmatrix}
 b_0&\nu_AB&0&0\\ B&b_0&0&0\\
 0&0&b_0&\nu_AB\\0&0&B&b_0
 \end{pmatrix}.
\]

44出力は全て同名 certificate の finite DAG に保存した。必要な scalar input は
36個の実 (Q_0\) 係数、八個の \(\beta_i\)、実 \(\nu_A\) だけである。

native source の四 basis 倍という Q の構造により、全9 block で

\[
 Qa''=(b_0+B\tau_A)Qa'
\]

が正確に成立する。そのため cycle 上では

\[
 Q_0a''=b_0b_0'+\nu_ABB',\qquad
 \sum_i\beta_i c_i''=b_0B'+b_0'B.
\]

これは元 C0 の実積そのものである。したがって product は cycle となる。
さらに cycle 条件を課す前にも、残差の順序を \((\tau_A,\tau_B,\tau_A\tau_B)\) として

\[
 Vz''=
 \begin{pmatrix}b_0&0&0\\0&b_0&\nu_AB\\0&B&b_0\end{pmatrix}Vz',
\]

\[
 Q_0a''-(b_0b_0'+\nu_ABB')=\nu_AB\,(Vz')_{\tau_A}
\]

という普遍恒等式がある。従って even 成分の公式を非 cycle の任意入力へ
無条件に拡張したとはしていない。

## 3. 境界と quotient 上の乗法

受理済み完全比較定理により、cycle の D0-class が零であることと、その cycle が
\(\operatorname{im}X\) に入ることは同値である。境界の C0 代表元は D0 で零であり、
その積も零である。従って上の operator は H1 に降り、正確に実 nilideal の乗法を与える。
780個の境界係数を新しく解いたという主張は不要であり、行っていない。
chain-level の (a''\) が可換・結合的でなくても、D0 の積が同じなので、差は完全な
X-boundary である。28 Koszul を捨てた議論ではない。

## 4. 一入力 cubic operator と正確な次の gate

一つの cycle (z\) に対し

\[
 A=b_0^2+\nu_AB^2,\quad T=2b_0B,\qquad
 a^{[3]}_j=\mu_C(A+T\tau_A)a_j,
\]

\[
 c_i^{[3]}=(3b_0^2+\nu_AB^2)c_i
\]

と置く。この44出力も別の literal DAG root として保存した。cycle 上で回復する even 成分は
\(b_0^3+3\nu_Ab_0B^2\)、odd scalar は \(3b_0^2B+\nu_AB^3\) であり、元 D0 の三乗である。

char0 の三次 polarization により

\[
 \operatorname{Nil}(D_0)^3=0
 \quad\Longleftrightarrow\quad
 z^{[3]}\in\operatorname{im}X\text{ for EVERY exact cycle }z\in\ker V.
\]

未計算の全 kernel を free44 の任意 basis と同一視しない。一つの実 cycle に対する
三乗の非 membership は非零性の証拠になるが、有限 truncation の membership は零性の
証明にはならない。全域の零性には、完全な nil-cycle 生成元またはこの cubic map の
全像に対する正当な関係計算が必要である。

既受理の theorem はこの cubic image が原点有限長 (H^0(D_0)\) に入ることを与える。
しかし全 H0、旧 Tor63 image、cubic image を同一視しない。ここでは未決 gate を
未指定の積から実 source-bound finite operator へ置き換えた。
