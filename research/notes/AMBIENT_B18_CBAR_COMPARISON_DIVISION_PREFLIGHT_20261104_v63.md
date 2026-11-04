# Cbar 比較の正則除法と unit 列消去：source preflight

2026-11-04、v63。実 Cbar の三 target は \(\tau_A,\tau_B,\tau_A\tau_B\) である。
実最小分解の (d_1:\bar R^{24}\to\bar R^3\) は、三つの regular ideal

\[
 (z_7,r^3),\quad(z_7+c r^2,r^3),\quad(z_7,r^2)
\]

の八生成元を各 target へ送る。ここで \(\bar R\) は実 CA F9 座標の環である。

actual 比較行列 (V\) の各成分 (f_s\) はその sector ideal に入る。
中央 sector だけ (u_i=z_i+c_i r^2\)、他は (u_i=z_i\) と置き、七つの coordinate-prefix
Hadamard と残りの \(f_s(0,r)/r^{n_s}\)、\(n_s=3,3,2\) を取れば、
具体的な正則係数行列 (m\) が定まり、\(d_1m=V\) となる。
最後の monomial 除法は actual augmentation が零であることから正確に割り切れる。
元の (r/g\) 関数を特殊化した後に割り直す操作ではない。

この (m(0)\) の七座標係数は元関数の一次 (z_i\) 係数である。
最後の係数は (z=0\) または (z=-c r^2\) 上の (r^{n_s}\) 係数である。
実 (H_{20}\) から \(\mathbb Q[r]/r^5\) への両 map があり、保存済み有限値から求められる。
full physical (V:3\times56\) は (Q48\) の三非定数行と (-U8\) であり、
(Q48\) の値は既存 Artin105 の \(\delta,b_A,b_B,\nu_A,\nu_B\) から構成できる。
(U8\) は正確な七座標と actual (g_A\) であり、(g_A\) の三次係数 (125/24\) だけが
最後の (m(0)\) slot に寄与する。未計算の四次 tail を零と仮定する必要はない。

\(V\) が \(\operatorname{im}d_1\) を生成し、(d_1\) は最小24生成なので、(m(0)\) は全射である。
24列 (S\) を選び (M=m[:,S]\) の定数行列式が非零と確認すれば、

\[
 n=E_S M^{-1},\quad mn=I_{24},\quad Vn=d_1
\]

となり、以前の homotopy (h\) は正確に零と選べる。
列を並べて (m=[M,N]\) とすれば

\[
 W=\begin{pmatrix}M^{-1}&-M^{-1}N\\0&I\end{pmatrix},\qquad VW=[d_1,0].
\]

proper40列の source が実際に束縛された後は \(\ker V\cong\operatorname{im}d_2\oplus\bar R^{16}\)。
従って完全な (T\) 表示は、168本の (d_3\) と752本の actual cycle lift により
100生成・920関係まで小さくできる。full56列を使う段階では116生成・936関係であり、
さらに16本の既知 free-zero source unit を actual Res4 adapter で消す。
この preflight 自体はその行列消去や cycle lift を実行しない。

残る cycle (mX_j\) の持上げには、三 regular sequence の具体的 Koszul contraction を使える。
関係 \(\sum g_i a_i=0\) を七 coordinate で割った商で見ると
\(r^n a_7(0,r)=0\) なので (a_7(0,r)=0\) であり、その coordinate Hadamard で最後の成分を除く。
残りは七 coordinate の通常 Koszul contraction で解ける。
これは実正則関数に対する具体的 operator であり、未定の matrix membership 判定への置換ではない。

実 (m(0)\) と unit minor の一回の評価は、後続の
`AMBIENT_B18_FULLV56_COMPARISON_UNIT_20261104_v63` で別に行う。
proper40 adapter、全 cycle 係数、348環での scalar 評価、(K\) の torsion はこの preflight の結論ではない。
