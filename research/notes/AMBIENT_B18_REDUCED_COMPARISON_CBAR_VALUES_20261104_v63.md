# actual reduced-comparison Cbar：既存 B11 値からの正確な cokernel

2026-11-04。元の比較 cokernel (C=\operatorname{coker}(D_0\to E_4)\) と
\(C'=\operatorname{coker}(D_0\to E_{\rm red})\) を区別する。
既知の \(\bar C=C/\ell C\) の三つの literal cyclic sectors に、二つの normalized nil
generator の実像を加えて \(\bar C'=C'/\ell C'\) を計算する。
新しい source・gauge・陰関数解・jet・点・GB を計算しない。

## Source availability と unknown coefficient を落とさない理由

独立監査された native B11 packet は、actual proper 九座標
\((a,d,x,v_3,v_4,v_5,t,v_7,v_8)\) を全て coefficient algebra C44 の元として保存している。
これは source96 の実値と全30 gauge 方程式を有限 B11 で正確に満たす値であり、
first jet だけを座標関数の定義にしたものではない。

\[
 B_{11}=k[z_0,\ldots,z_6,r]/
 (r^3,z_iz_j,z_ir^2,z_2,z_3,z_6),
\]

基底は \((1,r,r^2,z_0,z_1,z_4,z_5,z_0r,z_1r,z_4r,z_5r)\)。C44はその上で
\((1,\tau_A,\tau_B,\tau_A\tau_B)\) を basis とし、保存された actual \(\nu_A,\nu_B\)
による二つの monic を使う。全1936 basis 積を別の monomial arithmetic で検査する。

既知の actual \(\bar C\) は、標的 \((\tau_A,\tau_B,\tau_A\tau_B)\) に対応する

\[
 k[r]/r^3\ \oplus\ k[r]/r^3\ \oplus\ k[r]/r^2
\]

であり、z作用は第一・第三で零、第二で \(z_i=-c_ir^2\)。
\(c_0=9835/864,c_4=29875/3456,c_5=-4565/864\)、その他は零である。
B11 の全 defining relations はこれら三環で零になる。特に欠けている z2,z3,z6 は
**全 sector の annihilator に入る**。z_i*rや全三次以上も同じ理由で消える。
従って B11 から三 sector への map は正確な環準同型で、未知高次係数を零と推測しない。
全121 B11 basis 積について三つの map の積互換性を確認する。

## 八つの実 nil-image と module closure

same proper coordinates から

\[
 h=a-d+x-2t,\quad c=(d+x)/2,\quad l=(d-x+a)/3,
 \qquad n_1=dh,\quad n_2=dcl
\]

を C44 内で計算する。これらの四 basis 倍は normalized nilideal 全体を R0-module として
張る。C44 の四成分のうち component0は元 D0 の unit image に含まれるため捨て、
三つの残り成分をそれぞれ上の sector に写す。これにより八つの actual 8-vector を得る。

\(\bar C\) は環ではなく module なので、その上で Hensel 計算を行わない。
有限 C44 環で計算した後に、実 scalar coefficient を module へ送るだけである。
八 vector の k-spanだけを像とはせず、保存済み八 R0-coordinate actions の全一次・二次
倍を加える。\(\mathfrak m^3\bar C=0\) のためこれで R0-submodule の全像になる。
得られた image の action-invariance と quotient の commuting actions を検査する。

この cokernel base change は right-exact なので

\[
 \bar C'=\operatorname{coker}(N_E/\ell N_E\longrightarrow\bar C)
\]

を正確に与える。kernel/H1 を mod-ideal で置き換える Tor shortcut は使わない。
全 vector、image の固定pivotとnormal form、quotient basis、全八作用、Hilbert function は
同名 report に保存する。原点の \(\mu_{R_0}C'=3\) と整合することも確認する。

この結果だけから元 nilideal の全生成数、K/JK、generic DVR rank、torsion exponent、
nil三乗の零性を推論しない。\(r\) は \(\ell\) で切った後の CA-F9 parameter であって、
比較 DVR parameter \(\ell\) とは異なる。

## 保存された実値からの結論

八つの nil-image と全 R0-multiples の image は正確に四次元で、8次元 \(\bar C\) の
座標 pivot は \((1,2,5,7)\)。従って

\[
 \boxed{\operatorname{length}\bar C'=4,\quad H_{\bar C'}=(3,1).}
\]

literal quotient basis は \((\tau_A,\tau_B,r\tau_B,\tau_A\tau_B)\) で、全 z_A 七座標は零作用。
残る作用は

\[
 r\tau_A=-\tfrac{11}{12}r\tau_B,\qquad r(\tau_B)=r\tau_B,\qquad
 r(r\tau_B)=r(\tau_A\tau_B)=0.
\]

よって新しい basis \((\tau_A+\tfrac{11}{12}\tau_B,\tau_A\tau_B,\tau_B,r\tau_B)\) により

\[
 \boxed{\bar C'\simeq k\oplus k\oplus k[[r]]/(r^2)}
\]

が S-module として正確に従う。有限値計算は全1936 C44 basis 積、全363 sector 積互換性、
360 image-closure columns の normal form、全八作用の不変性と quotient 上の可換性を検査した。
原点の3生成は保持され、追加 unit 消去ではなく higher module relations による長さの減少である。
generic DVR rank、torsion exponent および元 nil 三乗への追加の外挿はしない。
