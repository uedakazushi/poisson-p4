# 実 raw cokernel の閉ファイバー作用：Artin20 評価の source・費用 preflight

2026-11-04、v63。本 packet は方法と入力の検査であり、新しい有限代数値 source solve は実行しない。

直前の実 conductor 積の閉ファイバーを
\(A_{20}=R_0/(I_A I_B,\ell)\) とする。実 CA F9 座標
\((z_{A,0},\ldots,z_{A,6},r_A)\) の全乗法表が確定し、
\(\mathfrak m^5A_{20}=0\) である。raw physical 座標関数と両 normalizer の係数をこの
有限代数で直接求めれば、未保存の全九変数三次 jet を補う必要がない。

## 同じ実 source に対する単位問題

凍結 joint105 は CA の57変数系（pure48 と実 F9 の九式）と CB の pure48 から成る。
その九 target を

\[
 (z_{A,0},\ldots,z_{A,6},r_A,r_A)\in A_{20}^9
\]

とする。raw physical 九座標は入力として別に推測せず、この105未知数の一部として解く。
元の両 pure48 と既受理の normalizer free4 packet の pure48 が、保存済み pruning map
を通じて未知数・式の両方で完全に同じであることを確認した。従って得られる
\(\nu_A,\nu_B,a_{A,j},a_{B,j},b_{A,j},b_{B,j}\) は同じ actual normalizer 係数である。
全 source 変数の中心は零、必要な六つの逆元の中心は全て1である。

原点 Jacobian の逆行列は既存の固定105逆行列であり、非零係数は866個である。
固定原点逆行列による補正は、一般に誤差次数を一つ上げる。自動的な倍増とはしない。
残差が \(\mathfrak m^d\) にあれば

\[
 u\longmapsto u-J_0^{-1}(F(u)-\text{target})
\]

で残差は \(\mathfrak m^{d+1}\) に入るため、零初期値から四回の補正で足りる。
その後、全105式を (A_{20}) で正確に確認する。
単位逆元は nilpotence を用いる五項の幾何級数で正確に評価できる。

## 実測した入力と上界

全 source 41,315 nodes のうち、必要な105式・12組の係数・両 conductor を覆う閉包は
18,584 nodes である。積10,608、和7,669、符号151、変数105、定数45、逆元6。
参照回数で値を解放すれば同時に生きる node 値は414個であり、未知数105個を別保持しても
20成分ベクトルの係数上界は10,380個である。

20個の基底 monomial の400組の積のうち、非零は80組だけであり、各積の support は1である。
従って四回の全 source 評価の単純上界は約340万個の有理係数積である。
四段階の未知数・式の全ベクトルを保存しても16,800係数で済む。
source 閉包・定数 Jacobian・20表・必要 roots の全 hash と順序を証拠に保存した。
今後の一回の評価は60秒・2 GB・出力32 MB/150,000係数を guard とし、失敗時に自動拡張しない。

## raw26 と実 cokernel への接続

別担当の raw26 座標 packet は

\[
 C_{\rm raw}/(\ell)=
 (A_{20}/I_A)\tau_A\ \oplus\
 (A_{20}/I_B)\tau_B\ \oplus\ A_{20}\tau_A\tau_B
\]

の実26基底と、二つの20→3商写像・その代表持上げを構成する。
各実 physical difference
\(F_j=\delta_j+b_{A,j}\tau_A-b_{B,j}\tau_B\) の作用は、その座標で

\[
 \begin{pmatrix}
 \delta_j&0&-\nu_Bb_{B,j}\\
 0&\delta_j&\nu_Ab_{A,j}\\
 -b_{B,j}&b_{A,j}&\delta_j
 \end{pmatrix}
\]

となる。成分間の積は実商写像・代表持上げを伴う。
持上げ依存性が消えるのは \(I_AI_B=0\) と (b_A\in I_A,b_B\in I_B\) による。
12個の26×26作用行列が得られれば、右完全性により実 (C/(\ell)) はその像の和で割った
cokernel である。特殊化した kernel を実比較 kernel と同一視しない。

この packet は Artin 値 solve、12作用、実 (C/(\ell)) の次元、一般 DVR rank、
比較 kernel の torsion・深さのいずれも計算したとは主張しない。
