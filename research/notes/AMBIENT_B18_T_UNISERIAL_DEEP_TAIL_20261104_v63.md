# 実比較加群 T の深い部分は一つの r-chain である

2026-11-04。対象は実特殊比較加群 (T\)。係数環は実切断

\[
 S=R_0/(\ell)\simeq k[[z_{A,0},\ldots,z_{A,6},r_A]],\quad \mathfrak m=\mathfrak m_S.
\]

以下 (r=r_A\) と略す。これは (T\) の coefficient cut に残る actual CA-F9 座標であり、
比較 (K\) を有限加群とする DVR parameter \(\ell\) ではない。
\(e_4\) は実最小五生成表示の最後の生成元で、source V56-row55 の標的である。

既受理の実二次 rank は \(\operatorname{HF}_{\mathfrak m}T=(5,3,2,\ldots)\) を与える。
新しい実三次関係を含む40生成の既知初期部分加群には、独立監査された module GB がある。
その商の次数3以上の標準単項式は一つ \(e_4r^d\) だけで、別の実 base 証明は
\(\mathfrak m^{10}T=0\) を与える。既知初期部分加群を全初期加群と同一視しない。

## 1. 七つの固定 reduction が与える追加情報

標準単項式が一つという事実だけでは、\(z_i\) の作用が零とは限らない。
そこで実40生成を使い、七つの単項式

\[
 e_4z_{A,i}r^3,\qquad 0\le i\le6
\]

を position-over-lex 順序で別の sparse Fraction reduction により確認した。
全て零まで還元され、その40生成に関する厳密係数 trace を保存した。
これを (r^{d-3}\) 倍すれば、全 (d\ge3\) で既知初期部分加群の商に

\[
 z_{A,i}(e_4r^d)=0
\]

が成り立つ。従って実 T では

\[
 z_{A,i}(r^3e_4)\in\mathfrak m^5T.
\]

これは全次数で (z_i\) の作用そのものが零という主張ではなく、深い部分での初期作用の零性である。

## 2. Nakayama による cyclic な深い加群

\(M=\mathfrak m^3T\)、\(g=r^3e_4\) と置く。
\(M/\mathfrak mM\) は g の class 一つで張られるので、Nakayama より (M=Sg\)。
この class が零なら (M=0\) であり、後の (s=0\) の場合になる。

次数4の実 graded piece は \(rg=r^4e_4\) で張られるので

\[
 \mathfrak mM=rM+\mathfrak m^2M.
\]

有限 module \(\mathfrak mM/rM\) に Nakayama を適用すると

\[
 \boxed{\mathfrak mM=rM},\qquad
 \mathfrak m^jM=r^jM\quad(j\ge0).
\]

よって \(M=kg+rM\)。T は有限長なので、この式を有限回反復すれば M は
\(g,rg,r^2g,\ldots\) の k-linear span となる。従って

\[
 \boxed{M\simeq k[[r]]/(r^s)},\qquad 0\le s\le7.
\]

上界7は \(r^7M=\mathfrak m^{10}T=0\) による。(s\) は M の長さである。
全 S-submodule は r に安定なので、この部分は S-module としても uniserial である。

## 3. 残る七作用と正確な Hilbert 関数の形

\(s>0\) なら、各 \(z_{A,i}\) は r と可換な M の endomorphism だから、
一意の \(h_i(r)\in k[[r]]/(r^s)\) による乗法である。
第1節と \(\mathfrak m^2M=r^2M\) より

\[
 h_i(r)\in r^2k[[r]]/(r^s),\qquad 0\le i\le6.
\]

\(s=0\) では M は零で、この作用は全て自明である。
どの (h_i\) の係数も本稿では計算しない。これは元 T 全体を一変数加群にした主張でもない。

以上により、未知の全 Hilbert 関数は一つの整数 (s\) だけで決まり、

\[
 H_T(t)=5+3t+2t^2+t^3\frac{1-t^s}{1-t},\qquad
 \boxed{\operatorname{length}T=10+s},\quad 0\le s\le7.
\]

すなわち深い layer は1が連続してから零となり、途中で消えた後に再出現しない。
正確な長さの次の判定は、実元 \(g=r^3e_4\) の r-annihilator を求めることである。
この元の追加 relation を、原点の matrix の零性や ambient truncation の membership から
推測しない。新しい source coefficient・jet・GB探索は行っていない。

(s\) は DVR \(k[[\ell]]\) 上の K の torsion exponent、(H^0(D_0)\) の長さ、
元 D の nil 指数とは別の量である。これらの同一視や nil-cube の判定は主張しない。
同名 report に exact source gate、七 reduction、標準単項式と全証明を保存する。
