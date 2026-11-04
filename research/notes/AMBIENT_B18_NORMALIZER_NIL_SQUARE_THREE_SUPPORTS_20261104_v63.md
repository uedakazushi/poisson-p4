# normalized E4 の nil 二乗：三つの支持と正確な reduced annihilator

2026-11-04。対象は既受理の実 proper-source 正規化比較環 (E_4\) だけであり、
元の (D_0\) または全交叉 D の nilideal と同一視しない。新しい source・jet・点・CAS
計算を用いず、whole E4 の五直線模型から nil-square module を特定する。

## 1. core と whole E4 の積を区別する

\[
 S=k[[u_1,\ldots,u_5,a,d,x,t]],\quad
 J=(u_1,\ldots,u_5,Q_0,Q_1,(a+x)(x-a-d)^2),\quad E_4=S/(dJ).
\]

core (T=S/J\) は既受理の一次元 CI である。whole E4 の associated primes は

\[
 (d),\quad \mathfrak p_N,\mathfrak p_C,\mathfrak p_R,
 \mathfrak p_P,\mathfrak p_L
\]

の六つで、N,P は embedded である。記号 N 等は既存 source の座標直線ラベルであって、
旧 physical strata への新しい同定ではない。

N,P で (E_4\) の nilradical は ((d)/(dJ)\) であり、乗法 (d\) は module 同型

\[
 T_{\mathfrak p}\ \xrightarrow{\ d\ }\operatorname{Nil}(E_4)_{\mathfrak p}
\]

を与える。しかし **環同型ではない**。この同型の下で

\[
 \operatorname{Nil}(E_4)_{\mathfrak p}^2\ \leftrightarrow\ dT_{\mathfrak p},
 \qquad
 \operatorname{Nil}(E_4)_{\mathfrak p}^3\ \leftrightarrow\ d^2T_{\mathfrak p}.
\]

従って core で (d\ne0,d^2=0\) なら、whole E4 では nil 二乗は**非零**、三乗は零である。
core の (d^2=0\) を whole E4 の nil 二乗零へ読み替えてはいけない。

## 2. nil 二乗の全 generic 長さ

既受理の実 Laurent source 模型をそのまま用いる。

| E4 の associated prime | generic source | nil-square の長さ |
|---|---|---:|
| N | (T_N=\kappa[[v]]/(v^3),\ d=(4/t)v^2\) | 1 |
| P | (T_P=\kappa[[v]]/(v^2),\ d=v\) | 1 |
| C | (d\) は単元、core HF (=(1,2,1)\) | 1 |
| R | (d\) は単元、core は双対数 | 0 |
| L | (d\) は単元、core は被約一点 | 0 |
| ((d)\) | whole E4 の主成分の一般点は被約 | 0 |

N では (dT_N=\kappa v^2\)、P では (dT_P=\kappa v\) である。
C では Artin core の maximal ideal の二乗が HF の最後の一次元部分である。
全ての長さは対応する generic residue field 上の長さであり、新しい数値点を使わない。


\(M=\operatorname{Nil}(E_4)^2\) は E4 の部分 module だから
\(\operatorname{Ass}M\subset\operatorname{Ass}E_4\) である。上表から正確に

\[
 \operatorname{Ass}_S M=\{\mathfrak p_N,\mathfrak p_P,\mathfrak p_C\}.
\]

原点 prime は含まれない。従って M は一次元 Cohen--Macaulay module であり、
有限長の原点 submodule を持たない。

## 3. annihilator と重複度

三直線の actual 座標は

\[
 N=(0,0,0,t),\quad P=(2t,0,-2t,t),\quad C=(0,t,t,0)
\]

である。次の有理線形座標変換は行列式 (1/2\) を持つ。

\[
 c=d,\quad p=a/2,\quad n=t-a/2,\quad H=a-d+x.
\]

N,P,C はそれぞれ n,p,c 軸となる。従って三つの reduced line primes の交叉は

\[
 I_{NPC}=(u_1,\ldots,u_5,H,cp,cn,pn)
       =(u_5,\ a-d+x,\ ad,\ d(2t-a),\ a(2t-a)).
\]

最後の (u_5\) は五変数 ideal の略記である。上表で M の各 generic localization は
長さ1なので、それぞれの maximal ideal が M を殺す。従って (I_{NPC}M\) は三 generic
点で零となる。非零ならその associated prime は M の上記三点の一つにならざるを得ず、
矛盾する。よって (I_{NPC}M=0\)。逆に M は三点全てで非零なので、その annihilator は
各 line prime に含まれる。従って

\[
 \boxed{\operatorname{Ann}_S\operatorname{Nil}(E_4)^2=I_{NPC}.}
\]

M は reduced 三直線環上の torsion-free module で、各枝の generic rank は1である。
三直線はいずれも原点重複度1なので、associativity formula により

\[
 \dim M=\operatorname{depth}M=1,\qquad e_{\mathfrak m_S}(M)=1+1+1=3.
\]

これは M が三直線環そのもの、free module、または canonical module と同型という主張
ではない。M の最小生成数や具体的な枝間 gluing はここでは求めない。

## 4. 元 D0 の上界8に必要な corrected gate

元の nilideal (N_0\) の像から生成される E4-ideal

\[
 U=E_4\cdot\operatorname{image}(N_0^2)\subset M
\]

を考える。上の完全な associated-prime 記述から、

\[
 U=0\quad\Longleftrightarrow\quad
 U_{\mathfrak p_N}=U_{\mathfrak p_P}=U_{\mathfrak p_C}=0.
\]

この三つの**実 source 比較写像の像**の平方零性を証明できれば、全域で
\(N_0^2\subset K\) となり、既受理の (K^4=0\) から (N_0^8=0\) を得る。
Qb の Z4 branch を C と source 同定しても、それだけでは N,P の二つの gate は残る。
N,P の target nil-square 自体は上表の通り非零である。

本稿はこれらの source-image gate を閉じず、上界8を主張しない。
既受理の無条件上界12、下界3、元 nil 三乗の未決性を保持する。
全域 associated closure、成分数、元 D0 の nil 生成元や深さについて新しい主張はしない。
