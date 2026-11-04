# actual ell の source 一次形、正規化の regularity と零次数 nilideal の DVR-finiteness

2026-11-04。共通 G19 factor を除いた零次数実環 D0、その normalizer tensor E4、
\(N=\operatorname{Nil}(D_0)\)、\(N_E=\operatorname{Nil}(E_4)\) を対象とする。
DVR は \(D=k[[\ell]],\ell=5(s_A-r_A)\)。\(r_A\) はこれとは異なる F9 の座標。
新しい source 解、jet、点、CAS 計算は行わず、保存された二つの9×9一次変換だけを使う。

## 1. 同じ physical ell の proper source 上の一次形

独立監査済みの actual F9-to-physical9 行列により

\[
 x_{112}=\tfrac12r_A-s_A+O(2),\qquad
 x_{117}=-\tfrac38r_A-\tfrac34s_A+O(2).
\]

他の F9 七座標との混合も含む**全9×9行列**を解くと

\[
 \ell_1=-\tfrac{15}{2}x_{112}+\tfrac{10}{3}x_{117}.
\]

同じ monic source による physical R0-nine の regular C-nine への readback と、accepted
proper nine Jacobian を通す。全て同じ原点・同じ二つの trace basis における source であり、
得られる proper 一次形は

\[
 \boxed{\ell_1=-\tfrac56a+\tfrac53d+\tfrac{15}{2}x-\tfrac{10}{3}t.}
\]

五 pivot 方向の係数は全て零。proper u座標は pivot から accepted rational core roots を
引いたもので、roots は次数2以上なので、この一次形を変えない。

accepted 五直線の tangent に対する値は

| proper line | (a,d,x,t) tangent | ell一次値 |
|---|---|---:|
| N | (0,0,0,1) | -10/3 |
| C | (0,1,1,0) | 55/6 |
| R | (-2,1,-1,-2) | 5/2 |
| P | (2,0,-2,1) | -20 |
| L | (1,1,-1,-1/2) | -5 |

全て非零。したがって各 actual proper support line 上の全関数 ell は order1で、
その completed parameter ring は D と同型である。一次計算から未知高次係数を零とした
のではなく、非零一次係数による formal inverse theorem を使う。

main hypersurface(d)についても、ell一次形は d の倍数ではないので ell はその prime に
入らない。whole E4 の全 associated primes は(d)と五直線で尽くされているため

\[
 \boxed{\ell\text{ is a nonzerodivisor on BOTH whole }E_4\text{ and }E_{4,\rm red}.}
\]

reduced targetだけの regularity を whole E4へ推測したのではなく、六 associated primesを
全部確認した結論である。特に embedded N/P を落としていない。

## 2. normalized nil の正確な free DVR rank

N_Eは E4の部分 module。generic lengths は既存 full-source local modelsより

\[
 N:3,\quad C:3,\quad R:1,\quad P:2,\quad L:0,\quad(d):0.
\]

N/P の nilideal は core moduleへの d倍であるため長さ3/2であり、nil二乗の長さ1と
混同しない。Cでは長さ4のlocal algebraの maximal idealなので長さ3、Rでは双対数の
nilidealなので長さ1である。

support は一次元で、ellは全ての support curve 上の parameter。complete moduleの有限性
判定により N_Eは D上有限、かつ E4上の ell-injectivityにより D-torsionfree である。
各 curve のD-degree1と上表から

\[
 \boxed{N_E\simeq D^9}\quad\text{as a D-module, noncanonically}.
\]

同じ議論は既受理の nil-square moduleにも適用でき、\(N_E^2\) は D-free rank3となる。
これは S9 上の既知のCM1/multiplicity3結果と整合するが、source ellの transversalityを
別に証明して初めて得られる DVR上の結論である。

## 3. 元の nilideal の有限性と非標準的 module 分解

\(K=\ker(D_0\to E_4)\)、\(C=\operatorname{coker}(D_0\to E_4)\)、
\(C'=\operatorname{coker}(D_0\to E_{4,\rm red})\) とする。実 integral comparisonの核はnilで、
nil 元の preimageもnilだから、正確に

\[
 0\longrightarrow K\longrightarrow N\longrightarrow N_E
 \longrightarrow C\longrightarrow C'\longrightarrow0
\]

が成り立つ。N_EとKは D-finiteなので NもD-finiteである。
受理済み actual ranks \(\operatorname{rank}_D K=4\)、\(\operatorname{rank}_D C=4\)、
\(\operatorname{rank}_D C'=2\) により

\[
 P:=\operatorname{image}(N\to N_E)\text{ has rank }9-4+2=7.
\]

PはD-free module N_Eの有限部分なのでD-free7。従って

\[
 0\to K\to N\to P\to0
\]

は D-moduleとして分裂し

\[
 \boxed{N\simeq K\oplus D^7,\qquad\operatorname{rank}_D N=11.}
\]

これは **非標準的な D-module 分解のみ**。R0-module、D0-module、idealの積、または環の
分解として主張しない。具体的な七つの nil lift はまだ構成していない。

Kについて既受理の \(\mu_D K=L-4\)、\(\mu_D\operatorname{tors}K=L-8\) を使うと

\[
 \boxed{\mu_D N=L+3,\qquad\mu_D\operatorname{tors}N=L-8.}
\]

新しい actual third-layer source peer は \(11\le L\le17\) を与えるため、その区間は14〜20と
3〜9である。また \(\operatorname{tors}_D N=\operatorname{tors}_D K=H^0(D_0)\) は実 inclusion
の下の等号である。torsionの全長・指数は計算しない。

## 4. honest nil-comparison への利用と範囲

ellが Ered上 regular、Nが DVR-finiteという、honest nil-comparison specialized T_nil に
必要な二つの前提が閉じた。実 source chain section/left inverse が完全に構成されれば

\[
 0\to N/\ell N\to T_{\rm nil}\to C'[\ell]\to0,
 \qquad\operatorname{length}T_{\rm nil}=L+5
\]

となる。ここで \(\dim C'[\ell]=2\) を使う。ただし新しい104×948 source presentationを
このノートで作ったとはしない。既存 finite imageやmod-J homologyをhonest T_nilとする
Tor shortcutも行わない。

NのD-finitenessを \(I_AI_BN=0\) と取り違えない。Rの nil部分はその conductor 積の外にも
あり得る。N_EがD-freeであること、また上のmodule分裂から、元 nil三乗や正確なnil指数は
決まらない。全正次数nilを含む full graded ring 全体がD-finiteという主張でもない。
