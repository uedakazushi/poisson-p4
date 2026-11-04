# conductor transfer の実512 nil-cycle と有限 K/JK pairing の正確な残余

2026-11-04。元 B18 の零次数実交叉 (D_0\) について、新しい source 解・jet・点・
scalar coefficient 展開を行わず、既存 normalizer 比較から実 nil-cycle を構成する。
全 nilideal の生成元を求めたとはしない。

## 1. K/JK の有限性と、なお必要な honest presentation

\(K=\ker(D_0\to E_4)\)、\(J=I_A+I_B\) とする。既受理の上界7の証明は

\[
 JK^2=0,\qquad N K^2=0,\qquad N=\operatorname{Nil}(D_0)
\]

を与える。従って実 pairing \(K\times K\to K^2\) は両引数について K/JK に因子化する。
\(J\supset\mathfrak m_{R_0}^2\) なので \(\ell^2\in J\)。一方、実比較の DVR
\(D=k[[\ell]]\) 上で \(\operatorname{rank}K=4\)、\(\mu_DK=L-4\le13\) だから

\[
 \dim_k K/JK\le\dim_k K/\ell^2K
 =8+\sum_{i=1}^{\mu_DK-4}\min(2,e_i)\le26.
\]

ここで \(e_i\) は未計算の torsion exponent であり、値を仮定しない。
\(K^2\) は D/\(\ell^2\)-module で、D上の K の13生成元の対積で張られるため
\(\operatorname{length}K^2\le2\binom{14}{2}=182\) も従う。
この長さ上界および二乗 maximal ideal の annihilation は K² とその部分 N⁶ のみであり、
**全 H0 や N³の長さ・annihilator の上界ではない**。

ただし K/JK の完全表示は未構成である。完全な \(Q_{48}\)、\(\operatorname{Syz}Q_{128}\)、
\(UP=QB\) は target と比較複体を指定するが、\(\operatorname{Syz}V\) の全 source lift を
与えるわけではない。(A=R_0/J\)、比較 cokernel を C と書くと、free comparison complex
\(\mathcal F\) について edge sequence は

\[
 \operatorname{Tor}_2^{R_0}(C,A)\longrightarrow K/JK
 \longrightarrow H_1(\mathcal F\otimes A)
 \longrightarrow\operatorname{Tor}_1^{R_0}(C,A)\longrightarrow0.
\]

左端の connecting image を捨てない。現在の有限 mod-J homology を K/JK と同一視する
ことはできない。実 full cycle lift とこの connecting map の完全な補正、または honest な
K-presentation が必要である。(J=(z_7,(r,s)^2)\) から \(\operatorname{length}R_0/J^2=31\)
となるが、その厚みだけで K/JK を回復できるという Artin--Rees strictness は証明していない。

以下はこの未解決 gate を仮定せずに作れる、実 nil 乗法の別の source operator である。

## 2. 全 conductor による transfer

実二次包含 \(A_i\subset B_i\) と \(a\in I_A,b\in I_B\) に対し、乗法の tensor で

\[
 \phi_{a,b}:E_4=B_A\widehat\otimes_{P_0}B_B\longrightarrow
 D_0=A_A\widehat\otimes_{P_0}A_B
\]

を作る。これは D0-linear であり、\(f=ab\)、\(\psi:D_0\to E_4\) とすると

\[
 \psi\phi_{a,b}=f\,\mathrm{id}_{E_4},\qquad
 \phi_{a,b}\psi=f\,\mathrm{id}_{D_0},
\]

\[
 \phi_{a,b}(x)\phi_{a,b}(y)=f\phi_{a,b}(xy).
\]

写像は組 (a,b) により定める。同じ積 f の異なる分解で transfer が等しいとは仮定しない。
whole E4 の nilideal を \(N_E\) とすれば \(N_E^3=0\)。従って

\[
 W:=\sum_{a\in I_A,b\in I_B}\phi_{a,b}(N_E)\subset N,
 \qquad \boxed{WN^2=0},\qquad \boxed{(I_AI_B)N\subset W}.
\]

前者は D0-linearity による
\(\phi_{a,b}(x)nn'=\phi_{a,b}(x\psi(n)\psi(n'))=0\) から従う。
後者は \(abn=\phi_{a,b}(\psi(n))\) から従う。従って W³=0であり、元の nil 三乗写像は
各引数について **N/W** に因子化する。N/W は (I_AI_B\) に annihilate されるが、有限長
とはしていない。W、N、K、全 H0、旧 Tor63 image の同一視は行わない。

## 3. source-bound な512生成元

CA の実 ideal basis を \(\beta=(z_{A,1},\ldots,z_{A,7},g_A)\) とする。
CB は八つの raw zero odd rows

\[
 j=(18,19,28,29,60,67,72,63)
\]

の係数 \(b_j\) を使う。row63 と (g_B\) の間の既受理 unipotent 変更により、これらは
全 (I_B\) を生成する。旧 gA/gB と r/g functions の source は変更しない。
normalized nilideal の実二生成元は

\[
 n_1=dh,\qquad n_2=dcl.
\]

これらは actual radical source の第5・第7 native generator（0開始）である。
E4 は D0-module として \(e_\gamma=(1,\tau_A,\tau_B,\tau_A\tau_B)\) で張られる。
従って W は

\[
 \phi_{\beta_i,b_j}(n_\alpha e_\gamma),
 \quad 1\le i\le8,\ 1\le j\le8,\ \alpha=1,2,\ \gamma=0,1,2,3
\]

の **512個**で D0-module/ideal として生成される。最小性や非零性は主張しない。
有限表示の nil 元を512個選んだだけという意味ではなく、上で定義した transfer 部分 W の
全生成元である。一方で、元の N 全体の生成元という意味ではない。

## 4. literal C0 representative と complete H1 cycle

実 physical 差の最初の八本は、同じ C-basis で

\[
 F_j=\delta_j+b_{A,j}\tau_A-b_j\tau_B
     =\sum_{g=0}^8 A_{jg}\,dR_g.
\]

符号と順序は accepted `actual_F12_C4`、`native_new_A12_by9_nodes` に束縛する。
(p=n_\alpha e_\gamma=p_0+p_A\tau_A+p_B\tau_B+p_{AB}\tau_A\tau_B\) と置く。
transfer 元の実 C0 free9 representative は、even 成分

\[
 E=\beta_i\bigl(b_jp_0+\delta_jp_B+\nu_A b_{A,j}p_{AB}\bigr),
\]

および **第i odd 成分だけ**

\[
 H=b_jp_A+b_{A,j}p_B+\delta_jp_{AB}
\]

を持つ。実 C0 の像を U とすると、coefficient algebra C 内で正確に

\[
 U(E,H)=\beta_i b_jp+\beta_iF_j(p_B+\tau_Ap_{AB}).
\]

この四成分恒等式を独立した小さい普遍 polynomial interpreter で確認する。
従って complete radical comparison の36係数を、九つの C-coefficient として

\[
 a_g=\beta_iA_{jg}(p_B+\tau_Ap_{AB})
       +\mathbf1_{g=g_\alpha}\,\beta_i b_j e_\gamma
\]

と置けば \(Q_{\rm rad}a=U(E,H)\) が全 source 上で成立する。
四成分 restriction した (a_{36}\) と sparse (c_8=H e_i\) は正確に
**44成分の実 cycle** であり、消去した even 成分も \((Q_{\rm rad})_0a=E\) に戻る。
512本全ての root ID と source row/normalizer frame、全 nine C-coefficient を保存する。

780列への boundary 表示を推測したり、これらの cycle が N を尽くすと主張したりしない。
関係 \(F=A(dR)\) と36列の source coefficient interface は独立監査済みであるため、
これらの operator は未指定の未知係数を含む recipe ではなく、既存 actual algebraic
functions に対する finite typed DAG である。

出力には既知の512転送 cycle と未知の K/JK presentation gate を区別して保存する。
元の nil 三乗、W以外の nil 生成元、K²、正確な指数について数値的結論は出さない。
