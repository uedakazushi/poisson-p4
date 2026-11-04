# 実 conductor と全 nil-cycle の even 係数による nil 指数上界7

2026-11-04。対象は元の実交叉 (D=CA\cap CB\) の B18 代表点 (P_0\) の完成環。
新しい source・jet・点・CAS 計算を行わず、既存 cached source と完全比較から

\[
 \boxed{\operatorname{Nil}(\widehat{\mathcal O}_{D,P_0})^7=0},
 \qquad 3\le\text{nil index}\le7
\]

を得る。旧上界12の保存 packet は変更しない。三乗零・三乗非零・正確な指数は未決。
正規化 E4 の embedded N/P で nil 二乗は非零であり、それを零とする議論は使わない。

## 1. 同じ physical 底上の正確な ideal

零次数実比較は

\[
 A_i=R_0\oplus I_i\tau_i\subset B_i=R_0\oplus R_0\tau_i,
 \qquad \tau_i^2=\nu_i\in R_0,\quad i=A,B
\]

から得る。ここで (R_0\) は同じ九つの physical 函数による正則局所環であり、
\(I_i\) は実 zero-odd ideal。二族の trace を同一視しない。

\[
 D_0=A_A\widehat\otimes_{P_0}A_B,
 \qquad E_4=B_A\widehat\otimes_{P_0}B_B,
 \qquad K=\ker(D_0\to E_4).
\]

既存 actual joint105 source と独立監査は、CA の unit 座標
\((z_{A,1},\ldots,z_{A,7},r,s)\) 上で全次数の等式

\[
 J:=I_A+I_B=(z_{A,1},\ldots,z_{A,7},r^2,rs,s^2)
\]

を証明する。CB restriction の全二次係数を張る minor は625/768で、Nakayama を
\((r,s)^2/J\) に適用した厳密同定である。従って

\[
 \mathfrak m_{R_0}^2\subset J.
\]

## 2. trace 係数の小さい actual source gate

元41315-node source の実根は \(\nu_A=20086,\nu_B=40865\)。受理済み joint105
の CA zero7 平面 \(z_A=0\) 上の保存 jet を使い、各 root の一次係数だけを再解釈する。
新しい陰関数解を計算しない。二つを合わせた依存 closure は11 nodeで、結果は

\[
 (\nu_A(0),\partial_r\nu_A(0),\partial_s\nu_A(0))=(0,0,0),
\]

\[
 (\nu_B(0),\partial_r\nu_B(0),\partial_s\nu_B(0))=(0,0,0).
\]

joint source prefix が元 source と byte-literal に同じであること、unknown binding、
同じ physical 底と元 root ID を確認し、全105式を独立監査した既存 receipt に束縛する。
したがってこれは別々の形式座標を混ぜた零性ではない。

各関数は \((z_A)+\mathfrak m_{R_0}^2\) に属する。前節と合わせて全次数で

\[
 \boxed{\nu_A,\nu_B\in J.}
\]

有限 jet から関数そのものの零を推論していない。必要なのは二乗 maximal ideal を
含む既知 ideal への所属だけなので、一次係数で判定が完結する。

## 3. あらゆる実 nil 元の even 係数は J に入る

actual normalizer coefficient algebra (C\) は同じ (R_0\) 上の free4 環

\[
 C=R_0[\tau_A,\tau_B]/(\tau_A^2-\nu_A,\tau_B^2-\nu_B)
\]

であり、同じ原点で \(\tau_A=\tau_B=0\)。その正則九座標は accepted proper physical
interface により指定されている。reduced target の関係は

\[
 dR_9=d(u_1,\ldots,u_5,h,cr,cl,rl).
\]

\(d\) と九つの \(R_9\) 函数は全てこの同じ原点の maximal ideal に入る。
従って36列全部、すなわち各 (dR_j\) の四 basis 倍も \(\mathfrak m_C^2\) に入る。

任意の (f,g\in\mathfrak m_C\) を四 basis で書くと、constant 成分は

\[
 (fg)_0=f_0g_0+\nu_Af_Ag_A+\nu_Bf_Bg_B+
             \nu_A\nu_Bf_{AB}g_{AB}.
\]

\(f_0,g_0\in\mathfrak m_{R_0}\) なので

\[
 \operatorname{pr}_0(\mathfrak m_C^2)
 \subset\mathfrak m_{R_0}^2+(\nu_A,\nu_B)\subset J.
\]

ここで原点、trace basis、component projection は実 coefficient interface と同じである。
ゆえに actual (Q_{\rm rad}\) の even row 全36係数は J に属する。

完全 radical comparison は **全ての** (y\in N_0:=\operatorname{Nil}(D_0)\) を
cycle (z=(a_{36},c_8)\) で表す。その実 C0 代表元は

\[
 y=b_0+\sum_{i=1}^8c_i\beta_i\tau_A\pmod{\text{actual D0 relations}},
 \qquad b_0=(Q_{\rm rad})_0a\in J,
\]

\(\beta=(z_{A,1},\ldots,z_{A,7},g_A)\)。旧 (g_A\) も実 (I_A\) の生成元である。
特定の候補 cycle だけでなく全 nilideal について、この even 所属が従う。

## 4. 二段核と ENTIRE conductor

一般二次包含の任意基底変換核が平方零であることは旧上界12 packet で証明済み。
先にAを正規化する写像の核を (K_A\)、先にBを正規化する核を (K_B\) とする。
どちらの順序でも第二段の核が平方零なので

\[
 K^2\subset K_A\cap K_B.
\]

各 inclusion の conductor には **全 ideal** (I_iB_i\subset A_i\) が入る。
任意の \(c\in I_iB_i\) に対し、乗法 \(c:B_i\to A_i\) は (A_i\)-linear であり、
inclusion との合成は c 倍である。そのため任意の tensor base change 後の inclusion
核も c に annihilate される。有限 complete module でも同じ連続写像である。

従って (I_AK_A=I_BK_B=0\)、ゆえに (JK^2=0\)。加えて各
\(\beta_i\tau_A\in I_AB_A\) は (K_A\) を殺すので (K^2\) も殺す。
前節の **全 nil 元** の表示を用いると

\[
 \boxed{N_0K^2=0.}
\]

ここで scalar ideal (I_A\) のみで止めると odd 項を殺せないため、whole conductor
による factorization が不可欠である。Tor を捨てた単射の仮定は使っていない。

## 5. 上界7と範囲

whole (E_4\) の既受理の nil 三乗零より

\[
 N_0^3\subset K,\qquad N_0^6\subset K^2,
 \qquad \boxed{N_0^7=0}.
\]

また K は nil ideal なので前節から (K^3=0\) も従う。K の正確な指数は主張しない。
既受理の actual stable tail、全正次数 torsion 零、完成閉性の議論により、
全 nil の二個以上の積は零次数に限られる。従って全横断 D の nil 七乗も零である。
actual (G_{19}\) の正則群因子と projective scalar chart もこの指数を保存する。
Qb family の coherent nil-square 閉包による下界3と合わせ、元 B18 germ で

\[
 3\le\text{nil index}\le7.
\]

本証明は E4 の N/P/C における像を零と推測せず、五枝の物理ラベル同定も使わない。
深さ19や全 (H^0\) の長さ・annihilator は不要であり、旧 Tor63 image を全 H0 としない。
元 nil 三乗の零・非零、正確な指数、全 nil 生成元と全域分類はなお未決である。
