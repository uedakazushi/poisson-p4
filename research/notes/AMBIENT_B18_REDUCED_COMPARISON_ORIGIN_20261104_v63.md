# reduced target 比較の原点には追加 unit cancellation がない

2026-11-04。対象は actual \(C'_{\rm red}=\operatorname{coker}(D_0\to E_{\rm red})\)
の (R_0\)-module 原点 fiber。全 nil 生成元、nil 三乗、K/JK は計算しない。

accepted proper-to-physical interface の Jacobian を、同じ R0 の全九函数を零にした
二つの trace 方向へ制限する。proper source の retained9 座標順は
\((a,d,x=v_2,v_3,v_4,v_5,t=v_6,v_7,v_8)\)。保存された実一次係数は

\[
 a=\tau_A-\tau_B,\quad d=-2\tau_A-\tfrac12\tau_B,
 \quad x=\tau_A,\quad t=\tau_A
\]

であり、五つの pivot \(v_4,v_3,v_5,v_7,v_8\) は trace の一次部分を持たない。
accepted actual core elimination の五つの rational root の numerator は全て次数2以上で、
denominator は D/P の原点単位である。従って proper
\(u_i=v_{\rm pivot,i}-\operatorname{root}_i(a,d,x,t)\) の trace 一次部分も全て零。

同じ coefficient algebra は

\[
 C/\mathfrak m_{R_0}C=k[\tau_A,\tau_B]/(\tau_A^2,\tau_B^2).
\]

この fiber の maximal ideal の三乗は零であるため、\(\mathfrak m_C^2\) に入る native
generator を評価するには上記一次部分だけで十分である。高次 source 項を零と仮定した
のではなく、それらの寄与は三乗 maximal ideal へ入る。

\(h=a-d+x-2t\) の一次部分は \(2\tau_A-\tfrac12\tau_B\) なので

\[
 dh=-4\tau_A^2+\tfrac14\tau_B^2=0
 \quad\text{in }C/\mathfrak m_{R_0}C.
\]

\(du_5\) は零であり、残る \(dcr,dcl,drl\) は三乗 maximal ideal に入る。
従って reduced target の全9 native generators \(dR_9\) はこの fiber で零となり、
その四 basis 倍、つまり全36 Qrad 列の R0-origin 四成分も全て零。

actual C0 の odd 八生成元の係数 \(\beta=(z_{A,7},g_A)\) は原点で零であるため、
unit1 消去後の実比較行列

\[
 V_{\rm rad}:R_0^{44}\to R_0^3
\]

の原点行列は零で、rankは**0**である。従って

\[
 \boxed{\mu_{R_0}C'_{\rm red}=3.}
\]

rank1を期待して二生成へ unit 消去する経路はここでは成立しない。
この結論は cokernel の right-exact fiber だけを使う。kernel や homology を安易に
base change したり、元 nilideal の生成数をこの3から推論したりしない。

次の有限 ell=0 quotient は依然
\(C'_{\rm red}/\ell C'_{\rm red}=\operatorname{coker}
(N_E/\ell N_E\to C/\ell C)\) である。既存 Cbar の長さ8から長さは8以下だが、
二つの normalized nil generator の四 basis 倍の**実係数値**はまだ求めていない。
原点で全係数が零でも、この八つの有限 module class が零とは限らない。
