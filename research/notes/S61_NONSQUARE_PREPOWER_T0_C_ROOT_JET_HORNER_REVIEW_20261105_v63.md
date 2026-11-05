# T0のC valuation四：独立five-jet Horner受理

2026-11-05。元のprimary note/reportは変更しない。
ROOT独立report `S61_NONSQUARE_PREPOWER_T0_C_ROOT_JET_HORNER_REVIEW_20261105_v63`
（SHA256 `60be3df6763cbfb294727209de1d2fd0bfd672b21488bd0d6753a47a83f2140f`）を結ぶ。

旧root239のoriginal QQ7208係数、a-offset4を一度だけ除いたU0=T0/a⁴を使い、
producerのbinomial／R powers／native engineをimportせずにtruncated Horner五jetを計算した。
R=2+9e−18e²、stageごとの共通denominator D*100^(20-i)を保持し、
`Jnew_r=R*J_r+100*J_(r-1)`のcrossjet倍率100も含めた。

**下四層は全係数零、k=4の全646係数とdenominatorが完全一致**した。
したがってv_C(T0)=v_C(U0)=4を独立に確定した。
独立計算109781 products、実432 bits、owned39914。
高いk層、他因子・tail、商、GCDはこの監査では計算していない。

CはFULLU既知単元であり、この結果だけで不良集合やD12の補集合を消さない。
次の固定C四乗商は別SOURCE／code／budget監査と一回GOを要する。
Gaussによる元D197での整数化をその準備に使えるが、まだ商値は未計算である。
