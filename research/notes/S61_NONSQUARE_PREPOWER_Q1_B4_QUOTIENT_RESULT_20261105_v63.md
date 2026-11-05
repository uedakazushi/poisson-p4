# 原共通分母を保つ固定 B⁴ 商：2317項の確定結果

2026-11-05。一次報告
`S61_NONSQUARE_PREPOWER_Q1_B4_QUOTIENT_ENTRY_20261105_v63`
SHA `24370d35886931b1c636e5d4b421fe01bf67924684ba22ff13235b6386719cc5`。

\[
Q_2=Q_1/B^4=T_1/(a^3C^4B^4)=J_2/D
\]

を2317項の整数J2と元D197で保存した。このQ2は元S61の二次変数q2ではない。
四つの固定lex50ad除算は各全余り零で、項数は4087→3432→2843→2317。
固定B⁴は74項。原入力をimmutable fileから再読し、全係数の \(B^4J_2=I_1\) を確認し、
全保存係数の読戻しも通った。

peak30404/52826 slots、integer247922・LC50逆乗算12679、RAW1192875 bytes、約0.070秒。
J2の最大係数269 bit、整数L1 273 bit、box(12,25,14)、全次数38、各変数の最小指数0。
出力body SHA `4eef3dad39c0906804508505c9a867a5671f7537a1471ccb272d2b7d71d5defe`。

独立 ROOT closed126-composition multinomial検証
`S61_NONSQUARE_PREPOWER_Q1_B4_ROOT_MULTINOMIAL_REVIEW_20261105_v63`
SHA `4fad8d57f18c57173386a8dd38ae27bc93cf9034608d38f87bb35b4aae4e2806`
でも元4811全係数とD197の一致を確認した。これはlex除算/native factor engineを使っていない。

既知FULLU単元B⁴の除去なので十分開の意味は変わらない。
新しい因子・GCD、二tailの互いに素やunit ideal、実bad Fitting等号、補集合消滅、全FULLU57/239は主張しない。
