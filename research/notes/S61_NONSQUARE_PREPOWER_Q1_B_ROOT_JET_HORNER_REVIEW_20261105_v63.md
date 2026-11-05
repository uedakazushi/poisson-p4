# Q1 の B付値4：独立 ROOT 五jet d-Horner 検証

2026-11-05。独立報告
`S61_NONSQUARE_PREPOWER_Q1_B_ROOT_JET_HORNER_REVIEW_20261105_v63`
の SHA は `a5266d75e23f5643a950771e0b3ebac30eb6fc355575ca7925cbba9b87838042`。

元のcanonical4811整数から下降d-Horner五jetを計算し、低次四層全零、第4層1107係数、
元D197とsymbolic b1¹⁸の一致を確認した。一次nativeのS-cacheやbinomial式を使っていない。
859755 products、最大399 bit、39448 owned slots。
これにより普通のQQ[a,e,d]で正確な \(v_B(Q_1)=4\) を受理する。
b1をFULLU単元とはせず、b1=0の成分を捨てない。

一次資料とSOURCEは変更しない。B⁴商はまだ未計算で、新しいworker/source/whole auditとONE GOが必要。
他因子・GCD・二tailの互いに素、補集合消滅や全FULLU57/239は主張しない。
