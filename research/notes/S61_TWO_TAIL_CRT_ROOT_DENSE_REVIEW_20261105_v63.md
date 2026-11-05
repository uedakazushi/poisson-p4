# 固定二tailの最大公約因子：独立整数検証の受理

2026-11-05。一次CRT結果を変更せず、ROOTの別演算による検証をここに記録する。
独立報告は `S61_TWO_TAIL_CRT_ROOT_DENSE_REVIEW_20261105_v63/report.json`、
SHA `6d7f852148000830f26d4f9ed709a7e7a5bc205e8c2678e4613648ea82b23bc7`。
検証コードSHAは `25f32c74659add2dd559d25cc9b9db61602ef24e5bec0ada4bb442fd67dbba8d`。

元の整数分子をN0,N2とすると、保存されたprimitiveな94項Hについて

\[
N_0=H A_0,\qquad N_2=H A_2
\]

を独立に確認した。producerの除算・CRT・GCD・binary readerを使わず、Hとcofactorの
疎な順序を逆にし、三変数のdense tensorへ整数積を加えた。
元4500項に対応する全10659係数、元2317項に対応する全5070係数が、それぞれ零係数も含め一致した。
計236786積、実最大293 bit、owned56966、約0.23秒。これは二恒等式の全係数検証である。

さらに既存の独立報告
`S61_TWO_TAIL_FIXED_KRONECKER_QUOTIENT_REVIEW_V2_20261105_v63`、
SHA `9041a2b792e537070184a5c2b3ef33c9c35b96401a3e9b81199de2286b83652f`
を用いて最大性も受理した。固定写像は `a=t^646,e=t^19,d=t`、既存primeは32003。
Hの先頭係数はこのprimeで非零であり、そのmonic像は既受理の次数3899のgcdに全係数一致する。
二cofactor像の先頭係数も残り、次数は6645と4076。
既存の独立証明でこの二商のmodular gcdは1である。

仮にA0,A2に非定数のprimitive共通因子があれば、Gaussの補題により整数商を取れる。
そのe,d次数は共通に高々16,10なので、同じKronecker写像は支持を混同せず非定数性を保つ。
両cofactor像の先頭係数が残ることから、共通因子像の先頭係数も残る。
したがってmod32003の二商に非定数共通因子が生じ、既受理証明に矛盾する。
よって `gcd_QQ[a,e,d](N0,N2)=H`、ただし非零有理数倍の違いを許す。
新たな素数・GCD・因数分解は使っていない。

Hは94項、変数次数(6,9,4)、全次数11、primitive、正のlex先頭係数を持つ。
Hの既約性はここでは証明していない。実際のbad Fitting idealとの等号、全FULLUの
cubic57/quartic239、補集合消滅、一般成分の尽くしも主張しない。
独立検証前の一次ノート・reportは履歴としてそのまま保持する。
