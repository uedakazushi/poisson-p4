# T0/C⁴ 商の独立 ROOT 整数長除算検証

2026-11-05。一次計算の有限和 engine とは別の固定25項 C⁴ 整数長除算が通った。

独立報告は
`reports/S61_NONSQUARE_PREPOWER_T0_C4_ROOT_LONGDIVISION_REVIEW_20261105_v63/report.json`、
SHA `4d17fc39925b4c35e3b636c368c16083b582da5bd876201b66ff8cc9476feb72`。
元の `root239.poly` の7208 QQ body から整数分母を取り直し、外部 offset4を
二重に適用せず、先頭 a⁴ 係数100⁴で通常の整数長除算を行った。
finite-series、R-power、100²⁰ scaling を検証 engine に使っていない。

全4500商係数、全余り零、元の D197との同一性を確認した。
実際112500 products、最大320 bit、30690 owned slots、約0.152秒。
これにより一次結果の

\[
T_0=a^4C^4J_0/D_0
\]

という全係数恒等式を別原理で受理する。商は元の D197を保ち、内容や GCD の探索はない。
原 primary、SOURCE、worker custody とその値未実行という準備時点の文言は変更しない。
他因子・tail や、実 bad Fitting ideal、補集合の消滅、全 FULLU57/239は主張しない。
