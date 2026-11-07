# Hの厳密二分解とS4絶対既約性：独立受理

2026-11-05。ROOTの独立報告
`S61_TWO_TAIL_H_EXACT_QQ_ROOT_DENSE_REVIEW_20261105_v63/report.json`、
SHA `b7a6b60c04a288e54128f0df3df8a601152e8d58f576f9fbc0d6b0eb94fc8625` を受理した。
APIやproducer packed readerを使わず、二factorの順序を逆にして三変数dense tensorへ整数積を足し、
元H94と全350係数が一致した。整数330積、最大40 bit。従ってH=R7*S4は独立に確定した。

S4はaの二次式で、a先頭係数125d²、a一次係数はd=0で−4000。
従って代数閉包上でもa-contentはunit。a-discriminantはe次数3、先頭係数54000d³なので
C(d)(e)の平方ではない。GaussによりS4はC[a,e,d]で既約である。
この分解とS4の絶対既約性にFLINT APIの完全性契約は使っていない。

R7のQQ既約性・絶対既約性はこの独立報告では未証明である。
Hは既に非自明なQQ分解を持つため絶対既約ではないが、geometric componentが正確に二つかは
R7の絶対既約性を別に証明するまで結論しない。
actual bad Fitting/full57/239/尽くしの主張は含めず、一次packetも書き換えない。
