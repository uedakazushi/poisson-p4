# F-Hasse valuation零：独立整数Horner監査の受理記録

2026-11-05。元のprimary result noteとreportは変更しない。
ROOTの独立report
`reports/S61_NONSQUARE_PREPOWER_Q1_F_ROOT_HORNER_REVIEW_20261105_v63/report.json`
（SHA256 `3d4b5fb36b922aac2714a3fb25327b3c1f884e75e64bd3eba885d409b49bc6b5`）を結ぶ。

Q_C=T1/(a³C⁴)はROOTの独立integer longdivisionで得た4811係数、元の共通denominator D197である。
F=500a+90e²+3e−24、R=24−3e−90e²とする。
producerのbinomial／Hasse／R powers engineをimportせず、inputをa次数の係数I_iへ分け、

```
H=I16;
for i=15,...,0: H=R*H+500^(16-i)*I_i
```

という別の整数Horner算術を使った。得られた全665係数と一個のdenominator
`D*500^16`がproducerのk=0層と完全に一致した。
独立計算の記録は20856 products、実最大415 bits、owned28206、RAW546844 bytes、0.035秒。
追加の点・標数、他tail、GCD、商計算は行っていない。

k=0層は非零であり、**v_F(Q_C)=0を独立に確定した**。
completed-zero-orders.txtは空であり、未計算のk≥1層を零と主張しない。
したがってFを除く商のSOURCE／workerを新しく用意しない。

この結論は既知単元の因子が無いという多項式の性質である。
T=0やD12=0の補集合、FULLU全点cubic57、quartic239、実bad Fitting ideal、
一般五点尽くしへの新しい推論は含めない。
既存の十分開 `FULLU D(D12)∩(D(T0)∪D(T1))` と、その未処理補集合を保持する。
