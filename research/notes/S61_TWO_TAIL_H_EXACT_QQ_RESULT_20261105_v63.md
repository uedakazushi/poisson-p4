# 確定Hの直接QQ分解：全係数掛け戻しの一次結果

2026-11-05。ROOTのfresh全source/code監査
`S61_TWO_TAIL_H_EXACT_QQ_ROOT_CODE_REVIEW_20261105_v63/report.json`、
SHA `b9be184d3a7bfab6490d1bc1f996794ec51c1dbe0601e319fe7bd86297fc51cd` 後、
一回だけ直接三変数fmpz_mpoly_factorと厳密整数掛け戻しを実行した。
一次報告 `S61_TWO_TAIL_H_EXACT_QQ_ENTRY_20261105_v63/report.json`、
SHA `835f6a6404958bf947efa5e1809f49bbe95872ce2033fd89bdab0d574c640a5b`。

primitiveで正のlex先頭係数を持つ二つの非定数候補を得て、元Hへ全係数で掛け戻した。
ここで局所的にR7,S4と名付ける。rawH rank7や旧Ka13 monicH4とは別の多項式である。

```
H=R7*S4,
R7:29terms, totaldegree7, (adeg,edeg,ddeg)=(4,6,2), max28bits,
S4:11terms, totaldegree4, (adeg,edeg,ddeg)=(2,3,2), max12bits.
```

小さい方は

```
S4=125a²d²−900aed+100ad−4000a−108e³d+36e²d
   −540e²+21ed+120e+2d+160.
```

R7の全29項とS4は `canonical-factor-candidates.json`、
SHA `4dea15b612a158ff9fedfac45b493c1ae3c1c0218f6b196842b242577a2bd979`。
全掛け戻し証拠は `canonical-factors-and-full-product-certificate.json`、
SHA `99775d61d9977765baaea45609c751953b23018b94255b9942f0adc75c3e4f98`。
元の94非零係数だけでなく350boxの他の零係数も一致した。

verification348積、観測最大40 bit、Python owned12834、native accessible12454、
約0.045秒、入力込みRAW333783 bytes。384 bitのstage gateは全て通り、停止・retryはない。
native内部の高さ/supportはopaqueのまま1.9GB/time制限で扱い、内部384 bit認証とは言わない。
外部canonical候補357 bitと、掛け戻し384 bitの硬い上限を守った。
初回のmissing-header compile-only失敗の全bytesも保存した。

二因子の全積は一次結果であり、ROOTの独立backproductを待つ。
baseのQQ既約性については完全FLINT factor APIの契約も使うため、独立scopeを別に確認する。
この分解を独立受理すればHはQQで可約であり、絶対既約でもない。
各baseが絶対既約かどうか、geometric componentが正確に二つかどうかはまだ結論しない。
実bad Fitting、actual full57/239、補集合消滅について新しい主張はない。
