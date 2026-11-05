# 固定四次商環でのD12余り：非零の一次結果

2026-11-05。ROOTの実行前監査後、一回だけ指定の固定軸商環workerを実行した。
一次報告 `S61_TWO_TAIL_H_D12_AXIS_QUOTIENT_ENTRY_20261105_v63/report.json`、
SHA `23646500355809a05d19af1774fa33ce75a8ec56f0e31151f52fa56ed76063c8`。

Fp32003上の四次商環のmonic modulusは低次順に

```
h=[17430,11013,30453,27809,1].
```

既存独立sourceで固定された整数C12の全144成分をこの環に移し、
全4096subsetのsigned determinant DPを完了した。余りは低次順に

```
detC12 mod h =[4884,25267,23474,25095],
D12_axis mod h=[27788,575,7510,2387].
```

後者はa5倍とLCD逆2378で正規化したもの。LCD residueは30563。
a逆 `[754,17498,1840,29225]` の積1も同じ可換環で確認した。
全subset値、全ring入力、a冪、modulus、cursorは
`complete-subset-and-ring-certificate.json` に保存し、
SHAは `f52f736ce74e17ab8ceb4fc0421a609b80a89ddfe6f5ca66ed202f2761a1c227`。

691393 scalar積、所有26996、field一時値30 bit、約0.784秒、入力込みRAW135478 bytes。
guard失敗・再実行・別axis・新prime・旧source再評価はない。

この非零余りは、sourceに記したGauss/合法制限/単元性によりHがD12を割らない一次証拠である。
独立値peer前には非除算の最終受理としない。
独立受理後はsquarefree Hとunit-coprimeによりV(H)∩FULLU∩D(D12)の非空性を導けるが、
それでも選択二tail十分条件のgapに過ぎず、実bad Fittingやfull57のrank低下を意味しない。
全FULLU57/239、補集合消滅、尽くしについての新たな主張はない。
