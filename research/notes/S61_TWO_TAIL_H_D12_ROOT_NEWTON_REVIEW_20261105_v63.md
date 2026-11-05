# HがD12を割らないことの独立Newton検証

2026-11-05。ROOTの独立報告
`S61_TWO_TAIL_H_D12_ROOT_NEWTON_REVIEW_20261105_v63/report.json`、
SHA `041d9234ccf704286eacd0da6e0f6161d8b268b27dfaa678e2080657f5f32193` を受理した。
コードSHA `d310cae53c79166b87e601acac45b95210219fb47662c87c23d6abbd3ae5f632`。

同じ有限free四次可換環S=Fp32003[a]/monic(H_axis)で、subset determinantを使わず、
12回のmatrix productとtraceによるNewton/Faddeev recursionを行った。
1から12は既存primeで全て単元。全12段階を保存し、最後のCayley–Hamilton行列も零。
producer codeをimportせず、体・etale性も仮定しない。

得られた余りは低次順に

```
detC12=[4884,25267,23474,25095],
D12  =[27788,575,7510,2387] !=0
```

で、一次結果と全係数一致した。584206 field products、owned26282、約0.196秒。
新しいaxis・prime・source oracleは使っていない。

Hは完全反転積U=adFBCJと互いに素であり、元D12はa,B分母だけで表される。
もし六因子環でH|D12なら、UFDで反転因子を消して、a,Bを払った多項式分子もHで割れる。
Gaussと先頭係数保持により既存primeへ還元でき、合法なe=d=0制限でもH_axisで割れる。
aはSで単元、B_axis=-40、LCDも単元なので、独立に得た非零余りと矛盾する。
従ってH∤D12を六因子BASE環で独立に受理する。
実bad Fittingやfull57 rank低下の主張は含めない。
一次結果や以前の条件付きnoteは書き換えず、本受理を別publicationとする。
