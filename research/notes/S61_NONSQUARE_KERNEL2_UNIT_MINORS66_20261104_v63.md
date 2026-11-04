# 残余154×2の指定66 minor: limited negative

2026-11-04。frozen partial20 Schur matrixを読み、前packetが列挙した66 pairsだけを
一回のbounded jobで検査した。20 Gaussian pivots、元22恒等式、fullquartic matrixを
再計算していない。候補の自動拡張・GB・Bézout・unknown inversionは行わない。

## 実際の結果

候補は、両列activeな最小term12行

`154,125,126,140,169,155,143,139,144,124,135,141`

の全66組である（行番号は元source174）。各entryをn_ij/d_ijと書き、
全denominatorが `a,d,F,B,C,J` と非零QQだけで割り切れることをexactに確認した。
行列式を、未知gcdやfactorizationを使わず

```
N=n00*n11*d01*d10−n01*n10*d00*d11,
D=d00*d11*d01*d10
```

で計算した。Dは既知unit。N/Dをunitと認めるのは、Nから同じ既知factorだけを
割り切って非零QQが残る場合に限った。

**6組は零、60組は非零で、unitと認められる組は無かった。**
全66のliteral N/D、既知factorのnum/den powers、残ったpolynomialを保存した。
非零Nのdegreeは9..15、既知factor除去後のremaining polynomialは非定数であり、
各polynomialのtermsは最大54。これらをさらにfactorizeしていない。

従って、この66組からglobal P2/P22または固定282 submatrixは得ていない。
前packetのexact P20、固定actualM284、残余154×2のFULLU全点rank2はそのまま保たれる。

## 限定と次のgate

全active row pairsは5565あり、今回の66はそのsubsetである。残るpairsの不存在、
global unit minorの不存在、Bézout splitの不存在を主張しない。
特にindividual entriesにunitが無いことはunitな2×2 determinantを排除せず、
今回66 determinantsにunitが無いことも、複数minorsのunit combinationを排除しない。

次に必要なのは、別途承認した小さい2-minor検査か、既存のwhole22独立性に使った
unit anchorsを**実際のSchur coordinateへ輸送した**Bézout constructionである。
`ell_x1=−16CJ/243` と
`3*Delta_x2−(7/2)*Delta_x4=−80ad²/81` を、そのまま残余entryと仮定してはいけない。
有効な証明には、残余154×2に対するliteral2左逆またはminorsのexact unit combinationが
必要になる。fullquartic239、C122rank64、残る43 kernel gateは未解決である。

## 保存と実行履歴

最初の実行は0.243秒、最初のdeterminant計算後のdegree metadata取得で
PolyElementに存在しない`total_degree()`を呼んで停止した。最初のoutcomeは
書込み前であり、そこからunit結果を主張しない。失敗source/hash/frontierを保存した。
metadataを `max(sum(exponents))` に直した一回のrepair runは0.560秒で全66件を完了。
合計約0.804秒、元の60秒/2GB/32MiB枠内である。

`code/s61_nonsquare_kernel2_unit_minors66_20261104_v63.py` と公開outcomes/reportが
実行済み証拠である。unitが見つかった場合のP22 circuit組立branchは今回実行されて
いないため、そのbranchが検証済みという意味ではない。今回の保存結果はlimited
negativeとpartial20 frontierの維持に限る。
