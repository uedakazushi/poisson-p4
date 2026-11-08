# 定数modular三完全源：固定三witnessの全ambient診断と停止範囲

2026-11-08。root承認の固定一batch、各whole source一つのみ。
`CONSTANT_MODULAR_THREE_SOURCES_EXACT_DIMENSION.md` の三閉包は各既約24次元だが、
今回の全ambient診断は三つともcomponent判定を閉じなかった。新成分増分は零。
全既知親への非所属、全sourceの接次元一定又は成分性は推論しない。

入力は `inputs/constant_modular_three_sources_fixed_witnesses_20261108.json`。
κ=1、全10/8/10係数をliteral integerで固定し、元source式から全係数を戻した。
normal formだけのpointをsourceのactual点と仮定せず、全parameterの読み戻しを元150座標上の
signed unit rowから確認した。full parametric source自体のdivY=1,Yf=fも再確認した。

各点をh-Casimir quadraticへ同次化し、全150係数、全元350 Jacobi、全126 actual Euler-Pfaffianを確認。
全actualDは順にh³zw、h³(zw+hu)、h³(z²+hu)と厳密に一致した。
新しい制限slice／normalizer条件をambient行列へ加えていない。

| whole sourceの固定一witness | rank D350×150 | Euler | 実GL＋scale＋Euler | 全source微分も追加 | projective tangent |
|---|---:|---:|---:|---:|---:|
| I rank二critical零 | 98 | 24 | 47 | 49 | 27 |
| II rank二noncritical | 98 | 24 | 48 | 49 | 27 |
| III rank一noncritical | 96 | 24 | 47 | 49 | 29 |

全D×(group/source directions)は零、全行列とQQ非零最大minorを保存した。
source微分を加えた49は既pure image dimension24と整合するが、full Dの上界代用ではない。
数学的に必要だったD101はどの点にも得られず、componentpositiveは零件。
選択pointを追加したり、別period／法線／jet／minor profileを探索していない。

実行は **6.023993秒／RSS62348KiB**、180秒／768MiBの契約内で全三点が完了した。
新caller全文とSHAはfinal reportへ埋め込み、同じ入力.pyも不変保存した。
既helperのbytesは変更せず、開始／終了の両方で入力hashを検査した。
各caseはfull matrixを個別checkpointへ保存し、final reportはそのhashをまとめる。

## 判定条件とpriorの正確な役割

D101が得られた場合だけ、既pure local dimension24下界とtangent24上界からsmooth unique24 componentが従う。
その後のprior同dimension監査では、既MPS24八件とCPLUS_EXT／dependent3H+Q24／Ga9cuspidal24は
generic actualDがreducedなので、全Dが閉multiplication像{h³Q}にある三Zとは一致しない。
残るknown double-rank3-quadric24のgenericD=h*Q²にはtriple linear factorが無く、同じ閉条件で別である。
三Z同士もresidual Qのgeneric rank2/4/3で異なる。
他dimensionのpriorはambient componenthoodが先に得られたときだけdimensionで別になる。

**今回はその前提D101が得られていないため、このprior比較を未知の大きい親への非所属とは呼ばない。**
CA/CB28、NEW25、M11など25成分やMPS大成分への真包含は、各actual factor条件が許す限り未決である。
rootのprior133監査／新countへ本batchをpositive resultとして渡さない。

契約通り、次手はwhole sourceのactual known-parent confluenceである。
sourceIIには別純証明 `CONSTANT_MODULAR_SOURCE_TWO_NEW25_WHOLE_COVERAGE.md` で全NEW25 bridgeを得たが、
その帰属は本ambient rankからの推測ではなく、独立の全actual map／normalizer／閉性による。
I/IIIの親はこのnoteでは未解決。
