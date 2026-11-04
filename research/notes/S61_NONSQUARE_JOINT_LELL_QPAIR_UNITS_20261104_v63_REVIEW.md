# qpair の固定14 unit pivots と44×122 frontier の独立監査

2026-11-04。対象は `S61_NONSQUARE_JOINT_LELL_QPAIR_UNITS_20261104_v63`。
入力58×136は、直前の overlap independent REVIEWで full q6・unit98・overlap5から全係数を再構成した
保存入力そのものである。今回もsource hashを照合し、新しいGaussian、candidate選択、core点、標数、GBは使わない。

## 固定操作だけの再生

独立な QQ Laurent dictionary arithmeticで保存済み14 pivotを順番通り再生する。
全pivotは literal

    19683/(1280*a²*d⁴)

である。未知因子を反転しない。
各stepで実際のpivotと保存値・既知unit factorizationを照合し、行全136係数を正規化する。
全40 elimination係数を **操作前**の行列と照合する。
authorのactive-column最適化を使わず、target行の全136列を毎回更新し、以前のpivot列が零のままであることも確認する。

最終的に checkpointの **全58×136=7888係数**を比較し、残る44行・122列の順序と全値を照合する。
残余は1223非零 entry、79165 numerator terms、最大391 terms/entry。
未完成stepはないこともreport/checkpointから確認する。

## selected minor と全域下界117

選択順に並べた元14行・14列のminorでは、後に選ぶ行は、それまで選んだ行の組合せを引くだけである。
保存操作により、このminorは上三角にされ、正規化を戻したdiagonalが14 pivotの列になる。
従ってdeterminantはpivotの積で、保存constantと a^{-28}d^{-56} に一致する。
選択順をそのまま使うため余分な行列式signはない。

このminorはFULLUの既知単元なので、前のfour-quad image103方向に加えて全点で14方向が独立。
**六quadの degree4 imageの全域下界は117**、target239のcomplementは122。
rank152の全点下界や残る35方向は証明していない。

残る1223非零係数について、a,dのLaurent monomialを外し、指定六因子 a,d,F,B,C,Jだけによる
厳密なpolynomial divisionで全てを調べる。全てに非constantの因子が残り、保存した
「個別の既知unit entryがない」という停止条件と一致する。
これはrow/columnの線形組合せでunitを作れないことや、cokernelが非零であることを意味しない。

## full presentationへの安全な帰結

保存された同じunit操作をfull lossless221→141へ作用させれば、overlap5とqpair14の
合計19 source/target pivotを消して、同じcokernelの **202→122** 表示を得る有限recipeがある。
ここでは13個の零overlap source行を保持して数えている。
112 genuine cubic倍行と33quarticは同じtarget projectionをlazyに受け、削除しない。
この大きな係数行列を本監査で展開したとは主張しない。

既知kernel38方向の数は変わらず、目標80まで残り42。
新しいtarget pivotに使うsource座標は、古いfree183表示のarbitrary B22/C108と混ざることがある。
P22とrestorationの合成がleft inverseを与えるだけでは、残るprojective quotientのglobal free basisは完成しない。
従って **free164はこのpacketから自動的には主張しない**。そのsource basis構成は別の証明書を要する。
全122/239全射、全kernel exhaustion、全点cubic57も未証明のまま保持する。

```sh
. /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_joint_Lell_qpair_units_review_20261104_v63.py --output work/marked-square-v63/joint-qpair-review-REPLAY
```
