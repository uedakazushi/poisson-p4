# Actual global P22：独立operator・source輸送監査

2026-11-04。FULLU上の frozen `S61_NONSQUARE_KERNEL22_EXTRACTOR_20261104_v63`
を監査する。作者producer・20pivot・66minor・Bézout探索を再実行しない。
実K174×22は既受理の元六quadに束縛済みのものを使う。

## 独立なexact検査

scalar2491/linear1685/atoms378の全nodeを、別のreview checkerで読む。
scalar nodeはcoreのみ、linear nodeは174入力とscalar倍の和だけであることを
topologyと型で確認する。全atom分母をa,d,F,B,C,Jの既知因子で正確に除き、
非零QQ定数だけが残ることを確認する。

作者の22列を同時に保持するevaluatorとは別に、**Kの一列ずつ**作用素を評価し、
全22×22=484成分をQQ(a,e,d)でI22と比較する。演算utilityは既受理のFLINT
rational headerを依存として明示する。したがって実R線形写像P22について
`P22 K=I22` が全FULLUで成立する。generic点や標本点の検査ではない。

さらに、元六quadと保存L,ellを別のFraction疎辞書でcert light座標へ作り直す。
その係数を新しいnative inputへ渡し、以下を追加で厳密に検査する。

* pair archiveの元memberと公開copyが同一。旧pair第2行は実C*q5、pivot5a、
  Schur factor−aC。全12実係数、元5係数と116非零項cofactorのBézout等号、
  保存pair dualの全4 pairingを確認する。
  元cofactor JSONには零placeholderが一つあり、保存record数117と非零項数116を区別する。
* ell=0を実際の線形代入から行い、保存した四quadの全84cleared係数と比較。
  old inverseのPをadへ項別に戻し、shifts(0,2,2,1)、
  m=(1,J²/25,J²/25,J/5)、g=−a⁴d²Jによる輸送を適用した全16 pairingを確認する。
* L_x1=−40C、ell_x1=−16CJ/243、g0=−80ad²/81とD0の作用を確認。
  pair dualからlambdaの全15係数を独立に作り、四bar quadの消失とellbar L→1を確認。
  global六quad dualについても元の同じ六quadへの全36 pairing=I6を確認する。

以上は実actual sourceの有限なquad/dual検査であり、239列の四次行列の展開ではない。
新22係数はP22の最後のoutputそのもので、全K列上で旧21列を殺し最後を1へ送る。

cofactor変換はq′=Aqに対しc′=A^(-t)cである。JSONのinverse tableはA^-1の
**rows**なので転置して使う。Koszulの戻しは曖昧な∧²表記で済ませず
`Sorig=A^t Sprime A`、すなわちminor A[rows(i,j),cols(k,l)]で確認する。
det A=−80は既受理の単元変換。全484の比較は最終basis順とこのorientationも検査する。

## 大域282表示と限界

受理済みM304=[A0,B0,C0]の実indexを、quad174・残るsourceの最初22・残108に
一致させる。`(k,z,h)↦(Kk+z,Pz,h)` の逆は
`k=Px−u, z=x−Kk, h=h`。PK=IとA0K=0により両合成は恒等である。
したがって同じcokernelを持つ **実global free282表示** は
`[A0+B0 P22,C0]` となる。これは係数oracleとして完全な有限recipeである。
旧282行を選んだ固定submatrixでも、全239係数の展開済み行列でもない。

P22は既知22方向をsplitする進展であり、新しいkernel方向を増やさない。
既知独立方向37・残る43という状態は変わらない。全239全射、C122 rank64、
欠損curve/pointの消失、affine exhaustionを結論しない。

元のgzip line-wrap reader failureとnew22 checkpoint/resumeはhashで保持する。
本レビューは旧source・prefix・operatorを新しい入力へ束縛し直したもので、
作者のnative jobや失敗jobの再実行ではない。review executableと非圧縮入力は
ignored workへ置き、ソース・圧縮入力・log・receiptだけを保存する。

review側の全native検査と282 index照合が通過した後、最後のreceipt作成だけで
relative work pathのmetadata例外が生じた。成功log・入力・当時source hashを
`finalization-frontier.json`へ保存し、pathをabsoluteに直したmetadata-only入口で
公開を完了する。native算術は再実行しない。全体のwalltimeを後から推定せず、
入力保存から成功logまでの実filesystem時間差と、publisher時間を区別して記録する。
