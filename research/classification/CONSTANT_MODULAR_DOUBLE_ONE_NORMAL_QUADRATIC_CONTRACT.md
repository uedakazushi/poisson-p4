# constant DOUBLE：保存全周囲行列を使う一normal二次支持の有限契約

2026-11-08、PRIMARY E_multiple_hyperplane。**準備のみ、root GO前、未実行。** 既存producer/report/helperを変更しない。新しい点、normal profile又はoperatorを選ばない。

## 分類への新しい必要依存

最高double-line全源は独立有限還元により、既知閉和K又は一つのforced actual orbitへ戻る。そのcanonicalは(t,x,y,r)順に

```
(1+x-2xy,4y²-y,(1-6y)r+x,0,-x²,1+2xy),
D=h²*(x²*y+x*h²+h³), N_h=-h*partial_t.
```

既存ONE fullambientはD350×150 rank101、Euler24、actual Euler/GL/scalar M150×51 rank48、kerD dimension49を保存した。射影actual orbit23、全接dimension24。固定hのT37中心二法と単独nilflowはconstant cofactorを作れず、同じ曲線の係数・次数変更反復は停止した。rootは、この未処理orbitの成分支持へ直接効く**一つの二次判定**を新PRIMARY依存として限定的に認めた。これは旧local jet系列の自動再開ではない。

## 真の切片とcoordinate inverse

保存したMの48独立columnをB、対応48row coordinate projectionをAとする。ABは保存nonzero minorで可逆。保存完全kerDの49basisから、Bに入らない最初の一column kを、

`R=k-B*(AB)^-1*A*k`

へ射影する。A R=0、D R=0、R≠0となり、残るnormalはkerD/imBの唯一の一次元である。別normalを試して障害を改善する手順は存在しない。

保存Dの101pivot columnのunit vectors e_jを

`U_j=e_j-B*(AB)^-1*A*e_j`

とする。AU=0、DUの保存101row minorが可逆。F=[B,R,U]は150×150の真正な全直和である。両inverseを、DUの可逆minorとAB及びRのunit座標から明示し、F*C=C*F=I150をQQで掛け戻す。従って48 coordinate cutの上には正確に一normal+101非核方向だけが残る。

## 二次余核の最弱criterion

全350 projective Jacobi二次mapをJとする。同じRのQ=J(R)を一度だけ計算する。保存Dの像membershipを既101minorで判定し、Q-DU*tの全350残差を見る。

非零なら、その最初の非零残差rowから全D150columnを零にする明示QQ functional lambdaを作る。全lambda*D=0とlambda*J(R)≠0を確認する。101選択式を形式陰関数で消去した後の一変数残式は、非零二次先頭係数を持つ。よってその根基は原点、全raw局所次元は48以下であり、48次元actual orbit/Euler/scale像と一致する。Euler24とprojective scale一を引き、**forced orbit閉包は23次元の全reduced component**となる。

これはone-variable Krull支持命題だけで、完全IFT graph、smoothness、正確nil指数や全component相違を主張しない。新成分数は全prior auditと別作者受理後に統合者だけが決める。

残差が全零なら、完全なQQ補正tと全350掛戻しを保存し、**二次余核零・成分支持未決で停止**する。cubic、別normal、別点、追加source、微小parameter変更は実行しない。

## 全方程式と再利用の射程

保存した同一raw150、fullD、Euler/GL/scale全列、完全kernel、全nonzero minorと列／行順を読んで再利用する。定義どおりのprojective Jacobiは全350係数であり、特定134式の一部への置換ではない。受理済みS61全134とのglobal polynomial coordinate同値は `S61_GLOBAL_COMPONENT_EVALUATION_FIBER.md` §3と別作者査読へ接続する。今回は全350の方を直接使うため、normalization chart又はbeta分母を新たに反転しない。

新しいD/operator/rank profileは生成しない。必要な旧matrix readbackは、同一canonicalにおける全linearized coefficient identityの確認であり、別点や新rank改善ではない。raw150と全Jacobi/Pf/cofactor、actualM51を独立sparse arithmeticで照合する。101/48のrankは保存完全kernelとnonzero minorsで確認する。旧producer又はcurrent sharedrunnerはimportしない。既独立polynomial/linear algebra helperだけを使い、新inline全bytesを結果へ埋め込む。

予算は120秒/1536MiB。fresh gateと明示root GOが前提。正負とも同じ一つのreportへ全source/hash、F/両inverse、R/Q、補正／functionalを保存し、失敗時も上書きせず失敗bytesを残す。
