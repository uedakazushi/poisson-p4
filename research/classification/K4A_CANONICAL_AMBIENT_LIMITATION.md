# K4a 全一軌道：正確な orbit23・ambient tangent24 と停止状態

2026-11-08、PRIMARY G/E。許可された同じ canonical point の一回の全ambient判定は完了。
producer 1.364秒、別作者alternate Fraction 0.377秒で同じ全行列/核/小行列を再構成し PASS。
**成分性は証明されず、新成分は零。** 真の法線一つを得たことだけで新成分を数えない。

入力の全source・同じ原始整数18Pは `K4A_CANONICAL_FULL_AMBIENT_CHECK.md` を参照する。
元350 Jacobi、全150係数微分D、Euler25/GL5 25/scalar1の全51列Mを残して

```
rank_QQ D=101, rank_QQ Euler=24, rank_QQ M=48,
D*M=0, dim kerD=49.
```

全疎行列、独立行/列、非零QQ determinant、自由座標が単位行列になる全49核basisを producer に保存した。
別作者も全350×150と150×51を新たに作り、全Dkernel/minors/actualPfaffianを確認した。
従って正確に projective orbit23、ambient tangent24。
全K4aは独立受理された一actual orbitへ入るため、全source像/閉包の次元は23である。
source閉包が全ambient成分であるとはまだ言えない。

actual quintic は18Pについて `-324*F`。F は primitive c一次の既約非正規式、
被約 singular support は `u=a=0` と `u=v=0` の二面である。
これは印刷されたoverall-scalar未定のFをactualPfaffianと同一視せず、全係数で決めた係数である。

## 同じ保存済み線形代数から読める唯一の normal character

stabilizer torus の射影に有効な重みを `(c,a,b,u,v)=(2,1,1,0,0)`、tensor scale=s^-1 と取る。
二次係数 entry `(i,j;e)` の character は `wi+wj-Σek*wk-1`。
保存済み rref kernel の各列は同じcharacterを持つ monomials だけからなる。
独立M列は Euler `wi-wj-1`、GL action `wi-wj`、scalar零であり、保存された全48 pivot列で数えると

| character | -3 | -2 | -1 | 0 | 1 | 2 |
|---|---:|---:|---:|---:|---:|---:|
| kerD | 2 | 8 | 15 | 14 | 8 | 2 |
| actual M | 2 | 8 | 14 | 14 | 8 | 2 |

よって quotient kerD/M は一dim、characterは正確に -1。
これは新しいJacobi/normal/jetの計算ではなく、同じ全行列の独立列と既torusからの読解である。
character だけで actual Poisson pencil 又は局所支持を推論しない。

高次normal、別point、追加minor、miniversal延長は実行していない。
次は同じ保存済みcocycleが literal compatible Poisson tensor になるかを手で認識するか、
既知24親の actual source 構成へ戻る。
canonical orbitが23だから既知Ga9/nil22long23との相違の手不変量は有効だが、
他の24親又は未同定の23/24成分の可能性は残る。
