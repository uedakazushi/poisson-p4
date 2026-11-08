# 被約非concurrent m三の残二全源：別作者Fraction監査

2026-11-08、PRIMARY E、`/root/generic_exhaustion`。
`PLANE_REDUCED_M3_WHOLE_SHEAR_FACTOR_EXCLUSION.md`（d08312cb…）を独立に読解し、元の二つの全必要源M/SF_C3_UNUSEDが既約被約actual quinticを持たないという全pointwise結論を受理する。新first symbol、M1、profile、点、jet、GBを作っていない。

## 全源と実変換の独立検証

旧全九源のcomplete constrained K1及び別作者完全性receiptを入力にする。今回producerのcode又はそのarithmetic moduleをimport/実行していない。異なる旧独立q4 Fraction engineの加算・乗算・微分・有理行列inverseだけを再利用し、parser、五変数source30/9、raw36/Euler6射影、有限push、Jacobi、signed Pfaffianは別に実装した。CAS、有限体、補間、数値rankは使わない。

全六実shear a_i→a_i+s_i u+t_i vを、六独立parameterのままfull tensorへ作用させた。**P0+全residual P1+全九P2**に対して、逆shearを掛けると元の全tensorへliteralに戻ることを検査した。raw36/Euler6の定数逆でweight一をcanonical30へ戻し、その全係数が保存complete K1に対する全7/8次元定数同型に一致することを確認した。保存した同型と逆は左右で積がidentity、そのdetは非零定数である。従って全kernel係数を六shearと一／二residualで全域逆変換できる。

さらにcanonical regauge後の全九P2係数を保存した。各係数は元の対応q_iに、shear/residualだけの多項式を加えた形であり、全九変数の線形部はidentityである。従って高次項を「P2に入るはず」と仮定せず、全九方向を持つ全域triangular同型として実検査した。Euler regaugeは全actual Dを変えないこと、finite shearのdet一に対するfull Euler-Pfaffian引戻し式も全parameterで検査した。

この実同型から、任意の元K1/P2係数を持つPoisson点は表示normalized sourceへ移る。leading amplitude零、shear零、homogeneous P2の零／非零を全て含み、unit localizationは無い。

## 元全350式と実因子

同じnormalized源にまず全九P2を独立未知量として加え、元全350式を定義から展開した。新しいM1は計算していない。既constant J2 minorとその逆を検証し、全inhomogeneous解と全homogeneous自由度を保持した。

Mでは九P2の解は全て零。代入後の全350式は全rで恒等零であり、signed Euler-Pfaffianは

```
D=(r/2)(bu+cv)(3au²+6buv+3cv²+r v³).
```

r零ならD零、r非零なら全点で正次数二因子を持つ。

SF_C3_UNUSEDでは八P2が零、Pab=q uvの全homogeneous qだけが自由に残る。代入後の全350非零式は−2rs,2rs,4rsなので、正確なcompatibilityはrs=0である。qは任意であり、compatibility前から

```
D=−c[2as u³−2br v³+(q+rs)u²v²]
```

が成り立つ。括弧が零ならD零、それ以外でも固定c因子がある。generic squarefree leadingのrs非零を使わず、この全proper factorだけで既約quintic枝を空にできる。

いずれも実PGL変換／Euler regaugeの逆を通して元の全源へ結論を戻せる。必要linear kernelの一chosen sliceだけの空性ではない。

## 受理範囲と実行receipt

新focused checker `classification_reduced_m3_whole_source_separate_review.py` はgate PASS後に一回実行し、0.315秒／25,340KiBで全PASSとなった。全six-parameter push、full9 triangular translations、constant inverse、全350／全Pfをreportへ保存した。旧producer/source/job/certificateは一切変更していない。

全root-coverの幾何的必要尽くしは、既独立 `PLANE_REDUCED_NONCONCURRENT_M3_ROOT_COVERS_SEPARATE_REVIEW.md` の全射程を読解して用いる。旧全九算術receiptの他七源は既proper factorで空であり、今回残二源をpointwiseに空とした。これらを接続すると、reduced nonconcurrent R・rank二plane・normal位数三の対応既約quintic全枝が空となる。concurrent、nonreduced R、他のnormal位数を本receiptから新たに閉じない。
