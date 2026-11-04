# exact unit98 Schur と既知38方向の輸送：actual free183→141

全域は FULLU。joint frame の full四二次は補正を持つ monic leading term
として使い、四単項式idealとはみなさない。
98本の二次倍行と98 target列の係数 square U は、保存された実係数により
unit upper triangular である。元319から221→141の損失のない Schur を作る。
18 overlap、58 q0/q5倍行、112 cubic倍行、全33 quarticを保持する。

source basis を変えるのは二次174本だけとする。cubic112の multiplier は
旧 v_j のままにして、その多項式だけを新 target 座標へ代入する。
これにより旧15方向の quartic11/cubic4 pivot と新 G20 separator は literalに
保たれる。全 cubic multiplier を新 light に変えるとこの性質を失うため行わない。

残る二次76座標を x とし、選んだ98 source を R、残る76 source を Eとする。
二次の残行の pivot target 部分を N_C とすると restoration は

```
beta_R=-U^{-T}*N_C^T*x,   beta_E=x.
```

これを旧二次 source へ戻し、既存P22を作用させて Pnew を作る。
Knew は旧K22を新二次 sourceへ変換して E に制限する。
純二次の実恒等式 M*K22=0 と unitU から restore(Knew)=変換した旧K22。
従って Pnew*Knew=I22 が成立する。source76から単に22座標を捨てるのではない。

旧15の正確な unit削除後は二次76＋surplus130。既存free281と同じ Bnew22、
Cnew108を使い、separator二座標を Cnew に保持する。
`[A76+Bnew*Pnew,Cnew]` の自由184 source と新kernelを構成し、
G20*v1+G20*v6 の既知単元で一つを削除する。actual free183→141 は旧319と
同じ cokernel を持つ。固定した旧source183行の部分行列とは主張しない。

core polynomial を展開しない scalar/linear DAG で、座標変換と逆変換、
unitU と target NF、二次 source の両方向変換、restoration、Pnew/Knew、
新kernel184とunit削除、全 source/target index を保存する。
新kernel separator は独立監査済みfree281の同じscalar根とliteralに一致する。
全33quarticの係数・183x141全行列・新しい点や rank は計算しない。

新しい実完了条件は、この正確な183→141写像が全域surjectiveであること。
全点rank141/full239・一般βでの尽くしは未解決。既知38方向、残り42も変わらない。
