# generic triple全九十source：六十二whole閉鎖と残四line＋十fixed

2026-11-08。producer exact QQ算術PASS、独立review未了。
全firstconnection有限化はTRIPLE_FINITE_NORMALFORMSとTRIPLE_CANONICAL_REDUCTIONの両noteを使う。
全G又は全m三を閉じたとは主張しない。

`CLASSIFICATION_PLANE_NONREDUCED_M3_TRIPLE_20261108.json` は
十五whole α＋七十五fixedについて、全flatness、全3fβ homogeneous inverse、globalK三未知、
P0全350Jacobi、fullP1三十/P2九、actualD3全二十四係数の正確なline条件を再構成した。
fullP1kernelのbad normal三行が零なら、P0/P2の全u-rowと合わせactualDにu因子があるとして空にした。
全α maximal-minor/kernel分母/P0の例外根を元全matricesへ戻す。
log generic classから除く四η値の全fiberは元正常形の二十一fixedへ戻しており、捨てていない。

全90labelのうち62wholeはglobal lift不整合、actualD3像零、bad三行零、又はprojective像二十四以下。
28labelはfalseのまま保持した。
この28labelにはwhole lineの特殊化で同じfirst symbolになるものが含まれる。
actual constant normal scaleでr=2dx又は−dxをr=dxへ戻せることを使い、
全未閉first symbolsを下記四whole lineと十fixedで覆える。

## 1. R2L四whole line

f=a²b、h=a^e b^(2−e)。p=(dlogf)/3+αdlog(a/b)、q=2p−dlog(f/h)。
全e=0,1,2のr=dx（x=b/a）が三whole line。
e=2のr零が一whole line。

これらのP0は全αでglobalに定義され、c非依存でK零。
三dx型は

`P0=∂c∧[a²∂a−2ab∂b+3αaz∂z+((6α+3e−4)aw−3bz)∂w]`。

r零型はe二で同じ式の−3bzを零にする。
元generic-normalizerで四η値を除いていたが、ここでは同じ明示P0の全αへ戻す。
η−1/零のdx fixedはそのまま、η一のr=d(x²)/x=2dx、η負二のr=x²d(1/x)=−dxは
normal ratioで同じP0へ戻る。
全special fiberの新kernelはwhole familyの元matrixをspecializeする際に必ず復活する。

## 2. R2L六fixed

- exact_B e零のα=−2/3、r零又はdlogxの二型。
  p=dlogb+dx、q=−2dloga+3dlogb+2dx。
- log η零の各eでr=dlogx+dxの三型。
  α=4/3−e、p=q、normal diag係数λ=4−3e、offdiag coefficient−3(a+b)z。
- log η零/e二でr=dlogxの一型。
  α=−2/3、p=q、normal diag−2a、offdiag coefficient−3az。

前二exact_B whole linesは他の全αで既空、r=dlogxのみのscalar源もe零/一では既空。
残六fixedだけを保持し、空だった他special sourceを復活させない。

## 3. R3L四fixed

f=a³、h=a²、x=b/a。

- φ零、r=dx。
- φ=x、r零。
- φ=x、r=x dx。
- φ=x²、r=dx。

各々p=dloga+dφ、q=dloga+2dφ。
P0はc非依存/K零で、∂c∧W0のW0は順に

```
(0,−3a²,0,−3az),
(0,−3a²,−3az,−6aw),
(0,−3a²,−3az,−6aw−3bz),
(0,−3a²,−6bz,−12bw−3az).
```

旧fulloperatorのfalse粗上界は保持し、normalizerだけをG閉鎖とは解釈しない。

## 4. 次の構造経路と実行bytes

残十四は全P0がc非依存の∂c∧W0である。
次はfullP1必要kernel＋FULLP2九をEuler gaugeして全Q/W sourceへ移し、
strictQJacobi/actualD3 line/actual normalizerを用いる。
同じcoarse kernelsをminorや長いjetで反復しない。

90batchは157.923秒/562596KiBで完走し、未完checkpointを最終reportへ統合した。
実import helperはb42fdf58で、fullbytesはpreappend archive内に保存している。
旧df9 helperからは全oldFunctionDef AST不変のcompatible extensionで、
current helperの後続追記を実行時入力と混同しない。
全runner/契約/数学証明書は旧bytesを保存してから現bindingを更新する。

m二、reduced R、全pointwise Eは別残義務。
