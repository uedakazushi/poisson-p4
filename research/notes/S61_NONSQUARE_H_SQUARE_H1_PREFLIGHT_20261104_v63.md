# 新global122 targetのH² pivotと完全H¹ image preflight

117二次pivotはH次数を混ぜない。残る122 target は pure98、H*light²23、H²1。
前段の追加19pivotのうち一つはq0*Hなので、H次数1の数は24ではなく23である。

actual free164には実quartic77,131が literal に残る（元319 index303,318）。
いずれも旧15の削除対象でなく、Bnew22はcubicのみ、新separatorのsource変換も
cubic二座標だけに作用する。従って既存の正確なsource変換

```
chi=u*Q77+v*Q131, psi=-b131*Q77+b77*Q131,
u*b77+v*b131=U=(50aBC)²
```

は同じ既知単元determinantを持つ。joint座標でHは変えず、117二次pivotにも
H²列はない。H² unitを消せば自由163→121、同じactual319cokernelを得る。
targetはpure98+H¹23。純98係数をここで展開しない。
chiの全重み4多項式DAGとexact係数抽出関数を保存し、必要な純係数もlazyに
定義する。有限jetの一致を全次数へ読み替えていない。

H¹については既存924個のquartic係数を全て再利用する。
**T34のraw係数をそのまま使わない。** Tn/Qhatへ移すにはH*x_i*x_j係数を
`D_i*D_j/D_H`倍する。H²係数は不変。その後jointlight逆座標とNF117を作用させる。
genuine cubic16のH-linear112係数は元々Tnの凍結済み57sourceから取り、追加の
adapter因子を掛けない。

H²消去後、全imageのH¹射影は、残るq5*H一行、全112cubic倍行、全21quartic
（psiと残る20行）の134x23冗長presentationで正確に与えられる。
旧15 unit削除とP57 complementは全imageを保存し、除いたcubic4も元の正確な
関係から復元できる。残るquadのpure行はH¹を持たない。
従って112行だけの十分条件へ落としているのではない。

全H¹係数はscalar DAGとして保存し、構造的support上界と実行費用だけを調べる。
非零DAG根を関数の非消滅とみなさず、H¹ rank・Gaussian・pure98行列展開は行わない。
全点rank121/full239・pure98の実像は未解決。既知38方向・残り42は不変。
