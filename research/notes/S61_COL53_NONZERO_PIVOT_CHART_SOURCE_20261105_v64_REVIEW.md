# col53 非零成分による追加ピボット chart：独立 SOURCE 検証

独立検証は PASS した。報告 SHA256 は
`8911d7071992759ca5544463fa6d098f97887e0bbf50dcb8abe169a3a016059e`。
実二面の独立符号付き source replay `d1708919…` と、係数側 owner・
現行 full167×14 の凍結 manifestを固定した。新しい行列・RHSの値は求めていない。

指数 `(-3,0,6,-1,10)` の係数 `241815/134217728` は、
全1,203 SSA node・root1202=`new_1013_sum` の有限 Laurent式 F11 の係数である。
411項の H=−3 面を F11 全体と取り替えない。
非零のこの係数により、F11 は同じ parameter関数体で非零である。

独立な十二記号の整数多項式として、元第九段の二つの式を代入した
反対称恒等式を別実装で展開し直した。
第九段 entry=N9/(U F9) と `U F10=C10 N9[58,56]` を使うと、
固定第十段 entryは
\[
p=\frac{N9[i,j]N9[58,56]-N9[i,56]N9[58,j]}{U F9 N9[58,56]}
 =\frac{P7\,F11}{U^2F10}
\]
となる。P7のliteralは `5q0²+6q0+1=(5q0+1)(q0+1)=Q R` で既存単元である。
F9は式で約されても、歴史的 chartと補集合には残す。

従って旧 chart内の D(F11) では pが単元で、逆元は
`U²F10/(P7 F11)` である。有限個の旧非零因子とF11の積は同じLaurent domainで非零なので、
この開集合は非空である。旧monomial units、Q/R/E、F9/F10、detA*detB、
一回のP290と全旧・新補集合を保持する。因数分解や互いに素性を仮定しない。

現行 manifestの full167/top101/selected14/columns14、追加RHS87、零制約RHS66を照合した。
元49単元＋制約10単元＋既存typed10単元で69個あり、今回の存在証明は
**typed第11個／元83列の第70個**に対応する。
row6とcol53の単位pivotを使う形式的なSchur後は full166×13、top100×13、selected13×13である。
全87追加RHSと66制約RHSは残る。selected正方行列をfullのcokernelとは同一視しない。

各残行の更新式 `Snew_ij=Sold_ij−Sold_i53*Sold_6j/p` と同じRHS更新、
`x53=(hold_6−ΣSold_6j*xj)/p`、旧69単元への逆代入が正則に定義できる。
FULL行列での単位pivotなのでFULL cokernelが保存される。
これは形式的な存在・復元手順だけの証明であり、実際の更新係数・RHS・
全逆行列・尽くしは計算していない。

独立 SOURCE所有上界49,491、実行0.0364秒、55秒・2 GB内で終了した。
再現コードは `code/s61_col53_nonzero_pivot_chart_source_review_20261105_v64.py`。
