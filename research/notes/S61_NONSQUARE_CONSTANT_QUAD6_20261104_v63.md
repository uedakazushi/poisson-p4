# 全aC開での六quadの定数階数6

2026-11-04。これは22 quadratic syzygyの独立性をFULLUの外へ延ばす主張ではない。
元の六quadそのものが定数上独立であることを、旧証明の残るcaseも含めて証明する。

既受理のZW商の四quadについて、五変数linear cofactor mapは全aC開で
20×35、rank19である。もし四quadに非零な定数関係があれば、それに五つの
独立な変数を掛けた五方向がkernelに入り、kernel次元1に反する。従って
ZW商の四quadは全aC開で定数上独立である。

元の六quadを `(q0,q1,q3,q5,ell*Z,ell*W)` とする変更のdeterminantは−80、
分母はaの冪だけである。ZW制限で上の四式が独立なので、残る問題はell≠0だけである。

- J≠0なら、既受理の `ell_x1=-16*C*J/243` が非零。
- J=0かつd≠0なら、`Delta=ell|ZW-2J*L/1215`に対する既受理の
  `3*Delta_x2-(7/2)*Delta_x4=-80*a*d²/81` が非零。
- J=d=0ならe=−2/3。既存の実ell全係数を別Fraction演算で特殊化すると
  `ell_x0=3200*a/243`となり、これも非零。

従って全aC開でell≠0であり、`ell*(bZ+cW)=0`はb=c=0を強制する。
これで六quadの定数階数は全域で6となる。

実weighted chart aFBC開でも単位frame変更を戻せるから、degree2の29方向
（七lightのSym²とheavy一方向）の商は正確に23次元である。旧Z2欠損は
この開で空となる。full cubic57、四次239、22方向のd0/J0での独立性は
結論しない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python \
  code/s61_nonsquare_constant_quad6_20261104_v63.py --output work/quad6-new
```
