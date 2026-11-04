# u4,target15：独立な記号 mod32003 多項式監査

2026-11-04。先の QQ `P/a^k` 検算は48秒ガードで停止し、その失敗資料は変更しない。
新たに認可された別の方法として、同じ一成分だけを記号環
`F32003[a,e,d,a^-1]` で検算する。a,e,dには値を代入しない。
新しい点・素数・成分・小行列・因数分解・rank探索は使わない。

独立に源対応を検査した768演算入力の同一SHAを固定する。実行前に全QQ分数の
分母が32003で非零であることを確認する。各実source原子の分母はQQ単位とa冪であり、
今回の closure にB分母はない。新しい nmod_mpoly 演算器は源768演算を係数ごと
記号的に再計算し、公開された全7022 QQ係数の還元と照合する。

actual pair2から3 source cofactors、全5 Bhat 成分を再構築する。
さらに producer の Adj4 contraction と異なり、同じ actual bordered5 行列の
120 permutation determinant を独立に計算して照合する。
P4の固定24 permutation determinant、`D4=a^4 C^3 R434`、
`delta2=-(6/5)B p²`、`Delta5=delta2 delta3` もこの記号環で確認し、
a,p,D4,R434,Delta5 と全named分母多項式が非零であることを必要条件とする。

結果の記号多項式が非零で、すべてのQQ分母がp-unitなので、標数0の actual Saux も非零。
受理済み actual hybrid scaling により generic residual位置(0,9)の一成分の非零性が従う。
これは mod p における実source恒等式と標数0非零性の独立証明である。
**QQの全7022厳密係数恒等式を独立に検証したという主張ではない。**
その exact QQ 値は引き続き primary の計算として区別する。

型付き chart の領域は `FULLU∩D(R434 Delta5)` のまま。残6行の階数6、
global cubic57、quartic239、他成分、chart補集合の除去は主張しない。
実行は60秒・2GB・32MiB・150000 live係数を超えない固定検査に限定する。

再現：`.venv/bin/python code/s61_nonsquare_prePOWER_u4_t15_modp_polynomial_review_20261104_v63.py`
