# 第33段階後の scalar unit の独立検証

## 対象と結論

`S61_BINARY_NONSQUARE_AFTER33_SCALAR_UNIT_20261103_v62` の証明書を、探索を再実行せずに独立に検証した。第33段階の実際の `24×58` 行列の第7行（0始まり、残存する8本の積の行の一つ）を `r`、保存された58成分のベクトルを `c` とすると、元の多項式環で厳密に

\[
r c=aP^2(5P-6e-4)=aP^2J
\]

が成り立つ。ここで native 形式の第3変数は `P=ad` であり、置換前の `d` ではない。したがって `aPF B_0CJ` を反転した既存の開集合で、次の単位の行が得られる。

## 独立な計算

検証コードは `code/s61_nonsquare_after33_scalar_unit_review_20261103_v63.py` と `code/polynomial_scalar_integer_review_20261103_v63.cpp`。前者は410MBの元行列をストリームで hash し、直前の `SECOND_SCALAR_PIVOT` の出力 hash と一致させ、保存行の全58成分を実際の元行と比較する。cofactor は33成分が非零、合計975項、次数7、係数の分子・分母の最大 bit 数3299である。

後者では行全体の分母の最小公倍数 `D_r` と cofactor 全体の分母の最小公倍数 `D_c` を取り、各成分を整数多項式に変換する。**整数多項式演算**により

\[
\sum_j(D_r r_j)(D_c c_j)=D_rD_c(5aP^3-6aeP^2-4aP^2)
\]

を検証する。producer の有理多項式による内積・因数除算とは異なる受理計算である。右辺は明記した3項と直接比較し、未知の因子を反転しない。小さな別の fixture `(e,ae-1)·(a,-1)=1` は受理し、右辺を2に変えた fixture は拒否する。

再現例（既存 Cloud runtime を activate する）:

```
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/s61_nonsquare_after33_scalar_unit_review_20261103_v63.py \
  --output work/after33-review-report --work work/after33-review-native
```

保存結果は `reports/S61_NONSQUARE_AFTER33_SCALAR_UNIT_20261103_v63_REVIEW/report.json`。全 input hash、元行との対応、整数演算の結果と fixture の結果を含む。

## 線形代数と範囲

`A c` を新しい列として追加しても cokernel は変わらない。この列の第7成分は既知の単位 `g=aP²J` なので、行消去とこの単位による scaling で行を一つ消せる。単に列を追加して pivot する場合の次の形は `23×58` である。`23×57` にするには別の可逆列変換の証明が必要で、この証明書からは主張しない。

この検証時点で巨大な展開済み native cache に対して実行済みの pivot 数は**33**のままである。第34段階の巨大な展開は再実行していない。後続の厳密な回路形式による Schur 証明書が作られる場合は、その別証明書で数える。本証明書は全点での積41本・全 cubic57本・quartic の被覆を主張しない。次数3の modular 探索失敗も、QQ 上の非所属や任意次数の障害には読み替えない。
