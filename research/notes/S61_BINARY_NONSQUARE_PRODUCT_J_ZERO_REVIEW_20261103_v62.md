# 独立review：J=0上のproduct欠損の有限性

2026-11-03。`S61_BINARY_NONSQUARE_PRODUCT_J_ZERO_20261103_v62.md`、
同名producer、全input/generated hashとQQ minor sourceを確認した。
48.82秒の四minor計算自体は再実行していない。

**受理する結論は、aFBCd≠0、J=5ad−6e−4=0上でproduct-rank41欠損はQQ上次元0以下。**
全三coreの全点階数や、四次・ambient冪零性は含まない。

J=0上のd=2(3e+2)/(5a)とB=2Bs/(5a)は正しい。
ℓの指定係数−32(3e+2)/27はこの開で単元なので、modℓの六変数表示は正確である。
分母とrow contentの除去、およびnativeの12unit pivotは、列挙された既知単元だけを反転している。
未知のcore因子を反転していない。

残る12×44行列の四つの実際の12minorはQQ多項式である。
Bareissの全exact divisionを確認するsourceと、報告されたgcd
a^12 C^3(3e+2)^31（非零有理定数倍）を照合した。
既知単元を反転した二変数UFDでは、全四minorを含む高さ1素idealは存在しない。
従って全rank欠損集合は次元0以下となり、体を代数閉包へ拡大してもこの上界は保たれる。
独立QQ数値determinantの二点×四minor照合は追加検査であり、gcd命題の代わりにはしていない。

modℓの24行が独立な点では、42行の線形syzygyをmodℓへ落とすと
四つの係数はciℓとなる。ℓで割ってZ=W=0へ制限し、既存pair-coverの四二次式独立性からci=0。
最後にuW+vZ=0のsyzygyはWZ−ZWの一つだけとなるので、42−1=41である。
この段階でℓ=0上の商だけの結論をambientへ持ち上げてはいない。
