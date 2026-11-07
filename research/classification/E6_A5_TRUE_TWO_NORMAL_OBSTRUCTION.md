# E6/A5の真二法線：全二次余核による最弱支持判定

2026-11-07、PRIMARY E。既canonical点を変更せず、未計算だった全二法線の二次obstruction mapを一回だけ扱う。
高次jet、別anchor又は部分方向の列は用いない。

## 全実群と真の線形切片

150係数のraw空間にはEuler垂直24方向の加法群、GL5及びscalarが作用する。
全350projective Jacobi idealはこの作用で保存される。
特にEuler∧linear fieldの加算は射影bivectorを変えず、全350式自体を変えない。
既受理の同じ点Πでこの実作用の微分像Fはdimension48、全Jacobi微分Dはrank100である。
sourceの23軌道とその全六磁気閉包は既に独立確認されている。

Fの48独立列と48pivot行を取り、その行の変動を零にする102次元の線形切片をCとする。
FとCはraw150の直和である。
Dの全kernelからFを引いてC内の二法線N1,N2を作る。
Dの100独立列に対応する標準方向もFへ射影してC内のU1,...,U100とする。
保存するunitは

```
F48 + N1,N2 + U100 = 全150の直和
```

の非零determinantである。D*F=D*Ni=0、DUは全100像を張ることも全行で保持する。

48個の実群parameterと切片Cのaction mapは、点Πでambient微分が可逆である。
従ってformal inverse function theoremにより、全raw Poisson germは形式的に
48次元smooth factorと真cutの積である。
実群作用が全idealを保存するため、この積分解は単にtangent spaceの分解ではない。

## 全350からの二次map

二法線係数をu,vとする。全projective JacobiをΠ+uN1+vN2へ代入すると、
constantとlinearは零で、二次係数は350×3列

\[
Q=(Q_{11},Q_{12},Q_{22})
\]

となる。uv列には正確に交差項を用い、二倍を落とさない。
全D像の100列・100行のunit blockをBとすると、全余核射影を

\[
\bar Q=Q-D_{\mathrm{cols}}B^{-1}Q_{\mathrm{rows}}
\]

で計算する。全350行を残し、選んだ三行だけで全mapを代替しない。
100implicit変数を解いたformal残式の二次leading termsは、まさにこの全Qbarの行である。
その全行が生成するbinary quadratic idealをI2とする。

rank三ならI2は(u²,uv,v²)。rank二なら二独立binary quadraticsの全homogeneous degree三Macaulay行列を使い、
そのdetが非零ならu³,v³の明示ideal-membershipを保存する。
これによりP1上の共通零が無いことを、無限遠を捨てずに確認できる。
rank一以下又はdet零では、v=1の共通因子とv=0の可能性を全て保存して、この二次routeを止める。

## m-primaryの場合だけの全成分性

I2の根基が(u,v)なら、全formal残式のinitial idealがm-primary idealを含み、真cutの被約dimensionは零となる。
全raw local dimensionは48以下で、実orbit下界48と一致する。
実orbit閉包は既約なので、そのprojective23軌道閉包が全ambientの既約componentである。
全六磁気source閉包もそれに等しい。

各全global componentはconnected実群によって保存されるため、点Πを通るcomponentは実orbit閉包を含む。
従ってこの支持判定はformal sliceだけの架空のcomponentを数えていない。
一方tangent dimensionは50なので、scheme smoothnessや精密nil厚みを主張しない。

共通P1根が残れば、この二次必要条件からcomponenthoodもactual変形も結論しない。
その根は真Poisson curveの存在証明ではない。高次への延長は行わない。
