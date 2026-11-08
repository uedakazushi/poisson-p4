# rank一long新三源：generic軌道23,23,21と商2,2,4

2026-11-08。全sourceのpure stabilizer補足。既一点のrankを全域へ外挿せず、共有source及び既証明書を変更しない。新point／minor／profileは実行しない。

## 1. 定数gの二共鳴源

全μ=1及び1/3のsourceは同じmodular場X=∂x+x∂yを持ち、[X,P]=0。Xは非零affine場であり、P4上の非零projective linear Poisson場へ延長する。従って全sourceでPGL5 stabilizerは少なくとも一次元、全軌道次元は23以下。

保存した各一固定点はEuler／GL／scalar行列のrankが48、Euler24なので実軌道は48−24−1=23次元。非零48-minorは非空開条件であり、各既約25次元source像のgeneric軌道は正確に23。generic Rosenlicht quotientは二次元となる。

## 2. g零scalar-cubic全source

q0=(z0,w0)、U=−Euler(q−q0)/3、Hは中心q0の任意binary cubic、Fは全二次式とする。d=−2/3、g=0で、元六括弧は

```
Pxy=−F−(y−x²/2)/3,
Pxz=−(z−z0)/3, Pxw=−(w−w0)/3,
Pyz=−x(z−z0)/3+H_w,
Pyw=−x(w−w0)/3−H_z, Pzw=0.
```

Xに加え、全sourceで

```
Yz=(z−z0)∂y, Yw=(w−w0)∂y
```

もPoisson場である。一般のℓ∈span(z−z0,w−w0)について、L_{ℓ∂y}Pのxy成分は

```
ℓ*∂y(Pxy)−U(ℓ)=−ℓ/3−(−ℓ/3)=0.
```

yq成分の残る項はPzw=0により零、他の全成分も零である。FとHはこの恒等式のどこにも追加条件を課さない。

X,Yz,Ywは独立なaffine場で、h-component零の斉次行列としてEuler方向も含まない。従って全sourceでprojective stabilizerは少なくとも三次元、全軌道次元は21以下。

保存した唯一のSC点のM rankは46なので実軌道は46−24−1=21。これと全source上界によりgeneric軌道は正確に21。25次元成分のgeneric Rosenlicht quotientは四次元となる。

H squarefree条件は全25次元像／唯一double-Hを確認する非空開であり、上のstabilizer証明には不要。全binary重根境界及び全Fも保持する。全点のorbitが同じdimensionだとする主張、全orbit分類又はquotientのcompactificationは含まない。
