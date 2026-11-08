# Projected Veronese 全十五 tangent 源の SL2 exhaustion：独立確認

2026-11-08、PRIMARY E。**独立 PASS：指定全 tangent rank四・被約既約 quintic 枝は空。**
新componentは零。元全15源・whole Poisson cone上界を使い、選んだ Poisson point から全支持を推論しない。

## 全 source、実群作用、全 constant 算術

既全280×160 tangent条件と canonical Euler商により、全 source V は十五次元。
全Poisson cone X⊂V の affine dimension≤2 は、旧別作者が全54元Jacobiと
全10079 polynomial lifts、全242 leading monomialsから確認済みの入力である。
今回そのGB、Jacobi又はPfaffianを新しい点へ評価していない。
同じ surface 上の全rank零五源のPoisson支持は原点だけ、という別の旧独立定理も利用する。

非退化 projection center a²+b²+c² を保つ実 SO3≅PGL2 が、全 projected Veronese と V と X を保つ。
独立 checker は QQ(i) の Fraction 演算で、元 frozen15 basisから全 push action
`AP+PA^t-(Ax)·∇P` を直接作り直した。
三つの full150×15 action、source15×15 matrix、15×150 left inverse の全積が一致する。
実 quadratic embedding と center の保存も別に確認した。

保存された全 chain と両側 inverse により、全 source が正確に

```
V=Sym² ⊕ Sym⁴ ⊕ Sym⁶
```

へ分かれる。各 n=2,4,6、j=0,…,n で
`Hvj=(-n+2j)vj,Fvj=v(j+1),Evj=j(n-j+1)v(j-1)`、両端零を全て確認した。
Hの固有値だけから分解を推測せず、十五列全体とその inverse を保持する。

全 rank零五basis を全十五 source へ戻し、Sym4 block以外は全零、5×5 block両inverseが一致。
三つの実 fundamentalfield wedge を divergence-zero projector で canonical化すると、
全て Sym2 blockへ入り、その3×3 block両inverseも一致する。
従って両者は包含だけでなく全spaceとしての等号である。

physical coordinate変換と逆も確認し、H weights は `(4,2,0,-2,-4)`。
embedding上の座標は `(u²,uc,uv-2c²,vc,v²)`、u=a+ib,v=a-ibである。

## 全支持から Borel eigenline への矢印

Xはhomogeneous coneなので、非零点の射影支持P(X)はdimension≤1。
Gは連結で有限個の既約componentを各々保つため、任意のG-orbitもdimension≤1。
従って全点の stabilizer のidentity成分はdimension≥2。
PGL2の二次元連結真部分群はBorel共役、三次元なら全群なので、いずれも共役Borelを含む。

Borel不変line上の unipotent radicalのcharacterは自明だから EΠ=0、torusから HΠはscalar倍。
完全分解のker Eは最低weight−2,−4,−6の三lineで、それらのweightは相異なる。
したがって全Borel eigenlineは厳密に三lineのいずれか。
同weight混合、未分類のgeneric orbit又はsource parameterは残らない。
原点はrank四でなく、scheme厚みから別の幾何閉点は生じない。

## 三branchの全排除

Sym2は三次元Lie algebraのΛ²の像。
三次元spaceの任意の2vectorはdecomposableなので、実二fundamentalfieldのwedge、射影rank≤2。
canonical Euler補正はこの射影rankを変えない。

Sym4は全rank零五源そのもので、既全rank零Poisson定理により非零Poisson点を持たない。

Sym6の最低lineはH-character−6、Pfaffian mapはSL2にequivariantな二次map。
ambient determinantは一なので、actual quintic F の各physical degree五monomialはweight十二。
y0無しのweightは高々2×5=10であるから、必ず F=y0*F4。
F零ならrank四でなく、非零ならproper一次因子を持つ。
反対Borelではy4因子となるだけであり、全actualG共役とscaleでも因子性を保つ。
この議論はactual Pfaffianの新しい評価を必要としない。

以上で全三branchを尽くす。全smooth projected Veroneseの単一PGL5軌道と全逆を戻して、
全埋込みの指定 rank四・被約既約 quintic tangent枝を排除する。
可約rank四sourceの全既知成分帰属、他の曲面又は全分類完了は主張しない。

## 保存と有限確認

実行前gate PASS、ONE Gaussian Fraction batch、0.777秒／RSS17696KiB。
全91項目の行列・多項式恒等式を確認し、producer/CAS import、Jacobi/Pf point、新kernel/rank/GB無し。
実 inline source は jobとreportの両方へ固定した。

- job：`jobs/projected_veronese_tangent_sl2_separate_review.json`
- receipt：`reports/CLASSIFICATION_PROJECTED_VERONESE_TANGENT_SL2_SEPARATE_REVIEW_20261108.json`

旧dimension/zero-source/producer bytesを変更していない。
