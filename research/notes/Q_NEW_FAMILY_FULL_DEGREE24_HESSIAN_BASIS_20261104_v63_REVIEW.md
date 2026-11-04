# 全Inv24基底の独立係数検証

2026-11-04。同じ旧 ambient 三ベクトルによる形式平面だけを使う。actual NEW
点や別平面の評価は無い。元の16次元表現、Omega、元の4×4 Lie基底の Gram
逆を独立に読み、C,Dの全二次係数と元のmomentによるthetaを再構成する。

H16は著者の低rank式から作らず、受理済みの全117項 quartic F4 から
`H16=Omega^-1 Hess(F4/2)` を直接微分して作る。各formal quadratic係数行列で
旧式 `theta+Omega^-1 V Gram^-1 V^T` とも一致する。これは4×4 Dへの置換では
なく、元の16×16 Hamiltonian Hessianである。微分式より元の群共変性があり、
全trace wordは普遍不変式となる。

別C++ engineは(t,s)変数順、LEXのnmod多項式で全wordを再構成する。pencilは
Newton trace式を使い、著者のpermutation determinantを再利用しない。
全degree12の34列、degree16の96列、degree24の70列を保存済みの全係数と
比較する。その上で独立Python消去により9、20、67角の係数minorを検査する。

minorは各々20149、2363、26663 mod32003。係数は元のQQ polynomialの
p単元な還元なので、この非零minorはQQ線形独立性を示す。受理済みの
character上界9、20、67により、degree12/16/24の全不変式空間の基底となる。
degree24の選択は70列から62、68、69を除いた67本。

旧C/D-only rank8停止は、その候補・平面でのfrontierとして保持する。
普遍的なrank上界8とは使わない。本稿は全Inv24のbasisを供給し、NEW上の
kernelやdegree24方程式の存在は計算しない。70列には普遍関係3本があるため、
後続source検査には独立な67列を使う必要がある。

```sh
. /workspace/cloud-setup/activate.sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_full24_basis_review_20261104_v63.py --output work/full24-basis-review --work work/full24-basis-review-native
```

60秒・2GB上限。新しいsource wordの展開も高次数probeも行わない。
