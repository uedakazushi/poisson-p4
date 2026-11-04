# 群座標の有理逆recipe：transpose・scalar・Cayleyの独立レビュー

2026-11-04。対象は実第一jet像B_NEW15。原著producerを再実行せず、元の
24 covariant triple indexから有限pullbackを別に構成した。
原著とは反対のtransvection `g=I+t E01` を用い、全24×16 primitive basisを検査し、
二つの非可換symbolic transvectionで `U(g)U(h)=U(hg)` を確認した。
従ってpullbackのanti-action規約は実装と一致する。

C/Dは対称二次係数行列 `Q(b)=bᵀ S b` として扱い、原著の16個のβ変数の代入展開とは
別に `Uᵀ S U` を用いて全4352原係数を確認した。
元の二つのCartan characterについても全非零二次係数を検査した。
変換則は正確に `C(g*beta)=gᵀ C(beta)g^(-ᵀ)`、Dも同じ。
gそのものによる共役へ置き換えてはいけない。

## 有理intertwinerから平方根なしでgへ

独立受理済みの五quotient座標と実rational sectionからbeta_barを作る。
Rosenlicht generic orbitにより同じ幾何的orbit上にgが存在し、generic stabilizerはidentity。
ここでは最初から有理gの存在を仮定しない。
simultaneous commutantがscalarなので、32×16 intertwiner systemの解は
`T=t gᵀ` の一次元lineである。rank15のCramer chartと非零な一kernel座標の正規化は
有限の有理操作で定まる。全chartの辞書式選択も有限であり、genericな一chartが非空で十分。

標準Jは `J²=−I,J01=1`。symplectic条件から `g J gᵀ=J` でもあるので
`lambda=(TᵀJT)01=t²`。三つのcovariant indexを持つ実有限pullbackはscalar degree3であり、
`(Tᵀ)*beta_bar=t³ beta`。非零beta座標で割った `k=t³` により

```text
g=(lambda/k)Tᵀ.
```

平方根・立方根は不要で、Tのscalar正規化を変えてもこの式は不変。
これらのscalar冪は全16basisのsymbolic pullbackでも確認した。新しいβ点は使わない。

## Cayleyの十変数

`g+I` が可逆なら `A=(g−I)(g+I)^−1`。
普遍的な行列恒等式によりsymplectic条件と `AᵀJ+JA=0` が同値になる。
`A=JS`、S symmetricと書け、その十個の上三角成分が群の有理座標となる。
逆は `g=(I−A)^−1(I+A)`。chartはg=Iを含むのでdenseである。
五つのquotient座標を加えた十五個のfinite rational inverse coordinateを得る。

scopeはB_NEW15のdense chartに限る。全boundary、全25次元Poisson成分、
新しいfull S61 fibreの主張をこのpacketへ混ぜない。
