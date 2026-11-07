# 被約退化因子・rank零曲面の可換conormal補題：別作者独立確認

2026-11-07。reviewer `generic_exhaustion`、producer `root`。主義務G。
入力は`REDUCED_DEGENERACY_ZERO_SURFACE_ABELIAN_CONORMAL.md`。
本記録は手による局所幾何証明の独立確認であり、CAS又はplan checkerによる定理検証とは呼ばない。

**受理する射程。** smooth complex fourfold上のgenerically symplectic Poisson tensor πの
実退化因子Dが被約であるとする。rank零locusに含まれる任意のintegral surface Sでは、
Sの全smooth locus上のconormal Lie bracketは零である。
Sに沿うDの重複度二、Frobenius性、Sの大域smooth性、Dの既約性を仮定しない。

## 1. 独立な局所証明

rank零locusにcodimension一の既約成分Yが存在すると仮定する。
Yのsmooth一般点においてambientは正則、Yは一元gの零集合なので、
全六局所係数はgで割れπ=gρとなる。
したがってπ²=g²ρ²であり、nonzeroな実退化因子にYが二回以上現れる。
Dの被約性に反する。この議論はYが別の曲面成分と交わる点でなく、その一般点で行うため、
rank零schemeのembedded component又は非被約性を仮定せずに済む。

従ってrank零被約台のdimensionは高々二であり、Sはその既約成分である。
Sのsmooth一般点近傍を縮め、rank零被約台がSだけとなるopenを選べる。

任意の局所nonvanishing holomorphic volumeについて、Hamiltonian fieldのdivergenceは
局所modular vector field Zを定める。二階微分項はπの反対称性で消える。
Hamiltonian commutatorのJacobi恒等式とdivergenceのcommutator公式を用いると、
ZはPoisson bracketのderivationである。Hamiltonian fieldをどちらの符号で定義しても、
Zの符号が変わるだけで`L_Z π=0`は同じである。
この式はπの零idealを保存する。あるいはZの十分短い局所holomorphic flowがπを保存する。
上のopenではrank零被約台はSだけなのでflowもSを保ち、ZはSへtangentである。

Sのsmooth点でcoordinatesを(a,b,z,w)、I=(z,w)とする。
全π係数がIに入るので、base-normal係数のa,b方向微分はS上で零。
`{z,w}=l1(a,b)z+l2(a,b)w mod I²`と書くと、
`Z^j=Σ_i ∂_i π^{ij}`という符号を選んだ全normal成分は

`Z^z|S=−l2`, `Z^w|S=l1`。

したがってtangencyからl1=l2=0である。
任意のI内の二関数をz,wで展開しLeibniz則を使うと、残りのmixed項もI因子を二つ持つため、
`{I,I}⊂I²`を得る。これはSのdense smooth openでの証明であるが、
conormal bracketはSのsmooth locus上のvector-bundle morphismなので、Sのintegralityから全smooth locusで零。
Sのsingular点でのconormal moduleのtorsionについて、より強い主張は加えていない。

volumeをunit倍してもZの差はHamiltonian fieldであり、π|S=0からS上で零である。
normal成分の消滅は選んだvolumeに依存しない。

## 2. 使用出来る帰結と使用出来ない帰結

Sがplaneの場合は全点smoothなので、全normal-normal linear coefficientが消える。
従って以前m二に限定したweight零35・weight一30・weight二9のprojective sourceと、
そのJacobiだけから導くco-Higgs O(0)/O(1)分岐及び全Lax式を、m三以上にも適用出来る。
これは古い証明書のbytesを変更することでなく、新しい局所補題を新たな適用入口として追加するものである。

ただし旧m二normal formで除外したleading-Pf零、τ零、resonant row零、rank低下を、
m三以上で除外することは出来ない。
それらを復活させた全parameter source、全higher Jacobi、dimensionと補集合は別の義務である。
本補題単独でm三以上が空、G完了又は既知成分帰属とは結論しない。

また本証明は**局所modular fieldのnormal成分**だけを使う。
projective quadratic tensorの任意のEuler liftのdivergenceを大域stabilizerとして固定したり、
first symbolのmodular fieldを全higher deformationへ延長したりしていない。
以前に棄却した大域divergence shortcutを必要としない。

以上、producerの全仮定、proof、m三以上への入口、未解決補集合を独立に読解してPASSとする。
