# rank零plane・traceもtracelessも零・非零binary symbolは全rank二以下

2026-11-07、主義務G。producer手証明。全sourceの逆変換と残余を保持する。

**命題。** rank零plane P=(z=w=0)を持つquadratic projective Poisson tensorを考え、
conormal Lie algebraは可換とする。
weight零のtrace及びtraceless co-Higgs symbolがともに零で、vertical binary quadratic Kが非零なら、
全higher tensorはrank二以下である。
従って実被約既約quinticを持つGの枝には現れない。
この適用では全重複度のconormal可換性を新modular補題から得るが、命題自体はその可換性を仮定するだけでよい。

## 1. 全sourceと全weight一forcing

normal GL2とtensor scalarにより、全非零binary quadratic Kはzw又はz²へ直る。
元の全nonzero scalarsとnormal matrixの逆で元のsourceへ戻せる。
weight零は`Π0=K∂z∧∂w`。
full higher sourceは既存のEuler quotient三十/九、すなわち

`Π1=zL+wM+z²D∧∂w+zw(B∧∂z+E∧∂w)+w²C∧∂z`

と、全base-baseのnormal-quadratic Π2である。
L,Mは三つのbase-pairsそれぞれ全linear-base係数（各九）、B,C,D,Eは全constant-base vectors。
Euler_baseをEbとする。

全weight一Jacobiは次のforcing coefficientを持つ。
K=zwならbase²normal²のquad成分は

`z³ Eb∧D+w³ Eb∧C−z²wL−zw²M`。

異なるnormal monomialsを比較するとD=C=0、L=M=0。
Eb∧constant=0はconstant零を意味する。
残るB,E全六係数は全Jacobiを満たす。

K=z²なら同じ全quad成分は

`z³(−Eb∧E−L)+z²w(Eb∧B−M)+2zw² Eb∧C`。

従って全解は

`C=0`, `L=−Eb∧E`, `M=Eb∧B`、B,D,E自由。

base³normalの残るquad成分はEb∧L、Eb∧Mであり、この解では自動で零。
よってこれは全kernel九係数であって、smooth cubicのinjectivityを誤って適用して六へ減らしていない。
これらのforcingは`classification_plane_smooth_cubic_weights.py`の全symbolic式で
cubic fを零としたものからも直接得られるが、本noteはそのzero-cubicの場合を全sourceで比較した。

## 2. actual有限shearとweight二injectivity

actual可逆base change `x→x+zY+wZ` はΠ0を保ち、first weightへ

`K*(Y∧∂w−Z∧∂z)`

を加える。Π2への追加も全九係数のsource内に入り、逆changeはminus。
K=zwではこれはB=−Z,E=Yという全六kernelへ同型なので、全Π1を消せる。

K=z²ではcanonical Euler representativeへ直すとD=Y、E=Z、L=−Eb∧Z、B=M=0となる。
従って任意のD,Eを消し、全Π1をBだけの

`Π1=w Eb∧B+zw B∧∂z`

へ直せる。first-weightへのshear作用はlinearであり、元Π1への追加作用はweight二に入るだけなので、
このnormalizationは元の全higher sourceを保つ。

全weight二D2について、base²normal²成分には厳密に`−2KΠ2`が現れる。
よってK非零ではfull D2はinjectiveである。

## 3. 残った全sourceのrank

K=zwではΠ1=0なのでweight二がD2Π2=0となりΠ2=0。
Π0だけが残るからrank二以下。

K=z²では上のcanonical Π1はprojectively

`R=B∧W`, `W=2zw∂z+w²∂w`

に等しい。実際、その差はEuler_total∧wBの符号倍。
Bはconstant-base fieldなので[B,W]=0、従ってRはstrict Poisson。
さらにW(z²)=4z²w、div_normal W=4wなので
`L_W(z²∂z∧∂w)=0`。
従ってΠ0とRはstrictにcompatibleであり、Π1自身のprojective Jacobiも零。
weight二は再びD2Π2=0、injectivityからΠ2=0。

全tensorはprojectivelyΠ0+Rで表され、そのimageは
`span(B,∂z,∂w)`というdimension高々三のspaceに入る。
skew tensorのrankは偶数なのでrank二以下。
actual finite shearとnormal GL2及び非零scalarで戻してもrankは同じである。

## 4. 正確な補集合

T非零のscalar co-Higgsはこの命題へ含めない。
traceless symbol非零のglobalV源も別の義務である。
T=V=K=0ではΠ0=0でm四sourceが残り、このforcing/injectivityを適用しない。
全rank零plane又は全Gの完了とは主張しない。
