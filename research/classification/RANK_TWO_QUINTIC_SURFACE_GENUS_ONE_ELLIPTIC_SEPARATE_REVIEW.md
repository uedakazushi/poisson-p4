# 次数五・rank二・種数一 elliptic 正規化：全排除の独立確認

2026-11-07。producer `plane_collision_review`、独立確認 `root`。
対象は `RANK_TWO_QUINTIC_SURFACE_GENUS_ONE_ELLIPTIC_EXCLUSION.md`
（SHA256 `1ea3d95ad56fc5f8b1792232f49208e53d68ffc197d09e3b5b7abbd9daea45b9`）。
**記載された全 elliptic 枝の不存在を受理する。** 新しいCAS計算は無い。
次数五・rank二の canonical index二正規化は別の未解決枝である。

## 正規化源と非ample境界

Fukumaの既確認 Proposition 3(iii) は smooth quasi-polarized resolution に適用する。
common resolution上のpolarizationの一致から全section ringが一致する部分は
rank零を使わず、今回もそのまま成立する。従って elliptic 型では
`Y=P_B(E) -> X`、`H=f*L`、rank E三、degree五、E globally generatedを得る。
X自身を最初からprojective bundle又はsmoothとは置かない。

Eのdegree零quotientについて、global generationは非零`E -> O_B`の全射性と
1をliftするsectionを与えるので、`E=O_B^r + F`、`H0(F*)=0`までsplit出来る。
全有限cover上のH-degree零曲線も取り落としていない。
そのnormalization Tのtautological quotientはdegree零・globally generatedでtrivial。
evaluationの係数はcompact connected T上の定数であり、Fからの非零部分があれば
surjective `T -> B`を通して全evaluation kernelを消す定数functionalがBへ降り、
禁止された`F -> O_B`を与える。このため零次数曲線は全てzero-block sectionだけにある。
同じ証明でP(F)のtautological morphismはfinite、Fはampleとなる。

r零はfinite birationalとXのnormalityから同型。
r一は一本のelliptic sectionのみを一点へ縮める。
proper birational morphismのconnected fiberと、正次元projective fiberは曲線を含むことから、
このsection外の余計なfiberを許さず、exceptional locusはcodimension二である。
degree一生成又はFの非常にampleなembeddingを仮定していない。

r二では実五sectionをzero-blockへrestrictするとrank二である必要がある。
核となる三sectionはF上でbasepoint-freeでなければならない。
従って実像Dは共通line vertexとplane curveを結ぶconeとなり、
その二次元非正規支持はbase curveの有限個の非正規点上のplanesだけ。
次数五integral Sとは両立しない。rank零用cone命題の流用は無い。

## 全導手と留数の矛盾

r零・一ではYとXは全codimension一で同型なので、dualizing generatorのdivisorは
Y上で正確に`-C`、`C ~ 3H-p*det E`となる。
small exceptional curveを追加のWeil divisorとして数えていない。
Pic Yは`ZH + p*Pic B`であり、degree零twistを捨てることも無い。

length公式はdelta S一・他の非正規曲面無しを与える。
cuspならC二倍primeだが、そのH係数三は奇数で矛盾。
nodeではnormalization branchは元の環の商なのでPoisson bracketがregularに降りる。
Hamiltonian derivationがminimal primeを保つことから導手にtangentでもある。
この主張は係数体をPoisson中心にせず、全tangential/mixed bracketを保って得られる。
inverse-residue一形式alphaはclosedで、全prime導手上に高々simple poleを持つ。
そのresidueは定数であり、branchからSへのgeneric separabilityとrank二restrictionから非零。

数値fiberをFとし、`A_n=H+nF`の十分大きい連続する二つの整数nを使う。
基底curve上のdegree一line bundleの選択は全Picard twistを保つ。
相対O(1)に十分大きいbase twistを加えたA_nはvery ampleとなる。
一般の二切断曲線は、Yの三次元性により全codimension二bad locusを避けられる。
従ってcompact curveのresidue theoremをそのまま使える。

導手がone primeなら、そのA_nとの交点は`10+6n`で正。
従って非零であるべき定数residueが零となる。connected nonsplit nodeも含む。
split二primeの場合は各primeのH二乗degreeが五。
prime vertical divisorのdegreeは一なので、H係数は共に正整数である。
全classは交換を除いて
`Gamma1=H+p*M0`、`Gamma2=2H+p*Mminus5`、
`deg M0=0`、`deg Mminus5=-5`となり、全Pic零twistを保持する。
交点数はそれぞれ`5+2n`、`5+4n`。
連続する二つのnの留数係数行列のdeterminantは

`(5+2n)(9+4n)-(7+2n)(5+4n)=10`。

従って二residueとも零となり、rank二からの非零性に反する。
elliptic Yの非零holomorphic一形式を消す必要は無く、GKKPも使っていない。

## 判定

全nef/GG bundle、r零・一・二、全Pic零twist、cusp、single/non-split/split nodeを
上記の相互に完備な場合分けで処理した。
producerの全elliptic枝不存在命題を独立に受理し、Gの一つの全枝を閉じる。
canonical index二、他次数、全G/E及び成分総数の上界は本命題から主張しない。
