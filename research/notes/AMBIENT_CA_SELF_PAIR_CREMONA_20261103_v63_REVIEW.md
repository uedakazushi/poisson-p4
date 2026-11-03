# CA self-pair Cremona・有理葉・旧Log/LPB sourceの独立検証

2026-11-03。対象は `AMBIENT_CA_SELF_PAIR_CREMONA_20261103_v63`。
producerを再実行せず、生のCA、元座標上の有限線形作用素、逆方向tensor変換、
直接dlog二形式によって検証する。後続のN/CB所属packetは対象外とする。

## 実CAと一様な座標変換

生CAへ `f=−x²−3ahx`、`g=0`、Hamiltonian
`−xhw+h²z+z³/3−axzw` を代入し、divergence canonical projectionを計算する。
CAの全20基底も独立に作り、元paired-sourceのfixed20とcanonical150を照合する。
実際のテンソルは `∂y wedge V`。商座標 `(h/x,y/x,z/x,w/x)` のJacobianから
射影chart bracketを直接求め、x=1で

    ∂y wedge D,
    D=(1+3ah)∂h+(h+az)∂z+(h²+z²−aw)∂w

を得る。

Qは作者の固有単項式公式から再利用しない。元のh,zの次数2以下六次元空間において
`(D+a)Q=h²+z²` の6×6線形系を解く。その行列式は `840a^6` であり、
全a≠0で一意に解け、追加の例外parameterを除く必要がない。
`H=h+1/(3a), Z=z−h/(2a)−1/(2a²)` に変換したQが作者の保存式に一致する。
`W=w−Q` は `D(H)=3aH,D(Z)=aZ,D(W)=−aW` を満たし、
この三角多項式自己同型と逆のJacobianは1である。

## actual Poisson Cremona (2,2)

独立に求めた自己同型を斉次化し、両合成がx³/X³倍の恒等写像であることを確認する。
forwardの最後の座標のz²係数は `−1/(3a)`、inverseの最後の座標のZ²係数は
`1/(3a)`。残りの座標にはx²/X²があるので、全非零aで共通因子はなく、
bidegreeは正確に(2,2)である。affine・homogeneous多項式mapの全分母はaの冪だけである。

作者がforward pushforwardを確認したのに対し、独立コードはhomogeneous-degree-one
inverse mapの全5×5 Jacobianでdiagonal canonical tensorを引き戻し、元Pとの全係数一致を
確認する。Jacobian1とEuler保持も確認する。函数倍で置き換えたfoliationだけの同値ではない。

## 完全なCasimir体と葉の次数

`f=H/Z³,g=ZW,t=Z,Y=y` は函数体のbirational座標であり、保存された逆式を直接検証する。
この座標でPoisson tensorは `at ∂Y wedge ∂t`。
char0の有理函数の微分定数から、そのPoisson中心は正確に `C(f,g)` である。
c1c2≠0の葉のchartは `C times Gm`。

y方向を除いたP3曲線のparameter sectionの次数は `(4,1,4,7)`。
最高係数は `c1²(4a²+1)/(28a³)`、t6係数は零、t5係数は `c1/(5a²)`。
特殊化は数値点でなく `QQ[a]/(4a²+1)` で正確に行い、次数が `(4,1,4,5)` と確認する。
有限t≠0ではx=t、t=0では最後の座標c2が非零。無限遠ではそれぞれ非零t7/t5項がある。
従って各場合の同次sectionにはbasepointがない。`t=B/x` が有理逆なのでmap degreeは1。
pullback O(1)から像の次数は7、または4a²+1=0で5となる。
yを自由に加えたP4の葉閉包は同じ次数の有理coneである。
特殊Casimir値c1c2=0についてこの一般次数を主張しない。

## 実Log1112/LPB sourceとの一致

divisors `(A,B,x,C)` に対する `dlog f,dlog g` のresidueは
`(1,−3,2,0)` と `(0,1,−3,1)`。両方ともdegree vector `(1,1,1,2)` に直交する。
旧 `log_bivector` sourceの係数辞書をASTから抽出し、parameters
`(1/4,−3/4,7/4)` でこの二つのresidueの外積と全六係数が一致することを確認する。

直接 `ABxC*dlog(f) wedge dlog(g)` を計算して全10係数を照合する。
`ω=−ΩP/a`、外微分とvolume dualを取った全150係数で `Plog=−5P/a` が成立する。
符号はaffine `(H,Z,W)` のvolumeに対して
`i_D(dH wedge dZ wedge dW)=−aω` で確かめられ、三角Jacobian1で元座標にも戻る。
LPBは一定因子eyを旧center e4へ送る置換によって、yに独立な二次場との外積という
旧族の定義そのものに入ることを全150係数で確認する。

旧Log1112_28とLPB39の成分次元を証明する大計算は、元の受理済みcertificateを依存証拠として
再利用する。source/tangent/family記録と全hashを照合する。全a≠0でのactual所属と
多項式族のa=0での非零性、成分の閉性から、a=0を含む全直線も両閉成分に入る。
CB/N所属、全self-pair locus、conductor、新成分数については結論を加えない。

独立コード：`code/ambient_CA_self_pair_cremona_review_20261103_v63.py`。
証拠と全入力hash：`reports/AMBIENT_CA_SELF_PAIR_CREMONA_20261103_v63_REVIEW/`。
