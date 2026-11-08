# Twisted-cubic cone rank二の全T38/T39被覆：別作者によるFisher実族と全源の監査

2026-11-08、PRIMARY E、独立review **PASS**。
producer `CUBIC_CONE_RANK_TWO_POINTWISE_T38_T39_COVERAGE.md` と、その全source還元・nodalモデルを別に読解した。
一つの有限Fraction実行で全parameterのliteral構成を再計算した。
新しいpoint、Jacobi operator、rank、minor、kernel、GB又はlocal jetは計算していない。
新component数の増分は零である。

## 1. 再構成の独立性と結果

契約は `jobs/cubic_cone_rank_two_nodal_fisher_separate_review.json`。
active gate PASSの後、一回の実行は0.17553秒、最大RSS15,228 KiB、111個のmatrix/polynomial恒等式が全PASS。
実行sourceは契約と結果JSONにそのまま埋め込んである。
標準ライブラリFractionと、以前の独立checkerの疎多項式加減乗算・微分だけを使用した。
producer、SymPy、Fisher実装をimportせず、literal alternating行列と循環式から再構成した。

非対角上三角の十entryは

```
(-z3-e*z0, -z3, -z4, z2, z2, z3, -z1, -z1, z0, 0).
```

signed complementary Pfaffian vectorは正確に `-g`、かつ `Phi*(-g)=0`。
各 `Omega_ab` は残り三添字の三つの循環積
`(partial p)^t*(partial Phi)*(partial p)` を符号付きで加えて求めた。
producerに保存された十entryを単に読み戻したものではなく、再構成した全entryと比較した。

実coordinate mapは `z=M*x`、

```
M=[[5,-10,10,-5,1],[-4,12,-12,4,0],[4,-4,-4,4,0],
   [0,-16,16,0,0],[0,16,16,0,0]].
```

Fractionで逆行列を別に作り、正逆両積と `det(M)=16384` を確認した。
実tensor変換は `M^-1*Omega(M*x)*M^-T` である。
この全150係数は全eについて保存族と一致し、eに一次。
e=0では保存済みinteger canonical代表の正確に二倍で、Euler補正を必要としない。

さらに次を全eの多項式恒等式として確認した。

* canonical column divergence全五本は零。
* strict Jacobi全十本、Eulerを掛けたprojective Jacobi全350係数は零。
* 五つのsigned Pf minorをEuler収縮した実五次式の全126係数が一致。
* 実五次式は `det Jac_z(g)(M*x)` の正確に `-25/32768` 倍。
* e=0で五quadric spaceを結ぶ実basis matrixのdetは1073741824。
* 元五quadricのJacobian determinantは正確に `-2F`。

したがってtensorとその実Pfaffianを同時に認識している。
quintic又はzero curveの一致だけからtensor帰属を仮定していない。

## 2. Smooth-only Fisher定理と、nodal点への正しい特殊化

Weierstrass族は `y^2=x^3+x^2+e`。
generic `e*(4+27e)!=0` でsmoothで、infinityに唯一の滑らかな点を持つ。
`(1,x,y,x^2,x*y)` のpole ordersは0,2,3,4,5で、Riemann–Rochからdegree五の完全linear systemである。
その五quadricsがgである。
z0=1では元Weierstrass affine曲線に戻り、z4=1ではg1,g2,g4がz0,z1,z2を消去してg0,g3も含意する。
したがって余計なinfinity成分を持つ五quadricsの集合をsmoothモデルと呼んでいない。

e=0での全normalizationは

```
[s^4*t:s^3*t^2:s^2*t^3:s*t^4:t^5-s^5].
```

全五quadric pullback零と、Mを通した五section

```
(t-s)^5,
4*s*t*(t-s)^3,
4*s*t*(s+t)*(t-s)^2,
16*s^2*t^2*(t-s),
16*s^2*t^2*(s+t)
```

を別に全係数で確認した。正規化chartの有理分母はPGL mapに入っていない。
実mapは上のconstant GL5だけである。

既受理の `KNOWN_FIVE_SMOOTH_LOCUS_20261101.md` と独立reviewは、
Fisherの同じ循環式がsmooth elliptic normal quinticのq(5,1)を与えるというPolishchuk Theorem Aを入力としている。
今回この定理を使うのは上のsmooth generic eだけでよい。
Omegaはmodel係数の多項式であり、e又はdiscriminantで除していない。
e=0でのtensorは二倍の非零canonical代表なので、projective tensor族は0の近傍でも正則。
generic fiberが閉成分T39内なら、その中心もT39内である。
singularな曲線へsmooth-only定理を直接適用していない。
この実一次e族は全eでJacobiを満たしており、有限jetからの推測ではない。

## 3. Critical lambda=2/5 の全六源と、nil全fiber

旧rank零semisimple sourceのrawbaseは
`P02=q0,P03=2*q1,P13=q2`、他零。
そのdivergenceは `(-3*x0,-x1,x2,3*x3)`。
五変数のcanonical lift `rawP-E5 wedge div(rawP)/5` の負が新A_ss(2/5)baseに一致する。
cの一次項もcanonical条件から同時に一致する。

全六vertex spaceをcanonicalへ送るmapはinjective。
kernelならraw WはEulerの倍だが、各成分が(q0,q1,q2)に入り、integral coneはP3をspanするため、そのlinear倍は零。
新旧全kernelが六次元である既証明を用い、全spaceが一致する。
よって全六coefficientを保った旧rank零sourceへの実linear逆があり、既独立受理のwhole T38認識へ戻る。
一点でのdivisor一致や次元上界からT38を推測していない。

lambda=2/5と3/5それぞれで、nonzero binary quadratic Uのsemisimple・nil両actual orbitにおいて
元350×18全vertex systemのrankが十二、kernelが六である。
従って `C^3\{0}` 上の全kernelはrank六vector bundle。
baseもtotalspaceもirreducibleで、semisimple逆像は稠密。
全tensorへの係数mapは正則なので、全nil六fiberが対応するsemisimple全源の閉像に入る。
generic kernelの分母へnil代表を代入したのではない。

## 4. Lambda=3/5 のwhole sourceを代表へ戻す実逆

actualFの全六kernel係数は

```
(c+ell)*q0*q2 + a*q2*Cminus + b*q0*Cplus,
ell=(5*k4/2)*x0-5*k1*x1-(5*k2/4)*x2-(5*k3/2)*x3,
a=k0, b=k5/2.
```

四つのell方向と残二方向が独立であることは
`x0*Cminus=-q1^2-x2^2*q0`、`x3*Cplus=q1^2+x1^2*q2` と
相異なるirreducible quadricsのcoprimalityから従う。
よって固定baseの全vertex kernelからactualDへの写像はinjective。
実shear `c'=c+ell` はbaseとcanonical条件を保つため、Dのellを消した後tensor自体も一致する。

ab非零でactual Cartanとc-scaleの作用は
`(a,b)->(gamma*t^5*a,gamma*t^-5*b)`。
`t^10=b/a`、`gamma=(a*t^5)^-1` で全ab非零をa=b=1へ戻す。
complex十乗根の選択は有限で、反転するunitは明示したabだけ。
sourceを単一の標本に置き換えるのではなく、この実逆が全unit sourceを覆う。
§1–2のactual Fisher中心がこの代表なので、全unit sourceはT39へ入る。

semisimple全六fiberはirreducibleでab非零が稠密、T39は閉なのでa=0又はb=0も全て含む。
§3のkernel bundle橋からnil全六fiberもT39。
これらのzero係数の境界を、既約Dのopenだけで失っていない。

## 5. 残る全lambda・B・W・base零のcover

元のcomplete27sourceと五whole base normalforms、全十九wall、lambda無限遠を入力とする。
新しいsource探索を仮定しない。
Aのlambdaが2/5,3/5以外ならfull vertex kernelは四次元。
lambda=1/2以外では実shear四方向がactualDのc-linear項からinjectiveとなり、全kernelを尽くす。
全tensorをV=0へ戻すので、actualDはproper factorを持つ。
lambda=1/2も、二critical値を除いた全lambdaのproven constant-rank四bundle上で
proper-factor affine coneの閉性を使えば全fiberを含む。
有限minorの壁又はlambda=1/2の一点をgeneric basis代入で捨てていない。

W_ss/W_nilはbaseが `E wedge H` 又は `E wedge N`。
全vertex係数を保ったままprojective Eulerを引くとcommon `partial_c` wedgeなのでrank二以下。
base零も同様で、lambda無限遠はこの二源として別に含まれる。

B_nil全lambdaではcanonical Casimirと二linear fieldsの独立性から
`N(F)=0,H(F)=-10*lambda*F`。
degree五の全SL2分解のhighest linesは
`x0^5,x0^3*I,x0^2*J,x0*I^2,I*J,x0*Delta`。
対応するmodule dimensions16,12,10,8,6,4の和56が全空間を尽くす。
これらの相異Cartan weightsは全六eigenlinesを尽くし、全てproper productである。
従って全B、全jump wall、zero又は特殊係数にも被約既約quinticは無い。
特定のB parameterだけを有限標本で除いたのではない。

以上から、被約既約actual quinticを持つtwisted-cubic cone rank二の全sourceは
**T38又はT39**に入るというpointwise E定理を受理する。
Rのsquarefree・既約性、clean、normal node、vertex smoothnessを追加しない。
既存D100/M49、二次余核零、旧conditional global gapは変更せず、この陽性包含の前提にも用いない。
新component、smoothness又は全Poisson分類の完成は主張しない。
