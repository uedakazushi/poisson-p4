# m2・linear root rows・irreducible conic determinant：全source被覆の別作者検証

2026-11-07、独立 reviewer `/root/quartic_review`。**PASS（純証明）**。
producer note `PLANE_M2_LINEAR_ROOTS_CONIC_MPS.md`、SHA256
`e84908de976f6902db8081029d134aca2a8fd9a4f8d3f4f00ec7c172e4277339`
を全文読解した。root-line入口と既済 normal filtration/MPS全component scope を入力とする。
新しい線形代数、Jacobi、minor、GB、jet計算は行っていない。

対象は restriction rank2 plane、m2、generic squarefree normal quadratic、
被約既約actual quintic D、reduced非零 cubic f、degree flag(1,1)で
`G2=κh ℓm, h detM=f, detM` irreducible conic の**全source**である。
結論は、hが両root-row basepointsを通らなければ空、両方を通る残全sourceはMPS40。
split determinant、他flags、double normal quadratic、nonreducedR、m3等へ拡張しない。

## marked rowsを保持する全正規化

linear Mのimage WはMat2内の三次元hyperplane。
二次元以下ならdetは二変数以下のbinary quadraticで複素数上分解するからである。
annihilator Kがrank1ならdet|Wも分解する。
このrank判定自体には任意の行列基底を使ってよいが、source正規化にはleft GL2を使わない。
Kがrank2なのでnormal右GL2だけでK→Id、base GL3で

```
M=[[x,y],[z,-x]],  Q1=xu+yv, Q2=zu-xv, q=-x²-yz
```

に戻せる。root rowsを混ぜた偽のorbit除算は無い。
両rowsはprimitiveなので各Qiはrank4 quadric threefold。
verticesはp1=[0:0:1]、p2=[0:1:0]。hの全positionを保持する。
両点を通るhはxの非零倍だけ。root尺度・normal/base逆変換・非零振幅κは可逆なまま残る。

## closed residueの全域pole評価と二ruling-line関係

全normal filtration Π=Π0+Π1+Π2 のΠ0はPoisson、actualD0=κhQ1Q2≠0。
よって逆ω0はclosed rational2formで、三つの被約成分だけに高々simple pole。
Qiのgeneric residue αiはclosed rational1formである。
Qi∩QjはS=(u=v=0)とintegral cubic scroll Tのreduced和。
Tは `(-a st,a s²,-a t²,b s,b t)` の像で、両componentはgenericにtransverse、重複度1。
特にS/T/h-section genericでαiのpoleはsimple/log、residueは各irreducible component上constant。
頂点はQi内codimension3、他のbad curvesはcodimension2であり、以後の一般直線はこれらを避ける。

rank4 quadricの二ruling-plane families内で一般のP1直線を取る。
Sのintersection degreesは(1,0)、Qi上のQj-cutはdegree(2,2)、従ってTは(1,2)、Hは(1,1)。
これらはvertexを避けた一般直線上の実交点数であり、singular resolution上のpole boundを仮定しない。
h(pi)≠0ならH∩Qiはsmooth irreducible quadric、h(pi)=0なら二つのopposite ruling planesに分かれる。
base-only Hのためrank3-sectionという余分な場合は無い。

H∩Q1がirreducibleとするとα1のP1への制限の留数和は

```
rS+rT+rH=0,  2rT+rH=0.
```

従ってrS=rT=λ、rH=−2λ。
差 `α1−λ dlog(Q2/h²)|Q1` はclosed simple-pole形式で全divisorial residue零。
よってQ1のsmooth locus上regular。rational coefficientにcodimension2だけの新しいpoleは無い。
全一般ruling-line P1上regular1formは零、二ruling-plane tangent spacesの和はTQ1なので差自体が零。
一般のreflexive differential延長定理や頂点blowupのexceptional極を使っていない。

α2のS,T residuesはantisymmetryにより−λ。
H∩Q2がirreducibleなら同じ議論、split E+FならS+T=2H=2(E+F)と二ruling-line関係から
rE=rF=2λとなる。zero-residue差を同様に消すと

```
α1=λ dlog(Q2/h²),   α2=−λ dlog(Q1/h²).
```

従って少なくとも一方のH-sectionがirreducibleなら
`ωlog=λ dlog(Q1/h²)∧dlog(Q2/h²)` がω0と全Qi-residueを共有する。
差はP4上でhだけに高々simple poleを持つclosed2form。
h≅P3上のresidueはregular1formなので零、残るP4上regular2formも零。
ゆえにω0=ωlogはrank≤2となり、rank4 inverseという仮定と矛盾。

この全open-position空性は有限個のgeneric sampleからの推論ではない。
hがp1又はp2の少なくとも一方を避ける**全position**を排除する。
残るのは両方を通るh=xだけである。
tangent conic lineはconic上の相異なるp1,p2両方を通れないので、本flagのtangent case全体は空。

## 残h=xの全higher sourceからactual SNCへ

残Π0のdivisorは `κx(xu+yv)(zu-xv)`。
xは被約Poisson divisorなのでprojective hyperplane ideal(x)はPoisson-invariant。
raw x-row modxはEuler倍ρEだが、base Jacobian代表f=xqのx-rowは既にxで割れる。
従ってmodxでρy=ρz=0、ρ=0となり、mixed normal entriesも含めraw x-rowがxで割れる。
任意raw代表にdivergence条件を無断に課した議論ではない。

N≥2、weights(x,y,z,u,v)=(1,0,0,N,N)とする。
Π0のbase-base、base-normal、normal-normal全degree2 monomialを保持すると重みは非負。
xを含むpairだけは上のx因子を使い負重みを消す。
weight1の**raw36方向全部**が重み≥N−1（Euler6のquotient30も同じ）、weight2全9方向が≥2N−1。
したがって全Π1/Π2を任意に戻しても、

```
Πt_ij(X)=t^(-wi-wj) Πij(t^w X)
```

はtの多項式、t≠0で元とPGL同値、t0でもPoisson。
actualDは `Dt=t^(-Σw)D(t^wX)`、Σw=2N+1。

```
x(xu+yv)(zu-xv)=x²zu²+xyzuv−x³uv−x²yv²
```

の唯一最小relative weight0項はxyzuv。
normal order r≥3の全項はrelative weight≥(r−2)N−1>0である。
よって `D(Π0limit)=κxyzuv`、κ≠0、真のSNC log-symplectic pointを得る。
period・resonance・magnetic又はhigher係数の一部を固定していない。

各元源に対するこのactual irreducible曲線全体を含むambient componentを一つ取れば、
それはSNC点に会う。既受理MPS Theorem1.4の**SNCに会う全component** scopeによりMPS40。
元源もこの既知閉成分に入る。
一点のknown membershipから曲線全体を推論する逆境界shortcutでも、T01の一意性でもない。

したがってconic determinantのflag(1,1) whole sourceは空又は全MPSで尽くされる。
残るsplit determinantと他flagsの補集合を消さず、全rank2plane/全E/Gの完了は主張しない。
