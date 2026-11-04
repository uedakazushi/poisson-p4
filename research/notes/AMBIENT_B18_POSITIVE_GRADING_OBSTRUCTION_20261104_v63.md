# Generic B18 の群由来 positive grading が使えない理由

2026-11-04。基点は受理済みの実際の四項

    P0=e_y∧(x²∂h+hx∂z+(h²+z²)∂w), (h,x,y,z,w)順

である。x²∂z∧∂wという別の平方constant-plane tensorではない。対象は元の全非被約D=CA∩CBの接平面39を群軌道19で割った20次元空間である。Dred、N normalizer、R20の接空間へ置き換えない。

## 1. Diagonal stabilizer と零weight class

GLの座標weightを(H,X,Y,Z,W)、テンソルscalar weightをσとする。P0を固定する四つの条件は

    H+Y−2X+σ=0, Y+Z−H−X+σ=0,
    Y+W−2H+σ=0, Y+W−2Z+σ=0。

従ってH=X=Z=W=r、Y=r−σである。σ=0のpureGLだけなら中心scalarとなるが、受理済み作用は **GL5×tensorGm** である。非自明な固定torus

    gλ=diag(1,1,λ⁻¹,1,1), tensor multiplier λ

が存在する。quadratic monomialに含まれるyの個数をkとすると、physical coefficientのweightはbivector pairがyを含む場合k、含まない場合k+1。全て非負である。

しかしu63=P^{xy}のhz係数はweight0であり、全actual D接平面上で非零、全orbit19上で零である。したがって群因子を除いた実cotangentにも非零のweight0 classを与える。独立受理されたnil-cotangent基底の一つでもあるが、生のu63そのもののnilpotenceやその高次liftをここで計算したものではない。

## 2. 任意の群由来torusへの拡張

保存された同じ150 physical coefficientから全25個のGL微分を線型演算で作り、以前のorbit19を全係数で再現する。pureGL stabilizer Lie algebraは6次元で中心Iを含む。tracefree部分は5次元である。その一般線型結合M(t0,…,t4)について、tr(M^k)=0(k=1,…,5)、さらにM^5=0を全係数で確認した。よってtracefree部分の全要素はnilpotent matrixである。

H=(Stab_GL(P0)∩SL5)^0と置く。そのLie algebraはこの5次元tracefree部分であり、Engelの定理で同時にstrictly upper triangular化できる。connected complex algebraic subgroupのLie algebraがこのnilpotent algebraなので、Hはunipotentである。pureGL stabilizerの任意のtorusをSL5へ写す際のdeterminant方向は中心GL scalarであり、中心は全quadratic bivector表示に自明に作用する。従ってeffective pureGL stabilizerに非自明なtorusはない。

S=Stab_(GL5×tensorGm)(P0)とする。torusのtensorGmへの射のkernelはpureGL stabilizer内にあり、effective表示では有限である。したがってeffective Sのtorus rankは高々1であり、上の非自明torusによって正確に1である。cocharacterの像はS^0に入り、connected linear algebraic groupのmaximal toriはS^0内で共役である。

SはDをスキームとして保ち、T_Dとorbit tangentを保つ。従ってSによる共役はquotient T_D/T_orbit のtorus weight multisetを保つ。任意の非自明effective cocharacterは上のrank1 torusの非零整数冪と共役であり、weight0の重複度は消えない。trivial cocharacterなら全weight0である。よって **群作用から来る、P0を固定するどのC*でも、全20 transverse target coordinateをstrictly positive weightにできない**。

これは単に選んだdiagonal frameでの失敗ではない。一方、GL5×tensorGmに由来しない抽象的なformal local-ring actionや、零weight parameterを残す相対weighted手法は排除していない。元D全nil・socle・depth19/20の判定も行わない。ここで閉じたのは「群由来の全方向positive gradingによって有限weighted IFTへ進む」というmethod gateだけである。

全計算は線型weight・小さい5×5 polynomial matrix演算だけで、GB、新しい点、高次source jet、primary solverは使わない。
