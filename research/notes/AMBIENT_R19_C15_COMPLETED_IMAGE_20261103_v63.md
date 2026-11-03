# R19のC15軌道における全完成像環

2026年11月3日。これは被約R自身の局所環であり、全Poisson空間や実CA/CB scheme交叉の局所環ではない。
既存のproper spectral模型と全normal locus分類を使う。

## 定理

C15の任意の点で、被約Rの完成局所環は正確に

    A = C[[u1,...,u15,B,E,F,p²,p³,Bp,Ep,Fp]]
      ⊂ Bnorm = C[[u1,...,u15,B,E,F,p]]

となる。右辺はnormalizationである。
次元19、embedding dimension23、depth16、multiplicity2で、analytic domainである。
normalizationの一点fiberは長さ2である。

normalizer内のconductorは

    (B,E,F,p²) Bnorm

であり、target A内ではそのcontraction、すなわち

    (B,E,F,p²,p³,Bp,Ep,Fp) A

である。この二つのambientを区別する。
normalization quotientはC[[u15]]*pである。
この点で被約境界B18は滑らかで、完成環はC[[u15,B,E,F]]となる。

## 1. 四変数の有限spectral族

C15の代表点は

    P0=∂y∧D0,
    D0(h)=x², D0(x)=0, D0(z)=hx, D0(w)=h²−xz/6.

normal spectral sourceのchart d=b=1、e=0で、四変数(p,B,E,F)を次の通り置く。

    t=−p/12, r=11p/12, c=−2p²/3,
    G=h²+(B−p/6)hz+(E−1/6)xz+Fz².

quadric r²+c−25t²=0は恒等的に満たされる。
源のcanonical tensorは

    D(h)=x*((r−13t)h+x+c*z),
    D(x)=12t*x²,
    D(z)=x*(h+(−r−13t)z),
    D(w)=2t*xw+G.

これに有限GL変換

    g=I+p*(E_zw+E_xh),  g^−1=I−p*(E_zw+E_xh)

を作用させる。行列のnilpotent部分の二乗は0である。
全150 canonical係数は有限多項式であり、元のtensor・quadric・canonical恒等式を厳密に検査した。

target係数から中心値を引いたものをv_iとする。indexは既存150係数順序で0始まりである。
次が全て**厳密な恒等式**であり、jetの切り捨てではない。

    B = v123+v112,
    E = v127,
    F = v132,
    p² = v15,
    p³ = −3*v18/2,
    Bp = v108+v127*v15+v15/6,
    Ep = v112,
    Fp = v117.

逆に全150係数の各monomialは、この八生成元の環に入る。
奇数p冪はp³以上か、B,E,Fの少なくとも一つを持つことを全項で検査した。
従ってこの四源族の**全像代数**は正確に

    A0=C[[B,E,F,p²,p³,Bp,Ep,Fp]]

である。

## 2. full R像であること：十九次元source chartと二十三次元代数cotangent

指定点でprojective orbitの次元は15である。
projective chartをtarget係数u106=1とし、全GL25作用から独立な15方向を取る。
それらを積として実現する局所15次元GL chartをHとする。

四源族のB,E,F微分とこの15 orbit微分はtargetでrank18を持つ。
p微分はtargetでは0だが、sourceの旗は

    X: x−p*h,
    K: span(h,x−p*h,z−p*w)

と実際に動く。従ってp方向はsourceで非零である。
全19 source微分を、physical150係数と旗のPluecker係数を合わせて検査しrank19を得た。
normal sourceはこの点で滑らか19次元なので

    H15 × (p,B,E,F)-slice → normal source

はetaleで、完成源環を同定する。
全proper fiberが一つであることは前稿で示してあるため、このnormal sourceはtargetの近傍の有限normalizationである。
従ってRの完成局所環は、その源の完成正則環への注入として扱える。

ここで単に四源族を横断sliceと仮定していない。
Hの作用を含む全target係数はA0[[u15]]に入るため、target完成像はこの環の部分環となる。
一方、A0の八つのminimal algebra cotangent生成元

    B,E,F,p²,p³,Bp,Ep,Fp

の係数ベクトルと、15 projective orbit列を合わせた行列は厳密にrank23で、
そのminorは−209/162である。
したがってambient完成環からA0[[u15]]への写像はcotangent上全射である。
complete local algebraのNakayama（逐次近似による全生成）により、この写像は全射となる。
よって**実際の全target完成像がA0[[u15]]に等しい**。
これはtangentだけでlocal ringを推測する議論ではなく、既に全有限出力が入るcomplete代数への全射を証明する議論である。

## 3. 全有限presentation、conductor、depth

D=C[[u15,a,B,E,F]]、a=p²と置き、j=(a,B,E,F)とする。
AはD-moduleとして

    A = D ⊕ (a,B,E,F)*p

である。z0=p³,z1=Bp,z2=Ep,z3=Fpと書くと、全kernelは十六生成元

    j_i*z_k−j_k*z_i       (0≤i<k≤3),
    z_i*z_k−a*j_i*j_k     (0≤i≤k≤3)

である。第二群で全式をD+ΣD*z_iへ還元でき、残るkernelはregular sequence jのfirst syzygiesである。
それらは第一群のKoszul syzygiesで尽くされるので、このpresentationは全kernelである。
全十六式は有限源への代入で厳密に0となることも検査した。

normalizerはD⊕D*pで、差は

    Bnorm/A = (D/(a,B,E,F))*p = C[[u15]]*p.

よってconductorの二つのambientでの表示を得る。
正確な列0→A→Bnorm→C[[u15]]*p→0によりdepth A=16となる。
minimal generatorは八つなのでembedding dimensionは15+8=23である。
maximal idealのnormalizerへの延長は

    (u1,...,u15,B,E,F,p²)

で、regular normalizerでのmultiplicityは2である。
有限birational拡大はrank1なのでAのmultiplicityも2となる。
normalization fiberの残余環はC[p]/(p²)であり長さ2である。

最後に全sourceのboundary Cartier divisor t=0はこのchartでp=0である。
H作用はそのdivisorを保つので、被約B18の完成像はp=0への制限、すなわちC[[u15,B,E,F]]である。

## 証明書と未解決部分

* `code/ambient_R19_C15_completed_image_20261103_v63.py`
* `reports/AMBIENT_R19_C15_COMPLETED_IMAGE_20261103_v63/report.json`

全QQ計算0.70秒。有限150出力、八つの厳密なtarget生成式、23-rank minor、
旗を含む19-rank source minor、全十六presentation式を保存した。
B16上のquadric-vertex像、B17のdiscriminant locusの全完成環、実CA/CB交叉のboundary nilpotent structureは別問題である。
