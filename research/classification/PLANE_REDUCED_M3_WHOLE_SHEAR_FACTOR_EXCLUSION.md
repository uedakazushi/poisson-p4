# reduced nonconcurrent R・m三：残二whole源の実因子

2026-11-08、`/root/quartic_review`。新しいfirst symbol、点、normal jet、GBは使わない。
既存の全九源のうちM及びSF_C3_UNUSEDについて、**全ての実退化五次式は可約又は零**と証明する。
したがって被約既約quinticを仮定した両枝は空である。
これはG上界22/23からの帰属推測ではない。rootの別作者読解は本稿作成時には未完了。

## 1. 入力と全源の実座標変更

座標を(a,b,c,u,v)とする。既独立算術
`PLANE_REDUCED_M3_NINE_SEPARATE_REVIEW.md`の完全K1、全9 P2及び全350式を入力とする。
P1はraw36/Euler6の既定30基底、P2はbase pair三つと(u²,uv,v²)の全9基底である。
K1の列順は元reportのconstrained_weight1_certificate.complete_kernelに固定する。
以下のα_iはその完全kernel座標で、chosen部分空間の座標ではない。

六個の実PGL shearを

    a'=a+s0 u+s1 v,
    b'=b+s2 u+s3 v,
    c'=c+s4 u+s5 v,
    u'=u, v'=v

とする。逆は各sを負にするだけで、det=1、全s∈Cで正則。
full homogeneous tensorには通常のpush-forwardを行い、その後weight一Euler項だけを除く。
raw36への全6 Euler列と全30基底を合わせた36×36の定数逆を用いた。
shearのweight一作用はP0だけに依存する。P1の変換による追加項はweight二で、
全9 P2に入る。P0の変換によるweight二項も同様である。
従ってP1上の下記の定数逆が全域で存在すれば、全tensorの逆変換でも全P2は9のまま保持される。
Euler項を除く操作は実射影Poisson tensor及び実Euler-Pfaffianを変えない。

Mでは

    s0=α3/3, s1=(α5−2α4)/3,
    s2=2α4/3, s3=−α0/3,
    s4=(α0+2α6)/3, s5=−2α1/3,
    r=α2.

SF_C3_UNUSEDでは

    s0=α0/2, s1=α2, s2=−α1, s3=−α3/2,
    s4=−α7, s5=α4, r=α5, s=α6.

これらは「six shear dimensionsがある」という数え上げではなく、
元の全7/8kernel係数と(shear,residual)との定数線形同型そのものである。
full30×6 shear matrix、7/8×7/8 matrix、その逆及び非零定数detはreportに保存した。
αのどの値も除かず、r/s/leading amplitudeの零も含めて逆shearで次の源へ移せる。
shearがleading parallel spaceを保つのは、実D3のassociated-graded base-degree二部分を変えないためでもある。

## 2. Mの全残源

上記逆shear後、全P1はr倍の次の一方向だけである。

    ΔPac=−b v, ΔPau=v², ΔPbc=(b u+c v)/2.

これに**全9 P2を未知量として加えてから**元350 Jacobiを展開した。
既存J2の9×9非零定数minor det=531441/4で全9未知量が零と強制される。
この代入後、残る元350式は全て恒等零である。全normal quadratic bracketも含む全代表は

    Pab=2ac−b²,    Pac=2bc−r b v,
    Pau=2bu+r v², Pav=−3au/2−bv,
    Pbc=c²+r(bu+cv)/2,
    Pbu=cu, Pbv=−cv/2, Pcu=0, Pcv=3cu/2, Puv=3u²/2.

実五次式は全5 signed PfaffianをEulerで合わせた式であり、直接展開すると

    D = (r/2)(bu+cv)(3au²+6buv+3cv²+r v³).                 (M)

r=0はD=0、r≠0では二つの正次数因子を持つ。座標shear又はEuler regaugeでこの性質は失われない。
したがってMの**whole** irreducible-quintic locusは空である。
必要ならr≠0ではbu+cvはrank4 quadric、後因子はa一次でprimitiveなので絶対既約cubicである。
既受理のquadric+cubic全支持定理へ接続できるが、本稿の空性にはその外部結果を用いない。

## 3. SF_C3_UNUSEDの全残源

逆shear後の全P1はΔPau=r v²、ΔPbv=s u²だけ。
全9 P2に対し元J2のrank8定数minorを反転すると、8係数は零、
唯一Pabへのq uvだけが自由に残る。affine particularとhomogeneous kernelを両方保持している。
全代表は

    Pab=ab+q uv, Pac=−ac, Pbc=bc,
    Pau=2au/3+r v², Pav=av/3,
    Pbu=−bu/3, Pbv=−2bv/3+s u²,
    Pcu=−cu/3, Pcv=cv/3, Puv=−uv/3.

代入後の元350式の非零係数は−2rs,2rs,4rsだけである。
従って正確なPoisson条件はrs=0。いずれのbranchも、qは任意のまま保持する。
Poisson条件を課す前から実因子恒等式

    D = −c[2as u³−2br v³+(q+rs)u²v²]                     (SF)

が成り立つ。したがって全解のDはcで割れ、既約quinticにはならない。
さらにこのprofileのsquarefree cubic leadingはr,s両方非零を必要とするので、
指定SF leading locusそのものもrs=0と両立しない。
後者の強化を使わず、proper actual factorだけでも目的の全空性は従う。

## 4. 算術receiptと範囲

実行前job `jobs/plane_reduced_m3_whole_shear_source.json` はgate PASS。
新しいfocused runner `code/classification_reduced_m3_whole_source.py` はproducerをimportせず、
既別作者の疎Fraction基底を使い、係数環QQ上で実push、定数逆、全350式及び全Pfを再構成した。
symbolic係数の展開と定数線形代数にはSymPyのexact rational polynomialを用いた。
有限体、補間、数値rank、GB及び新しいfirst-symbol operatorは用いていない。

一回の全二源jobは2.033163秒、最大59,316KiB。
report `CLASSIFICATION_PLANE_REDUCED_M3_WHOLE_SHEAR_SOURCE_20261108.json` は
完全実行source、input hash、全kernel逆、全9P2解、元350全残差及び実因子を保持する。
MはP1の一点でなく全r、SFは全(r,s,q)とその全compatibilityを扱った。
既report/source/jobは改変していない。

この結論を全root-cover尽くしと合わせると、reduced nonconcurrent零cubic・normal位数三を持つ
rank2 planeが存在する**既約quintic**の残二枝が空になる。
concurrent、nonreduced R、他のnormal位数又は全Poisson分類の完了は主張しない。
