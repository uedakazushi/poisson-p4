# CA・CB原点の全完成型：conductor coverとintrinsic pushout

2026-11-03。CA・CBを**それぞれ独立に**同じ形式モデルへ移す。
これは全complete local ringのformal isomorphismであり、convergent analytic coordinate changeは主張しない。
CAとCBの新しい座標を互いに一致させるものでもない。relative embeddingとactual intersection Dは別問題である。

## 1. 実際のconductorと被覆の仮定

既に独立に受理済みの原点conductorはCM20・edim22・重複度4で、十九orbit因子と三枝curveの初期形を持つ。
三つのreduced枝はsmoothで独立なtransverse tangentを持ち、generic長は(1,1,2)である。
新しいbinary initial認識により、w=W、s=U、r=U−lambda*Vと置けば、初期idealは(wr,ws,sr²)となる。
lambdaはCAが1/3、CBが9である。

intrinsic normalization B、conductor quotient C、source conductor E=B/conductorは全てcanonicalなG-equivariant objectである。
smoothな全150係数ambientで十九formal group方向と補sliceを選ぶと、action mapはformal IFTでisomorphismとなる。
G不変なactual成分とconductorはこのisomorphismで十九因子を分離する。
有限被覆についても(g,e)→g*eと(g,e)→g^(-1)*eが互いに逆なので、normalizationとconductor coverを同時に分離できる。
Gの十九parameter sectionがsubgroupであることは不要である。
trace-zero lineはalgebra isomorphismで保たれるので、そのgeneratorの変更はunit倍だけであり、勝手なadditive correctionを含まない。

従ってabstract conductor/cover rigidityを、g19を除いたcurve coreへ適用できる。
その独立証明は `AMBIENT_THREE_AXIS_CONDUCTOR_RIGIDITY_20261103_v63` に固定する。
これはCM/no finite-length submoduleからnilradical²=0、nilmoduleがself枝DVR上free rank1、rsがそのgeneratorであることを導き、full conductorを

    C=C[[g19,w,r,s]]/(wr,ws,sr²)

へ移す。三次initialだけからhigher項を零に置いたものではない。

## 2. actual trace quadraticを分類する

q=nu mod Iはcritical W、Vの両simple枝で恒等的に零であり、q∈m_C²である。
actual odd二次jetにより、七actual regular odd constraintを解いた後のqの二次項はCAがUV/75、CBが16UV/9である。
これは保存されたnu-liftの二次値の負号を取り、actual IFT補正を適用したものとなる。
上の(r,s)で

    CA: q2=(1/25)*s(s−r),
    CB: q2=(16/81)*s(s−r)

であることを厳密QQで検査する。両係数とも非零である。

monomial conductorでqは一意にs²*alpha(s)+rs*beta(s)と書け、alpha(0)とbeta(0)はともに非零である。
formal square rootを取り、S=s*sqrt(alpha)、R=−r*beta/sqrt(alpha)と置けば、monomial relationsを保つformal automorphismでq=S(S−R)となる。
以後はこのS,Rを再びs,rと書く。
従って実source conductor coverは

    E=C[t]/(t²−s(s−r))

と正確に分類される。

## 3. regular normalizerへliftしてactual Aを復元する

actual regular normalizer Bは29次元、Eのembedding dimensionは23である。
Eのgenerator g19,w,r,s,tをBへliftし、kernelの独立な六linear coordinate z6を補う。
これらはBのregular system of parametersであり、minimal Cohen presentationの一意性から

    B=C[[g19,z6,w,r,s,t]],
    ker(B→E)=(z6,t²−s(s−r),wr,ws,sr²)

となる。
Cに入るEのgeneratorのliftは、intrinsic conductor pullback squareによりactual Aに入る。
z6もconductorとしてAに入り、t²のE像もq∈CなのでAに入る。

nu=t²、R=C[[g19,z6,nu,w,r,s]]と置く。
これは**新しく構成した**正則29次元subringであり、以前の共通R29と同じ座標であるとは仮定しない。
B=R directsum R*t、t²=nuである。

任意のfinite normalizationとconductorについてA=B×_E Cが厳密に成立する。
実際、E像がCに入るb∈Bはあるa∈Aとの差がconductorに入り、その差もAに属する。
このintrinsic pushoutにより、全actual completionは

    A=R directsum I*t,
    I=(z6,F=nu−s(s−r),wr,ws,sr²)

と復元される。ここでhigher formal項を捨ててはいない。

## 4. 完全99関係

b=(z1,...,z6,F,wr,ws,sr²)とし、十odd targetをY_i=b_i*tと置く。
全kernelは次の99関係で生成される。

* 七regular generatorの21 Koszul syzygy。
* 七regular generatorと三core generatorの21 cross syzygy。
* coreの二syzygy s*Y_wr−r*Y_ws、rs*Y_wr−w*Y_sr²。
* 55個のpair products Y_i Y_j−nu*b_i*b_j。

最初の44は全第一syzygyである。
nuをF+s(s−r)へ変えれば七regular sequenceとcore Hilbert–Burch resolutionのtensor productが全resolutionを与える。
pair productsで任意の式をR+sum R*Y_iへ還元し、BのR-basis(1,t)で残る関係がこの44 syzygyで全て生成されることが分かる。
従ってpullback零だけでなく、これが全kernelである。GB・Macaulay2は実行しない。

conductorのcore primary decompositionは

    (wr,ws,sr²)=(r,s) intersect (w,s) intersect (w,r²)

であり、generic長(1,1,2)を明示する。
critical quotientはU=s、V=r−sにより正確に(WU,WV,UV,U³)となるため、既知の全source critical modelとも一致する。

coneのactual imageはdim29/edim39/depth21/e2、projectiveではg18となりdim28/edim38/depth20/e2である。
本定理が分類するのは指定B18原点での個別CA・CBの全形式完成環であり、全global bad locus、相対CA/CB配置、actual intersection D、full Poisson schemeの完成環ではない。
