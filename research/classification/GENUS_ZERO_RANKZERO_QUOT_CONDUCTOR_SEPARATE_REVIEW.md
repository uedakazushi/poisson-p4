# Degree五・種数零・rank零の全T40特殊化：独立したQuot／部分導手の監査

2026-11-08、PRIMARY E、別作者 `exceptional_coverage`。**純証明 PASS**。
対象は `GENUS_ZERO_RANKZERO_ACTUAL_T40_CLOSURE.md`、
SHA256 `cc6b209ce1fc38c517fab13f00a078f5e5052694314923c3ce30bd54ead706e8`。
旧全source還元、genus一認識、unstable bundleの実T40変形も合わせて読解した。
新しいCAS、point、operator、minor、rank、GB又はjet計算は無い。

受理する結論は、被約既約actual quinticを持つdegree五rank零支持のnormalization genus零全源が
閉成分T40に入ること。113型は実際にはrank四で空である。
degree五rank零のgenus一既済分岐と合わせる場合にも、rank二や他degreeの結論を流用しない。
新component数の増分は零である。

## 1. 全source入口と、任意の五sectionを保持するpinching

既受理の `RANK_ZERO_QUINTIC_SURFACE_GENUS_ZERO_REVIEW.md` は、全normalizationを尽くす。
cone三型、113のplane導手Q、generic cuspのdegree五支持は空。
残るものは `Y0=P_P1(Etilde)`、Etildeが113又は122、plane導手がp+q又は2pである。
全projection五空間Vは `W=H0(Etilde(-c))` の全二次元を含む。
`h0(Etilde)=8`、`h1(Etilde(-c))=0` なので、二fiber又は長さ二jetへの全評価は全射である。
これらの全源は一般位置のprojection、reduced plane導手又はplane上rank零に限定されない。

nodeでは二fiberが同じtarget planeへ同型に写る。
各評価kernelは二次元Wなので、指定されたV/Wから両fiberの実同型phiが定まる。
phiで貼ったmoduleは完備node環
`A={(f,g):f(0)=g(0)}` 上で
`{(u,v):v(0)=phi*u(0)}=A^3` である。
したがってgluingは任意のphiで局所自由。cyclic、semisimple又は一般phiを仮定しない。

cuspではV/Wのfirst-jet像はconstant termへの評価が同型であり、任意のmatrix Nに対して
`u1=N*u0` のgraphとなる。
完備cusp環 `A=C[[z^2,z^3]]` 上でこのlatticeは `(I+zN)*A^3`。
I+zNは完備正規化上unitなので自由rank三で、pullbackも元のbundleに一致する。
この局所計算は代数的なconductor squareの有限jet gluingの局所自由性を確認するためだけに使う。
formal変形又はformal target mapのままにしていない。

両場合ともgluing条件は独立三条件。
従ってrank三bundle E0は `H0(E0)=V`、h0五、h1零、degree五を持つ。
gluing点のfiber評価もV/Wの同型なので全域でglobally generatedである。
tautological bundleのfinite正規化pullbackはampleなLであり、finite surjective mapに対するample下降から
`O_P(E0)(1)` もample。
`H0(E0*)=0` は正規化への単射とEtilde*の負summandsから直接従う。

これで `Z0=P_B0(E0)` はintegral Gorenstein、有限正規化がY0であり、
元の五sectionによる有限双有理mapは `Y0 -> Z0 -> D0` と因子化する。
中心Z0がnormal又はsmoothであるとは置かない。
planeのnode gluing／cusp first jetはZ0の環に既に組み込まれている。

## 2. Relative Quotの障害零は、本当に全quotientをsmoothする

指定された五sectionとframeは、locally freeなkernel Q0を持つquotient
`O_B0^5 -> E0` そのものとして固定する。
B0はGorenstein genus一なのでSerre dualityにより

`Ext1(Q0,E0)=H1(Q0* tensor E0)=Hom(E0,Q0)*`。

Q0がO5のsubbundleなので後者はHom(E0,O5)へ単射し、H0(E0*)零から零となる。
singular curve上の非局所自由sheafへ誤ってvector bundleのdual式を適用していない。

B0をplane cubic smoothing B/Tの中心に置く。
特異点で非零の一般cubicを加えればtotal Bはregular、Tを0を残して縮めると他fiberはsmooth elliptic。
相対adjunctionは `omega_B/T=O` を与える。
relative Quotでは固定quotientの障害が上記Ext1であり、全small extensionの障害が消える。
局所的にはlocally free quotientの行列をliftでき、overlapの差はHom(Q0,E0)の一cocycle、
そのH1零によって貼れる、という同じ証明でもよい。
Quotのfinite presentationとinfinitesimal criterionにより、これはactual smooth morphismの近傍を与える。

smooth morphismのétale座標で相対座標を中心値に固定すれば、
0を通るétale-local algebraic曲線上のquotient `O_B^5 -> E` が得られる。
Tを必要ならこのregular曲線で置き換えるだけであり、有限jet又は未代数化formal arcではない。
propernessとlocally free／ampleの開性で、0を保って全域のEをlocally free、tautological bundleをrelatively ampleにできる。
cohomology/base changeと中心h1零からp_*Eはrank五。
指定の五sectionは中心でbasisなので、さらに縮めれば全fiberの完全frameになる。
ゆえに任意の元Vと元target frameを保持したままの変形である。

## 3. Imageを先にflat hypersurfaceとし、部分導手を正確にbase changeする

Z=P_B(E)からP4×Tへのmapはproperで、各fiberのtautological line bundleがampleだからquasi-finite、従ってfinite。
各fiberはintegral三次元、五sectionが独立なので像は非退化。
`L_t^3=5=deg(g_t)*deg(image_t)` と五の素数性から、各mapは双有理なintegral quintic imageである。

total image Dはregular P4×T内のhorizontal integral prime divisorなのでCartier。
base lineを局所trivializeしてdegree五の一本のFを取る。
F_tが零多項式ならverticalなambient fiber自体がDのdivisor成分になり、horizontal prime性に反する。
各F_tは同degree五のintegral imageを含むので、schemeとして正確にそのimageである。
これでD/Tのflatnessと `omega_D/T=O` が得られる。
正規化のflatnessや特異点型一定を先に仮定していない。

total BはregularでEはbundleなのでtotal Zはregular。
中心Z0は非normalでもCM/Gorensteinであり、Z/Tもflat Gorensteinである。
有限同次元mapの局所環A=O_Dに対しB=g_*O_ZはMCM。
実際、有限fiberと同じ局所dimensionによりAのparameter idealは各B局所環でparameter idealになり、
CM性から正則列となる。
AがGorensteinなので `Ext_A^i(B,A)=0` for i>0、finite dualityは `Hom_A(B,A)=g_*omega_Z/T`。

finite birational mapの共通fraction field内で、Homはevaluation at oneにより
`c_g={b in B:bB subset A}` に正確に一致する。
これはAかつBのidealであり、Z上のinvertibleなomega_Z/TをO_Zに入れる実partial-conductor mapである。
base parameter tauに対しAとBはflatなので、Ext1零から

`Hom_A(B,A)/tau = Hom_(A/tau)(B/tau B,A/tau)`。

evaluation at oneもこの同型と可換。
したがって中心で得るmapは正確にZ0->D0の導手mapであり、余分なtau因子又は全零specializationを持たない。
**Y0のfull normalization conductorとZ0のpartial conductorを混同しない。**
total Zがnormalでも中心Z0はnormalでないことを、このHom base changeが明示的に処理している。

## 4. Vertical Poissonと、元の中心tensorへの一致

相対dualizing generator `omega_B/T=O` を選ぶと
`omega_Z/T^-1 = omega_Z/B^-1 = Lambda2 T_Z/B`。
partial-conductor inclusionは従って全fiber上のregularなvertical bivector PZを与える。
Z/Bはsmooth relative dimension二なので、どのvertical bivectorもPoissonである。
この論証はbaseがnode/cuspでも成立し、dualizing generatorをKähler differentialと同一視しない。

PZの係数はpartial conductorに入り、そのB倍もAに入るので、
ambient座標の全bracketはAの正則関数となる。
よってdg(PZ)は `Lambda2 T_P4|D` の正則sectionへ下降する。
planeを越えるtensor matchingはZ0の環と真のpartial conductorで保証され、
正規化上だけの正則性から推測していない。

中心PZをY0へpullbackすると、smooth projective bundleのrelative tangentは任意base changeで同型だから正則。
degree五支持Sのgeneric preimageではZ0=Y0であり、partial conductorはCSに沿って消える。
従ってpullbackは全source空間

`H0(Lambda2 TY0(-CS)) = H0(Omega1_Y0(CP)) = C`

に入る。conductor mapとvertical top bivectorはgenericに非零なので、このsectionは非零。
元の任意の中心tensorも同じ一線にある。
dualizing generatorを一つの非零constantで調整して、中心を元tensorに一致させられる。
plane restriction rank二、cuspの全N、導手のcodimension二境界は失われない。

`H0(Omega2_P4)=H1(Omega2_P4)=0` から、下降sectionのambient bivectorへのliftは全familyで存在し一意。
中心liftは指定されたpi0である。
JacobiはDのdense smooth locusで零、従って全Dで零となり、
`H0(Lambda3 T_P4(-5))=H0(Omega1_P4)=0` によりambientでも零。
実Pfaffianはc(t)*F_tで、中心c(0)非零なので0を保ってcをunitにできる。
全fiberがactual rank四、被約既約quinticの実Poisson familyである。

## 5. Smooth elliptic側の全bundleをT40へ入れるscope

generic B_tはsmooth elliptic、E_tはrank三degree五、ample/GG、完全五sectionを持つ。
Z_tはsmoothなのでg_tは全正規化であり、PZ_tは**全導手**に沿って零。
したがってnormalization tensorの全空間は `H0(Omega1_Zt)=C`。
ここにintegral degree五の一つの導手surfaceという追加仮定は要らない。

stable Eについて既genus一認識は、既知T40のstable normalization bundleと、
Atiyahのcoprime rank三degree五 determinant一意性・translationによってpolarized projective bundleを同一視する。
この部分はbundle-wideであり、他の導手の分裂又はrank零支持の既約性に依存しない。
完全frameと上の一線のtensor一意性から、任意のstable Eのactual tensorがT40になる。

unstable ample/GG Eの全型は `A(deg2) directsum F(stable rank2 deg3)`。
これはample/GG quotientの必要degreeとelliptic bundleのsplitを用いる全型である。
前に独立受理した `UNSTABLE_ELLIPTIC_RANKZERO_T40_SEPARATE_REVIEW.md` は、
非split `Ext1(A,F)=C` の逆extensionから全導手／全frameを保持するambient実族を作る。
その証明は入口後はEと `H0(Omega1)=C` のみを使うので、ここでもそのまま適用する。
Quot smoothingのgeneric Eがstableだという未証明の仮定は不要。

以上から全t!=0のactual tensorがT40に入り、閉性により中心pi0もT40。
任意のnode/cusp gluingと全V・frameを保持する構成なので、全sourceに適用できる。

## 6. 113型の独立した空性と最終範囲

113では `H0(Etilde(-2))` はO(3) summandだけから来る二sectionであり、
planeを定める二ambient linear formsのratioはpullbackでs/tとなる。
normalization tensorがverticalなので、このratioはD上の非定数Casimir。
projective field `ell1*pi#dell0-ell0*pi#dell1` はTP4(2)のsectionであり、D上零。
`H0(TP4(-3))=0` によりambientでも零となる。
rank四のsymplectic一般点に非定数rational Casimirは無いので矛盾する。
node/cuspのどちらのcでも同じであり、113の全rank四sourceは空である。

以上を全源のT40包含と113空性として受理する。
固定vector bundleの閉包だけ、未確認cyclic normal form、normalization-flatness、
一部のsection又は有限jetによる論証ではない。
新component、rank二、他degree、全E又は全Poisson分類の完成を主張しない。
producer／旧証明書／main／global／Gitのbytesは変更していない。
