# 個別CA・CBのformal bad loci：独立監査

2026-11-03。全形式モデルとseminormalizationの受理済み独立証拠を前提に、
`AMBIENT_ORIGIN_FORMAL_BAD_LOCI_20261103_v63` を再生し、別のmonomial計算と局所環論で監査した。
元のsource・ノート・flat reportは変更しない。

## 全primeでのlocal AB argument

A=R⊕I*t、A_sn=R⊕J*tの既知の模型を用いる。
p⊃Iの上の**fiber A tensor_R kappa(p)** はkappa(p)とsquare-zero odd部分の和である。
従ってその上のAのprimeは一つ、residue fieldもkappa(p)で、局所化は有限local R_p-module R_p⊕I_p*tとなる。
原ノートの「A/pAが体kappa(p)とodd部分」という箇所は、このfiberの意味で読む必要がある。
非極大primeの場合、未局所化A/pAそのものを体上のalgebraとは置かない。

R_pは正則で、A_Pと同じdimensionを持ち、pA_PのradicalがA_Pの極大idealである。
CMならlocal cohomologyのdepthはR_pについてもfullとなり、ABによりA_Pはfreeである。
そのrank1直和因子I_pはprincipalとなるが、Iの三minimal primeは全てheight9なので、どのp⊃IでもI_pのheightは9である。
principal ideal theoremと矛盾する。Jもheight9のprimeで、同じ議論が成立する。
閉点でdepthが2未満だと述べる議論ではなく、CM条件そのものを否定している。

conductor外では環はregular Bと一致し、その上ではB/Aが非零なので非正規である。
従ってSing=Nonnormal=NonCMはそれぞれのconductor支持に正確に一致する。

## Nonseminormal支持の独立確認

C=R/Iのcore monomial idealはK=(wr,ws,sr²)である。
monomial divisibilityを使い、全monomialに対する厳密なcolonを計算した。

    K:w=(r,s),     K:r=(w,rs).

前者ではwrまたはwsがwmを割ればrまたはsがmを割り、sr²が割ってもsがmを割る。
後者ではwrまたはwsがrmを割ればwがmを割り、sr²が割ればsrがmを割る。
逆の生成元包含も直ちに検査できる。有限次数までのsample検査ではない。
LCMによる交わりは(wr,ws,rs)=sqrt(K)であり、C内ではAnn_C(w,r)=(rs)=Nil(C)となる。

A_sn/A=(w,r)C*tは従ってconductorの全被約支持を持つ。
各局所化でproperな有限subintegral extensionが残るので、Aはその全支持でnonseminormalである。
外ではnormalだからseminormalであり、A_sn自身は既に全normalizer内のsquare/cube testでseminormalと受理されている。
よってoriginalの四bad lociは全三枝、A_snの前三bad lociはself枝だけ、A_snのNonseminormal locusは空である。

## 正規化fiberと数値

p⊃Iではodd elementのB内の像がpBに入るので、normalization fiberは正確にkappa(p)[t]/(t²−nu_bar)となる。
W、V、Gではnu_bar=0でdual numbers長2、self枝でs!=0ならnu_bar=s²で異なる二根±sを持つ。
A_snではself枝以外は既にBなのでnormalizationはidentity、self枝では同じ二次式を使う。
normalization fiberと、local AB argumentのA→Spec R fiberを混同しない。

各stratumのregular base dimensionとideal generator数から、独立に

* original W/V: (dim,edim,depth)=(9,18,1)、selfも(9,18,1)、Gは(10,20,2)。
* seminormal W/V: 正則(9,9,9)、selfは(9,18,1)、Gは(10,19,2)。

を再計算した。proper conductor点ではeven linear項をodd productsが消さず、edim=dim(R_p)+mu(I_p)となる。
depthはR_p⊕I_pと0→I_p→R_p→R_p/I_pから得られる。
各proper conductor点でodd generatorのsquareはp²に入るのでpA_Pは極大idealのreductionであり、rank2からmultiplicity2となる。
その外はregularでmultiplicity1である。

G一般点の横断数値に十九smooth因子を戻すと、original coneは29/39/21/e2、seminormal coneは29/38/21/e2となる。
projectiveでは十九因子が十八になり、それぞれ28/38/20/e2、28/37/20/e2となる。
既にsmooth方向を除いた横断9/18表からさらに1を引くものではない。
三support枝はcone20・projective19、共通intersectionはcone19・projective18である。

結論は各個別のformal completionでの被約閉支持と局所不変量だけである。
global bad locus、収束解析型、relative CA/CB embedding、actual intersection Dには拡張しない。
