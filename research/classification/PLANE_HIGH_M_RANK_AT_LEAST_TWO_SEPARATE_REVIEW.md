# rank零plane・全m・k一・rank M二以上の非優越：独立査読

2026-11-07、主義務G。producerと別作者による全源幾何・Fraction算術の独立確認。**PASS**。

**受理する結論。** 実退化因子Dが既約被約quinticであるrank四Poisson構造について、
rank零plane、k=1、rank M≥2の**全重複度**のprojective PGL像は次元24以下である。
rank M二の高重複度m≥3では全像23以下、さらにM(p)=0の全高重複度源は空である。
rank M三の旧除外τ零ではsemisimple型24以下、nilpotent型23以下となる。
既知のm二結果と合わせた全源のG非優越であり、成分帰属又は全E/Gの尽くしではない。

## 1. 固定入力、有限job、独立算術

旧不変packetは `inputs/plane_rank_zero_rank_two_high_m_manifest_20261107.json`、
新κopen packetは `inputs/plane_rank_zero_kappa_open_manifest_20261107.json`。
各manifestの全file/hash、旧rank三及びnonkernel report、全m modular補題と独立reviewを、
`jobs/plane_high_m_rank_at_least_two_review.json` に固定した。
旧初回reportの `all_G24_closed=false` と旧flag reportの
`all_three_whole_flag_exclusions_pass=false` を変更せず、今回もfalseであることを検査した。
それらは、旧linear bound又は定数rowspace criterion単独では全体を閉じなかったという正しい記録である。

独立checkerは `code/classification_plane_high_m_rank_at_least_two_review.py`、証明書は
`reports/CLASSIFICATION_PLANE_HIGH_M_RANK_AT_LEAST_TWO_REVIEW_20261107.json`。
180秒・1536MiBの契約をgate PASS後に実行し、実測5.682秒、peak RSS32,280KiBで全検査PASS。
gateは数学の証明ではない。
既存の**別作者独立reviewer**のFraction(τ)、疎多項式、Schouten、Euler商だけを再利用し、
producer、producer helper、SymPy又はCASをimportしない。
全operatorと全非線形Jacobiは、保存された行列を信じてrankだけを計算したのでなく、
全150係数と全basisから独立に作り、その後に全entryを照合した。

## 2. 全m入口とrank二の全high-m first symbol

全m local modular補題 `REDUCED_DEGENERACY_ZERO_SURFACE_ABELIAN_CONORMAL.md` と独立reviewにより、
被約Dのrank零plane上ではconormalが可換である。
従って全weight零35、一30、二9とco-Higgs共通線O(0)/O(1)入口を使える。
旧m二で除外したleading Pfaffian零のrowを、この入口で復活させる。

k一の非零constant vector U、p=[U]、n=(z,w)、En=z∂z+w∂wを用いると全weight零源は

`Pi0=T∧En+U∧M(x)n+qS(n)∂z∧∂w`。

M0=M(U)、R=S+M0/2、T tracefreeとすれば全Jacobiの必要十分なLax式は

`T(M)−tM+[R,M]=0`, `[T,U]=−tU`, `[R,M0]=0`。

この式はmに依存しない。normal二次Pfaffianは非零定数を除き
`qM(n)*det(Ebase,T,U)`。
rank M二ではqM非零だから、m≥3は後者のdeterminant恒等零と同値である。
二次元W=im Mの全normal共役型は、Killing orthogonal lineの非等方・等方に対応する
`WN=<F,E>` と `WB=<H,E>` の二つで尽くされる。

### M(p)非零

旧nonkernel全Lax源のs,ell,τを残してdeterminantを課すと、NDのM0 semisimple/nilpotentではs=ell=0、
Borel semisimpleでは4τ−3s=0,ell=0、Borel nilpotentではs=ell=0となる。
独立係数は順にsa(b−c),ell(b−c)²、sac,ell c²、(4τ−3s)ac,ell c²、sab,ell b²である。
parameterの特別値で一方の係数を消して別のrowを見落とす操作は無い。

全四源は、base(a,b,c)で次の通りとなる。

|型|U|M|T|K=qS|
|---|---|---|---|---|
|NDSS0|∂b+∂c|bF+cE|0|−(z²−w²)/2|
|NDN0|∂b|bF+cE|0|−z²/2|
|BSS0(τ)|∂b|bH+cE|diag(−2τ,4τ,−2τ)/3| (1−2τ)zw|
|BN0(τ)|∂c|bH+cE|2τb∂c|(1/2−τ)w²|

τは全Cである。元のwhole M、Uの可逆normalizationは旧nonkernel証明の全逆と零補集合を保ち、
今回新たに消すhigher parameterは無い。

### M(p)零

U=∂a、M=bM1+cM2と置く。全Lax解は
`R∈N(W), J=adR|W, T(a)=ta+rb b+rc c, T(b,c)=(tI−J)(b,c),3t=trJ`。
determinant零はtI−Jがscalarであることを意味する。
NDではR=αH、J eigenvalues±2αだからα零。
BorelではR=αH+βE、J eigenvalues0,2α、off-diagonal−2βだからα=β=0。
したがって全high-mは正確に

`R=S=0, T=(rb b+rc c)∂a, K=0`。

R零では旧row shearがrowを全く動かさないため、rb,rc両方を保持した。
NDではnormal diag(t,t⁻¹)、base b'=t⁻²b,c'=t²cが
rowを(t²rb,t⁻²rc)へ送り、有限swapを併用して非零row全体を(1,ρ)、ρ=rb rcへ直す。
row零は別型で、ρ零の非零axisも含む。
Borelではnormal exp(sE)diag(t,t⁻¹)、base b'=b,c'=−2sb+t²cが
rowを(rb+2s rc/t²,rc/t²)へ送る。
rc非零はt²=rc,s=−rb/2で(0,1)へ直り、rc零では(ρ,0)の全ρを保持する。
全変換はactual invertibleで逆も同じ式の逆、必要な平方根はC上有限choiceである。
従ってND(ρ)、ND row零、B rc非零、B rowline(ρ)の四全源を得る。

## 3. 全30/9 operatorと全wallの再構成

全30方向とEuler六方向を合わせたraw36のrank、全九のweight二方向を独立に検証した。
旧packetの全150basisと一致した。
七つの新Pi0を上記T,U,M,Kから作り、全350 projective Jacobi零と実Pfaffian零を確認した。
各々全350×30と350×9を再構成し、保存全entry、complete kernelの全行代入零・独立性、nonzero maximal minorを検査した。
variable源では選択minor積とτを掛けた例外多項式の全一次分解を検証し、全rational wallをspecializeした。
非有理の未検査wallは無い。

|M(p)非零whole源|affine symbol B|generic nullity和d|全wall後の像上界|
|---|---:|---:|---:|
|NDSS0|10|8|23|
|NDN0|9|8|22|
|BSS0(τ)|10|7|22|
|BN0(τ)|9|7|23|

BSS0の全wallは0,−3/2,1/6,1/2,3/2,−1/2。
順にd=8,7,7,8,7,7である。
BN0の全wallは1,0,−3/2,1/2,−1/2,3/2で、d=9,8,9,10,8,10。
NDSS0は旧NDSSのτ零tensorそのもので、旧excluded arithmeticを全350から今回も再構成し、d=8を再認証した。
affine Bは全(p,M)十、Borel一減、M0 nilpotentなら更に一減、自由τ一増から直接数えた値である。
fixed wallではBを一減し、plane六とtensor射影化−一を戻す。
任意P1でP2は同じD2のinhomogeneous fiberだから全higher dはkernel和以下であり、全四源の上界が従う。

R0-ND非零rowの唯一のwallはρ零。
Borel rowlineの全wallは0,−3,1/3,1,5/3,−1/3,3,−1。
全wallのoperator/minor/kernelを同じく独立に認証した。
巨大なlinear上界だけではR0源をG非優越とは呼ばない。

## 4. R0の全common-U及び二次因子の排除

ND全ρとrow零、Borel rc非零、Borel rowlineのgeneric openと0,−3,1/3,−1/3では、
完全D1 kernelの全tensorを再構成するとPi1=∂a∧V1、係数はaに依らない。
全D2 kernelも六次元で、全て∂aを持ちa-independentである。
Pi0も同じ形なのでstrict Jac(Pi1)=0、weight二のforcingは零。
従って任意の実Pi2はker D2へ入り、全Piは共通∂a因子を持ちrank≤2となる。
高weight非線形解を未証明のまま商にした結論ではない。

Borel ρ=−1では全D1 kernelのうち∂aを持たない方向はdb∧dcだけである。
全九Pi2もbase pairだけなので、その全未知係数を自由にしても
`Pi=∂a∧V+F2∂b∧∂c`。
実Euler Pfaffianはquadratic F2で割れ、F2零ならrank≤2、非零なら既約quinticに反する。
ここではPi2をlinear kernelと仮定していない。
残るのはρ=1,5/3,3だけとなる。

## 5. 三つの全非線形源、実Pfaffianとκopen

三wallのcomplete D1 kernelは14,14,15方向。
この全basisを前節の独立検証済みkernelから再構成し、旧flag reportの全150basisと一致した。
全Pi2九方向を加えた全源から、parameterを独立変数として全strict Schouten及びEuler wedgeを展開し、
全350のlinear/quadratic coefficient polynomialを再構成した。
保存された全350行列と全entryで一致する。
さらに正しいEuler収縮Pfaffianを同じ全源から直接作り、w零の全56spatial coefficientと
全parameter coefficientを照合した。

保存された全original350-row functionalを各々別に評価した。
z5以外の全nonzero modw coefficientはliteral Jacobi row combination、
z5 coefficientのremainderは正確に−drだった。
dはPi1の{a,w}のz²係数、rはPi2の{b,c}のz²係数であることも、全150係数から直接照合した。
従って任意の全Jacobi solution上で

`D|w0=κz⁵`, `κ=−dr`。

κ零ではwが非零実quinticを割り既約性に反するため、全対象はd r非零に入る。
この補集合を未処理として捨てていない。

actual diagonal PGL変更とtensor scalar

`(a',b',c',z',w')=(ra/d³,b/d,c/d,z,w)`, `sigma=d²/r`

はPi0を保ち、d'=r'=1へ送る。
全150係数の(d,r) weightsを独立に計算した。
coordinatesはa=(−3,1),b=c=(−1,0),z=w=(0,0)、scalar=(2,−1)。
Pi0の全weight零、全P1 basisのhomogeneous weightは(−1,0)又は(1,−1)、
全P2は(−2,0)六方向と(0,−1)三方向である。
従って全他係数もLaurent monomialで可逆に移り、逆はa=d³a'/r,b=db',c=dc'とtensor inverse r/d²。
全open源はGm²と全normalized源の積へ同型であり、d=r1だけの点計算ではない。

## 6. 三つのliteral Jacobi矛盾と全normalized maps

全350 rowの零始まり番号51,121,122は、順に
((a,b,c,z),bz³)、((a,b,c,w),bz³)、((a,b,c,w),bz²w)である。
A=J51,B=J121,C=J122とし、公開P1 kernel順でd=p8、m=p13、
ρ1,5/3ではell=p11とする。
rは全parameter番号20,20,21のP2係数である。
全他係数を零にせず、独立に得た全多項式は次の通り。

|ρ|A|B|C|
|---|---|---|---|
|1|−6r+(2p0−p3)ell+4dm|−3d ell|−4r+(2p0−2p10)ell−4dm|
|5/3|−8r+(8p0/3−4p3/3)ell+8dm/3|−11d ell/3|−6r+(8p0/3−p3/3−8p10/3)ell−14dm/3|
|3|−12r−4dm|0|−10r+6dm|

ρ1,5/3ではd非零とB零からell零。
その後それぞれA+C=−10r、7A+4C=−80rとなりr非零に矛盾する。
ρ3では直接3A+2C=−56rで矛盾する。
全三源が空である。

補助的にproducerの全normalized affine elimination certificateも独立に評価した。
各段階のcurrent350-row functionalが保存affine equationへ一致すること、
各solved coefficientの全forward map、自由係数のinverse projection、
全350 residual equation、元の全350へ全mapを合成した恒等式を照合した。
最終のnonzero constant consequenceは正確に1。
ρ1,5/3は二段、ρ3は一段であり、次数は二以下、未知数の個数で先験的に有限である。
literalな三行証明だけでも十分で、GB、jet、ideal-multiplierの追加探索を必要としない。

## 7. rank三のτ零と全mの正確な結論

rank M三ではM可逆だからM(p)非零、全Lax解はR=τM0で、旧semisimple/nilpotent二型だけである。
leading Dはそれぞれ
`−4acτ(az²−2bzw−cw²)`、`c²τ(az²−2bzw−cw²)`。
従ってm≥3は正確にτ零である。
この二点のPi0をU,M,Kから別に構成し、全350×30/9 operator、nonzero minorsとcomplete kernelsを再認証した。
両方ともD1 kernel八、D2 injective。
全affine symbol次元はτ固定でsemisimple十一、nilpotent十なので、全像上界は
`11+8+6−1=24`、`10+8+6−1=23`となる。
旧m二 theoremの除外scopeを遡って書き換えてはいない。

rank二high-mはM(p)非零の全四源で≤23、M(p)零は全て空。
rank三high-mは≤24。
既受理のm二・k一・rank M二以上≤23を合わせ、**全m・k一・rank M二以上の全像≤24**が従う。
未知G成分の必要次元25〜27に反するので、この全branchは非優越である。
k零、rank M一以下、rank二restrictionのplane、全E membership及び全Gの尽くしは残す。
成分数を変更せず、909主TeX/PDF、global台帳又はGitは変更していない。
