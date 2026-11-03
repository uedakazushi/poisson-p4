# R19境界局所模型の独立監査

2026年11月3日。対象は凍結した `ambient_R19_boundary_local_models_20261103_v63.py` と
`AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_v63.md`。
被約R、被約B18のgeneric B17での二枝完成環と、全Rの

    Sing(R)=Nonnormal(R)=NonCM(R)=B17∪C15

を、以下の計算と全fiberの論証の範囲で受理する。
全Poisson局所環や実際のCA∩CBの境界schemeは対象ではない。
新しいC15 pinched完成環の別packetも、この監査の対象に含めない。

## 独立な全係数計算

producerのvector/action helperをimportせず、五変数の二次単項式の150係数から直接構成した。
全25個のLie作用は、二ベクトル方向の置換と−x_v∂_uを整数の疎行列として実装した。
十二次元spectral basisの全150係数が既存証明書と一致する。

指定generic B17点での二つのR branch、二つのB18 branch、その二つずつのunion、
B17族の接行列、およびC15のsource接行列と軌道行列の、全九行列を独立に再構築した。
保存された全成分と一致し、QQ上の階数と指定minorも一致する。
順にaffine階数20,20,22,19,19,20,18,19,16である。

全A,B,E,Fを残したone-Casimir chartでも、十九次minorが定数48であることを検算した。
その列を消去した全残差は、同じ一列の

    B, A/6+E, F, −B, −A/6−E, −2F

だけである。有限個の点の階数から全例外を推論していない。
任意のa,b≠0、任意のcと全六G係数を持つchainの有限座標変換も独立に計算し、
正規化とhx,x²項の除去を全係数で確認した。

## generic B17での全幾何fiber

指定点のcanonical fieldは

    D=(0,0,−hx,h²+2hx+3hz+4x²+5xz+6z²)

である。intrinsic factorの一意性を確認した。
定ベクトルvがU∧v∧D=0を満たすとし、x=0へ制限すると、Gの非零z²項からvはw方向に限られる。
次にD(z)=−hxを使うとv=0となる。従って元の五変数でのfactor lineはC∂yだけである。
線形Casimirはh,xの二次元空間Wだけである。

既存のproper spectral模型のt≠0像はR\Nに入るので、この点の全持上げはt=0上である。
sourceのXはW内の直線となる。
Xがbinary quadratic hxを割らなければ、D(ell)がXで割れる線形形式の空間はWの二次元だけであり、
三平面Kを作れない。Gの非零z²係数は、w成分をこの空間へ加えることも許さない。
従ってX=ChまたはCxだけであり、どちらでもK=span(h,x,z)に一意に決まる。
固定旗上の線形写像は単射なので、全幾何fiberは厳密に二点である。

propernessと有限fiberによりtargetを縮めると射は有限となり、normal birational sourceはnormalizationとなる。
両source点はquadricのsmooth点である。
projective接階数19により各点でunramifiedで、各fiber点の長さは1となる。
有限な完成normalizationの各直和因子はNakayamaでtarget完成環の商となり、各branchは正則19次元である。
単にsourceがnormalなことからtargetのnormalityを推論してはいない。

B17自身のproper triangular sourceでは、LとWに加え、Kがker(D²)から一意に戻る。
実際、D²(w)=−hx*(3h+5x+12z)は非零で、他のh,x,zはD²で零となる。
全fiberが一つでprojective接階数17なので、同じ有限・Nakayama論法でB17はここでsmooth17である。

binary hxの二根は相異なるため、formalまたはétale局所で両根を選べ、
両R branchは同じB17をschemeとして含む。
両接空間の交わりは17次元であり、scheme交叉はsmooth B17を含むので、
その局所次元と埋込み次元が17に等しく、交叉自体が正則17次元でB17と等しい。
この次元と接次元の議論が、支持の一致だけからscheme交叉を推論することを防ぐ。

被約excellent Rの完成は被約であり、normalizationの二点に対応する二branchを持つ。
二つのminimal primeの交わりが零なので、完成環は両正則branchの共通B17上のfiber productとなる。
従って

    C[[u17,a,b,c,d]]/(ac,ad,bc,bd)

という実際の完成環を得る。
同じ議論をt=0のsource Cartier divisorへ適用すると、B18の各branchは正則18次元であり、
共通17次元B17に沿うnode C[[u17,s,t]]/(st)となる。
fiber productの正確な列から、Rのdepth18、dimension19、multiplicity2、
normalization quotient C[[u17]]とconductor(a,b,c,d)が従う。

## B18\B17全点、特にC15での有限normalization

既存の全chain被覆でB17の外ならa,b≠0であり、Gはspan(hx,x²)の外である。
a,bのみを反転する有限座標変換で、全点は

    D=(x²,0,hx,A h²+B hz+E xz+F z²),
    (A,B,E,F)≠0

へ移る。
線形Casimirはxだけであり、sourceのXは一意である。
G mod x≠0なら、D(ell)がxで割れる空間がK=span(h,x,z)そのものとなる。
G mod x=0ならE≠0で、D=xN、Nは四次元regular nilpotentである。
その唯一の三次元invariant spaceはker(N³)なので、Kはやはり一意である。

intrinsic Lも一意である。前者ではx=0制限によりgeneric B17と同じ論法が使える。
後者でU∧v∧NU=0なら、U=v+s uの係数比較からNv∈Cvとなり、商U/Cv上のNは全直線を保つ。
従ってその商作用はscalarで、rank(N−lambda I)≤1となるが、regular nilpotent4のrank3と矛盾する。
このedge caseを除外していない。

以上から**全**B18\B17で、全geometric spectral fiberは一つの旗に尽きる。
sourceのquadricはd=1 chartでsmoothであり、proper有限化後の射は有限normalizationである。
C15についても、この有限性は接階数の主張に先立って成立する。

定数48と全残差の計算から、projective微分は
(B,F,A/6+E)≠0でinjectiveとなる。
唯一のfiber点がunramifiedならfiber長さ1であり、Nakayamaからtargetとの局所同型、従ってsmoothnessを得る。
残るB=F=A/6+E=0では、B17外なのでA≠0であり、wのscaleで全点が一つのC15代表へ移る。
軌道次元は16−1=15である。
ここではnormalization微分がinjectiveでないため有限normalizationは同型でなく、Rはnonnormalとなる。
この監査はC15の具体的なcusp/pinch完成環までは主張しない。

## 全locus等式と依存範囲

全R\Nのsmoothnessは、AF≠0だけの古いgeneric-double receiptではなく、
`AMBIENT_R_COMPLEMENT_20261103_v63_REVIEW.json` の、全六AF=0代表を含む監査を入力に使う。
proper quadric、全B18被覆、全RのR1は既存のspectral-model監査を入力に保持する。

generic B17のnonnormalityとnonnormal locusの閉性から、既約B17の全点がnonnormalとなる。
全RはR1であるため、そこで局所環がCMなら、その全localizationがS2を満たし、Serre基準でnormalとなって矛盾する。
従って全B17とC15はnonCMでもある。
これは各閉点のdepthが2未満という主張ではない。実際generic B17閉点のdepthは18である。
残る全点は前節と全complementの定理でsmoothだから、冒頭の三locusは厳密に一致する。

実装は `code/ambient_R19_boundary_local_review_20261103_v63.py`、
報告は `reports/AMBIENT_R19_BOUNDARY_LOCAL_MODELS_20261103_v63_REVIEW.json`。
