# degree五rank零支持・sectional genus一：canonical枝排除とelliptic scroll全源

2026-11-07、主義務G。新producer幾何証明。
外部入力の正確なscope、全normalization-sourceの完全性、既知T40への橋を分けて記す。
独立review前はglobal完了主張へ混ぜない。

**命題。** πをP4上のrank四Poisson tensorとし、実Dを既約被約quinticとする。
degree五integral surface S上でπが零となり、Dの正規化(X,L)のsectional genusが一と仮定する。
この全分岐は、既知T40=q_{5,2}の像、又はdimension≤24の射影像で覆われる。
S又はXのsmoothness、導手のnormality、ordinary-node性又はvertical divisor不存在を仮定しない。

## 1. 長さ・全normal-scopeのpolarized birational分類

入口のlength公式からδ_S=1で、S以外に二次元の非正規支持は無い。
正規化ν:X→D、L=ν*O_D(1)はnormal polarized threefoldで、Lはample・globally generated、L³=5。
`RANK_ZERO_SEXTIC_CONDUCTOR_EXCLUSION.md` §1と同じ二hyperplane/Bertini議論により
g(X,L)=1である。resolutionρ:Y→X上のρ*Lも同じsectional genusを持つ。

直接読解した一次資料はY. Fukuma,
*A lower bound for the second sectional geometric genus of quasi-polarized manifolds and its applications*,
Rend. Sem. Mat. Univ. Politec. Torino 69(1) (2011), 73–90,
<https://seminariomatematico.polito.it/rendiconti/69-1/73.pdf>、Proposition3(iii)、印刷76頁。
smooth quasi-polarized manifoldのg=1の全場合について、polarizationを保つbirational同値の下で
次の二つだけを与える。

1. normal Gorenstein・Q-factorial・terminal Vとnef big Cartier H、`K_V+2H∼0`。
2. smooth elliptic curve B上のP2-bundle Vと、fiberにO(1)をrestrictするnef big H。

同値の意味は同論文Definition1、印刷75頁のcommon-resolution上でpolarization pullbackが等しいこと。
Xがnormalなのでcommon resolutionのpushforward OはOであり、全mについてsection ringは等しい。
元のLがglobally generatedなのでHもglobally generatedである。
実際、common resolution上でpullback sectionがある点に非零なら、そのVへのimageでも同じsectionが非零。
ample model `Proj⊕H⁰(mH)` はXに等しく、V→XをfとすればH=f*Lである。
したがって「birational分類」をXそのものの二型と誤認していない。

## 2. canonical型のample modelはcanonical Gorenstein Fano

第一型ではpushforwardして `K_X∼−2L`。
従ってK_XはCartierである。
compatible canonical divisorを選ぶと `K_V−f*K_X` はf-exceptionalでnumerically trivial。
negativity lemmaの両符号適用により零となる。
fはcrepant、VはterminalなのでXはcanonical、特にkltである。
元々のXがterminal又はQ-factorialとは主張しない。

`K_V=−2H` とH nef bigからKawamata–Viehweg vanishingでH¹(O_V)=0。
terminal singularityはrationalなのでsmooth common resolution ZでもH¹(O_Z)=0、
従ってH⁰(Ω¹_Z)=0である。

もう一つの直接読解した一次資料はD. Greb–S. Kebekus–S. J. Kovács–T. Peternell,
*Differential forms on log canonical spaces*, arXiv:1003.2913v4,
<https://arxiv.org/pdf/1003.2913v4>、Theorem1.4、印刷3頁。
klt complex varietyのsmooth locus上のp-formは任意のlog resolutionへregularに延長する。
Xはcanonicalゆえkltなのでp=1を適用し、

`H⁰(X,Ω_X^[1])=H⁰(Z,Ω¹_Z)=0`。                     (1)

terminal-onlyのdel Pezzo分類をcanonical Xへ無証明で拡張してはいない。
degree五の特異V5分類又はtwisted-reflexive-cotangent消滅も使わない。

## 3. 閉residue一形式とδ一の全二型

Dのsmooth locusでπをrestrictしたrational bivectorをXへpullbackしてπ_Xとする。
generic transverse curveのδが一なので、全一般型はnode又はcuspだけである。
nodeならreduced conductor C_Xに沿ってπ_Xは零。
これはnode完成環の二branch商へのsurjectionとπ|S=0から従い、
`RANK_ZERO_SEXTIC_CONDUCTOR_EXCLUSION.md` の独立review済み局所証明を使う。

cuspならC_X=2Γ、Γは一つのprime divisorである。
generic完成環 `K[[t²,t³]]⊂K[[t]]` で、π|S=0とSのPoisson invarianceから
base同士のbracketは(t²,t³)へ入り、
`{a,t}={a,t²}/(2t)` はtで割れる。
従ってπ_XはΓに沿って少なくとも一次で零となる。
ここでK上のderivationを零とは仮定せず、base同士とmixed bracketの両方を保持する。

πの逆symplectic二形式σはreduced Dに沿って高々simple poleを持ち、Poissonによりdσ=0。
Dのsmooth locusでのPoincaré residueはclosed一形式α_Dである。
dualizing trivialization Ω_Dを選ぶと、四次元のinverse-Pfaffian恒等式から

`α_D=nonzero_constant * contraction(π_D,Ω_D)`。

Ω_DのXへのrational pullback Ω_Xはdivisor −C_Xを持つ。
したがってα=contraction(π_X,Ω_X)はclosed rational一形式で、
nodeならdivisorial poleを持たず、cuspなら唯一のΓに沿って高々simple poleを持つ。

cuspで閉性を使うと、Γのsmooth一般点では
`α=λ*dt/t+regular`、λはΓ上のconstantとなる。
一般の二ample divisorで切ったsmooth projective curve Aを、SingXと全codimension二bad locusを避けて選ぶ。
AはΓを五点（multiplicity込み）で横切り、α|Aの全residueの和は
`λ*(Γ·L²)=5λ` である。
compact curveのresidue theoremによりλ=0。
従ってαはcuspの場合もdivisorial poleを持たず、normal X上でreflexive holomorphic一形式となる。

(1)によりα=0。しかしDがπの実被約退化因子なのでπ_Dは一般にrank二で、αは非零。
これは矛盾する。よって第一型の全canonical ample modelは空である。
このsimple-pole結論はδ一のnormalization parameter計算に基づき、
vertical divisorを含む一般のblowup residueへsimple-pole仮定を課したものではない。

## 4. elliptic scrollのnef・globally generated bundleを全て保持

第二型をV=P_B(E)、H=O_V(1)と書く。Eはrank三、degree五、globally generatedである。
まずdegree零quotientを捨てずに処理する。
非零E→O_Bのimageはglobally generatedなO_B(−Z)なのでZ=0であり、
Eのglobal generationはquotientに1を与えるsectionを持つからO_Bがsplitする。
これを繰り返し `E=O_B^r⊕F`、H⁰(F*)=0と出来る。

このFはampleである。
もしtautological morphism P(F)にcontracted curve Tがあれば、fiber上のO(1)はampleなのでT→Bはsurjective。
Tのnormalization上でdegree零のglobally generated quotient f*F→Mがあり、Mはtrivial。
global evaluationからO_T^{h⁰(F)}→O_Tを得て、その係数はconstantなfunctional λとなる。
λはTの全点でevaluation kernelを消し、T→BがsurjectiveなのでBの全点でもkernelを消す。
従ってglobal F→O_Bが存在してH⁰(F*)=0に反する。
tautological morphismはfiniteとなり、O(1)はampleである。

r>0ならdegree一sectionのimage DはO_B^r方向をvertexとするconeである。
これは全section ringがFのsection ringにr個のdegree一変数を加えた形になることからも分かる。
元のν:X→Dはfinite birationalであり、唯一のdegree五非正規支持Sもconeとなる。
高次数cone排除（又はdegenerate支持補題）に反する。
したがってr=0、Eそのものがampleであり、H ampleゆえample model XはVそのものである。

全elliptic分岐はsmooth X=P_B(E)、E ample・globally generated、rank三degree五へ到達した。
PicXのL係数について `C_X=−K_X=3L−p*detE` は奇数。
cuspならC_X=2Γとなり矛盾するので、generic typeはnodeである。
したがって

`π_X∈H⁰(Λ²T_X(−C_X))≅H⁰(Ω¹_X)=C`。              (2)

各polarized bundleに対してprojective π_Xは高々一つ。
L ampleからH¹(E)=H⁰(E*)=0、よってh⁰(E)=degE=5。
νはその完全linear systemなので、projection subspaceを自由に選ぶparameterは無い。
また `H⁰(P4,Λ²T_P4(−5))=0` により、π_Xからambient πも高々一つである。

## 5. unstableな全bundleは一つのHN型で、像dimension≤24

ample・globally generatedなrank rのquotient QはH¹(Q)=0、h⁰(Q)=degQである。
global generationはdegQ≥r+1を強制する。
もしh⁰(Q)=rならevaluation O^r→Qがisomorphismとなり、ample性に反するためである。

Eがunstableなら、rank二destabilizing subbundleのdegree≥4はline quotient degree≤1を作るので不可能。
rank一destabilizerのdegree≥2を取ると、rank二quotientのdegreeは三以上なので、
厳密にdegree二line bundle Aとrank二degree三globally generated ample quotient Fになる。
Fがunstableならそのline quotientのdegree≤1となり矛盾するのでFはstable。
`Ext¹(F,A)=H¹(F*⊗A)=0` は、dual bundle F⊗A^−1がstable negative degreeであることから従う。
従って

`E=A_2⊕F_(2,3)`                                  (3)

が全unstable型であり、extension parameterを失っていない。
`h⁰(EndE)=1+1+h⁰(F*⊗A)=3` なので、actual polarized automorphism subgroup
Aut(E)/C*はdimension二。

Bのmoduliはdimension一、A∈Pic²Bは一、stable Fのmoduliは一以下で、全parameterは三以下。
Fの上界は `h¹(EndF)=h⁰(EndF)=1` というcurve上のdeformation計算でも得られる。
base translationsはこのbundle parameterをdimension一で動かす。
実際detEのdegree五なのでtranslation stabilizerは有限。
全target frame dimension24と(2)のprojective source dimension零を合わせ、
bundle automorphism二方向とtranslation一方向は同じambient πを与える実fiberである。
finite birational complete-L mapのtarget frameを固定してidentityとなるautomorphismだけがkernelで、これはidentity。
従って全unstable像は

`24+3−2−1=24`

以下となる。

## 6. stable bundleの全像はT40と一致する

ここで用いる追加一次資料はM. F. Atiyah,
*Vector bundles over an elliptic curve*, Proc. London Math. Soc. (3) 7 (1957), 414–452,
<https://math.berkeley.edu/~nadler/atiyah.classification.pdf>、Theorem7及びCorollary、印刷434頁である。
互いに素なrank・degreeのindecomposable bundleはdeterminantにより一意に定まる。
stable bundleはindecomposableなので、rank三degree五にこの正確な部分を適用する。

各elliptic curve Bについて、既受理のMPS/q_{5,2}構成は、stable rank二degree五bundle上のT40 tensorを与える。
その退化D_0の正規化は、Pic²(B)上のsmooth P2-bundle X_0=P(E_0)で、
L_0=ν*O_D0(1)はample・globally generated、rankE_0=3、degE_0=5である。
この全normalization/conductorの橋は `IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE.md` §4で既に証明・一次照合済み。
Pic²(B)はBとisomorphicなelliptic curveなので、baseをBと同定出来る。

まずE_0はstableである。
もしunstableなら(3)でAut(P(E_0),L_0)のdimensionが少なくとも二となる。
νの完全五sectionはこれらをP4のactual projective automorphismへfaithfully写し、D_0を保存する。
一方D_0の非正規支持は、stable rank二degree五bundle W上のsmooth elliptic scroll S_0=P(W)である。
そのpolarized automorphism groupはdimension零。
kernelはAut(W)/C*=0であり、base translationはdetWを保存するため五torsionだけとなる。
elliptic curveのorigin-preserving automorphismも有限である。
projective automorphismがD_0を保存すればS_0を保存し、S_0がP4をspanするので
S_0上でidentityならambientでもidentityとなる。
従ってAut_PGL5(D_0)もdimension零であり、二次元subgroupは不可能。E_0はstableである。

任意のstable rank三degree五Eについて、translation t_bを選んで
`det(t_b*E_0)=detE` と出来る。degree五line bundleへのtranslation作用はmultiplication-by-five isogenyでsurjectiveだからである。
Atiyahの上記一意性から `E≅t_b*E_0`。
従ってpolarized P(E)は既知T40のpolarized normalizationと同型であり、
その完全五sectionによるimage DもPGL5を除いて同じである。
(2)の一意なprojective bivectorとambient restrictionのinjectivityにより、π自身もT40と一致する。

既知q_{5,2}のrank二bundleと新rank三bundleを、dual-spanだけから同一のprojective targetだと仮定していない。
ここではnormalization、finite polarized automorphisms、determinantの一意性という実幾何の橋を使った。

## 7. 到達点

Fukumaによる全normalization-sourceの二分岐を保持し、canonical型は空、
elliptic型のunstable全源はdimension≤24、stable全源はT40とした。
これはg=1、degree五rank零支持の全G排除である。
g=0は別noteの全degree五+plane≤24、degree六は別noteの全排除で処理する。
rank二restriction、plane、quadric surface及び全G/Eの残余を本命題から閉じたとはしない。
