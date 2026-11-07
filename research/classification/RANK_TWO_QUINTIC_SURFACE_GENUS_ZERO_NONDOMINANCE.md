# 次数五rank二支持・正規化genus零：導手の分裂と全射影上界

2026年11月7日。rootによる新しい幾何証明案。別作者監査前はglobal台帳へ昇格しない。

**命題。** rank四Poisson構造の実Dが既約被約quintic、次数五の非正規曲面成分S上で
restrictionがrank二、Dの正規化のsectional genusが零なら、全PGL射影像は次元24以下である。
Sのnormality、clean又はvertical不存在は仮定しない。

## 1. rank二のnodal曲面上では二導手sheetが分裂する

一般横断型がnodeの曲面Sを考える。正規化導手のS上部分は、二つのbirational sheet、
又は一つのprime divisorからの二次被覆のどちらかである。
inverse symplectic二形式のclosed residueを正規化へ引き戻した一形式をαとする。
nodeの全branch商でPoisson bracketはregularであり、αはこの導手に沿って高々simple poleを持つ。
closednessから、各prime divisor上のresidueは複素定数である。

Sの一般点の分離拡大・完成後、D=(uv=0)、S=(u=v=0)と取る。
dualizing generatorの二branchへの制限は、それぞれ
`g*dv/v∧da∧db` と `-g*du/u∧da∧db` となる。同じambient単元gを使っている。
π|SのcoefficientをHとするとαの二residueは `gH` と `-gH`。
rank二restrictionから両方とも非零である。

もしS上の導手が一primeなら、その二次function-field拡大のinvolutionは二residueを交換する。
一方、そのprime上のresidueは複素定数λなのでinvolutionはλを固定する。
上の反対符号からλ=-λとなり、標数零でλ=0。これはrank二に反する。
従って必ず二prime Γ1,Γ2に分裂し、それぞれSへbirationalに写る。
係数体をPoisson中心と仮定した議論ではなく、global residueの定数性とbranch交換を使う。
他の導手成分又はcodimension二の特殊点はこの一般点の結論を変えない。

## 2. degree五・genus零の全normalizationとplane導手

既済のlength式からδS=1、他に唯一plane P、δP=1がある。
Höringの全normal g0定理とminimal-degree分類から、Dの正規化はdegree五scroll五型。
005ではDはline vertexを持ち、二次元非正規支持はplaneだけなのでdegree五Sは存在しない。
残る四型に対し、Y=P_P1(E)を標準scroll model、L=O_Y(1)、fiber divisorをFとする。
Eの型は113、122、014、023。
前三桁に零が無ければY=X、零が一つならY→Xはzero-summand section P1だけを収縮するsmall morphism。
exceptional divisorは無い。

L³=5、L²F=1、全導手のstrict transformはC=3L-3Fである。
実際、Dのdualizing generatorをYへrationalに引き戻すと、全codimension一でdivisorは-C。
small exceptional curveは追加Weil divisorを作らない。
S上部分C_SのL²-degreeは10、plane上C_Pのdegreeは二。

effective divisor aL+bFのdegree二条件は5a+b=2、a≥0。
H0(Sym^aE(2-5a))を全aについて調べると次表で尽くされる。

| E | 全C_Pの可能性 |
|---|---|
| 113 | 2F、又は唯一Q=L-3F |
| 122 | 2Fのみ |
| 014 | 2F、Q+F_b、2Q。Q=L-4Fは唯一、F_bは任意fiber |
| 023 | 2F、又は唯一Q=L-3F |

最大summand degreeはmax(E)*a+2-5aであり、113/122/023ではa≥2は不可能。
014ではa=2だけ追加され、その唯一sectionはhighest-coordinate平方、すなわち2Q。
a=1の014では全sectionはhighest-coordinateとdegree一base formの積なのでQ+F_b。
他のa=1は113/023でhighest-coordinateだけ、122では零となる。
2Fは二つのfiber又は二倍のfiberを両方保持する。

一般平面切断はrational quintic。Sがgeneric cuspなら五ordinary cuspを持つが、
同既済note §3のdegree八tangent Wronskians→dual degree≤3→bidual degree≤4に反する。
従ってSはgeneric nodeで、§1によりC_S=Γ1+Γ2、両者のL²-degreeは五である。

## 3. C_P=2Fの場合は空

Γ_i=a_i L+b_i Fとする。prime vertical divisorはfiber Fでdegree一なので、degree五のΓ_iはa_i≥1。
C_S=3L-5Fから(a1,a2)=(1,2)、degree比較から
`Γ1∼L`、`Γ2∼2L-5F`（順序を除く）。

122ではH0(Sym²E(-5))=0。
113では全sectionはO(3)座標の平方とdegree一base formの積。
014では全sectionはO(4)座標を因子に持ち、023でもO(3)座標を因子に持つ。
従って残る三型でも全2L-5F divisorは積に分解し、primeではない。
したがってC_P=2Fの全場合は不可能。
P上のrank又はnode/cusp、二fiber/二倍fiberを区別して捨てる必要はない。

014のC_P=2Qも不可能である。この場合C_S=L+5Fだが、分裂した二primeのL係数は
どちらも一以上でなければならず、その和を一には出来ない。

## 4. 残る全closed residue源

113/023のC_P=Q、014のC_P=Q+F_bだけが残る。
いずれもC_PのclassはL-3Fなので、C_S=2L、両degree五primeはΓ1,Γ2∈|L|。
全導手はreducedで、Pの一般型もnodeである。
全codimension一のbranch商とsmooth YでのHartogsにより、αはその全strict conductorに沿う
closed simple-pole一形式となる。small exceptional curveへ新たな極を許していない。

113/023では三residueをλ1,λ2,λQとする。

一般のample complete-intersection curve上のresidue theoremをL+nF（n≫0）の二乗に適用すると
`5(λ1+λ2)+2λQ+2n(λ1+λ2+λQ)=0` が全nで成り立つ。
従ってλQ=0、λ1+λ2=0。
同じresidueを持つ二形式の差は全divisorでpole無し、Hartogsによりholomorphic。
Yはrational scrollなのでH0(Ω1_Y)=0。したがってαの全vector spaceは次元一以下。
固定projectionごとのπ_X、及びH0(Λ²T_P4(-5))=0によるambient πのprojective choiceも高々一点。

014ではΓ1,Γ2のresidueはbranch交換からλ,-λ、Q,F_b上も同じ局所符号計算からμ,-μ。
Γ1,Γ2は同じL classなのでresidue curveで寄与が相殺する。
Q=L-4FはL²-degree一、F_bもdegree一だが、(L+nF)²との交点はそれぞれ1+2nと1。
従って全n≫0で2nμ=0、μ=0となる。
残るresidueはλだけで、H0(Ω1_Y)=0からここでも全α源のprojective dimensionは零以下。

## 5. 全projection parameterの必要Schubert条件

H0(L)は八次元、全三型でK=H0(L-C_P)=H0(3F)は四次元。
plane Pの二linear equationは全C_P上で零となるV内の独立sectionなので、
projection五空間V⊂H0(L)はdim(V∩K)≥2を満たす。

intersectionが二のstratumでは、V∩KのchoiceはGr(2,4)の四次元、
H0(L)/K内の三空間のchoiceはGr(3,4)の三次元、liftはHom(C3,C2)の六次元。
合計十三。intersection三・四のstrataはそれぞれ九・三次元で、同じ上界に含まれる。
finite/birational/conductor/Poisson descent条件はこのsupersourceを狭めるだけである。

113と023ではQは一意なので、C_Pの追加choiceは零。
014ではF_bのchoice一を戻し、全projection sourceは十三+一以下。
actual polarized automorphism subgroupはh0(EndE)-1+3から、113で13、023で14、014で16次元。
Y上の実作用はcomplete L-image Xへfaithfully降りる。base SL2のfinite central kernelは次元を変えない。
Qはそのunique divisor classの唯一の有効元で、014のF_b choiceも作用で保たれるため、
全parameter sourceの条件は実群作用に不変である。
frame付きfinite birational projectionを固定するautomorphismはdense openでidentityなので全体でidentity。
よって全像はそれぞれ

`113: 24+13-13=24`、`023: 24+13-14=23`、`014: 24+(13+1)-16=22`

以下となる。point-vertex coneをrank零coneの既済命題で排除したのではなく、
全relative conductor sourceとSchubert条件から直接評価した。

## 6. 正確な射程

この証明はdegree五rank二・genus零の全normalization五型を処理する。
genus一のcanonical/elliptic分岐、plane及び他degreeは別義務である。
全G/E又は新成分を主張しない。
