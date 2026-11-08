# Co-Higgs Borel・高multiplicity全tau非零源は不変hyperplaneを持つ

2026-11-08、PRIMARY E、log_manuscript。producer純証明／保存係数のreadback。
新operator、M1、rank、minor、point、GB又はjetは作らない。
1054主稿／frontierとglobal／Git／旧bytesを変更しない。
独立受理は別義務。

## 1. 全sourceと結論

rank零plane、O(1) co-Higgs、rank M二、M(p)非零、高multiplicityの全入口は、
NDSS0、NDN0、BSS0(tau)、BN0(tau)である。
全Lax還元／全actual normalizer逆は
[既全high-m独立review](PLANE_HIGH_M_RANK_AT_LEAST_TWO_SEPARATE_REVIEW.md) §2による。
ND二源とtau零二源は
[定数rank束の全rank二証明](PLANE_COHIGGS_NONKERNEL_ZERO_LAX_WHOLE_RANK_TWO.md)に送る。
本稿では残るBSS0／BN0の**全tau非零、全higher**を扱う。

結論は、全sourceに実Poisson hyperplaneが存在し、actual五次が一次因子を持つこと。
従って被約既約quinticを持つrank四源は全空である。
全D2自由核、非斉次particular及び全normal位数を保持する。

## 2. 元のweight零でw-rowは常にw倍

base(a,b,c)、normal(z,w)、M=bH+cE、En=z∂z+w∂wとする。
全weight零はP0=T∧En+U∧((bz+cw)∂z−bw∂w)+K∂z∧∂wである。

BSS0では

~~~
U=∂b, T=(−2tau a∂a+4tau b∂b−2tau c∂c)/3,
K=(1−2tau)zw.
~~~

従ってPaw=−2tau aw/3、Pbw=(4tau/3−1)bw、
Pcw=−2tau cw/3、Pzw=(1−2tau)zw。

BN0ではU=∂c、T=2tau b∂c、K=(1/2−tau)w²。
従ってPaw=Pbw=0、Pcw=(2tau−1)bw、Pzw=(1/2−tau)w²。
両型とも全tauでP0の全w-rowがw倍である。

## 3. 全P1に残り得るbad項を三つだけ同定する

旧全30方向の順序は

~~~
0..8:   z*(linear base)*base pair,
9..17:  w*(linear base)*base pair,
18..20: z²*(∂a,∂b,∂c)∧∂w,
21..23: zw*(∂a,∂b,∂c)∧∂z,
24..26: zw*(∂a,∂b,∂c)∧∂w,
27..29: w²*(∂a,∂b,∂c)∧∂z.
~~~

従ってw-rowをw倍でなくし得るのは、全P1のrows18,19,20だけである。
P1のbase-base、他のnormal成分、全P2を零とは置かない。
特に**FULL九P2は全てbase pairsにnormal quadraticを掛けた項**であり、
w-rowに一切寄与しない。kerD2又はparticularのどちらでも同じである。

## 4. 保存済みの完全kernelはbad三rowを全て零にする

不変入力はCLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_HIGH_M_20261107.json。
BSS0_tau／BN0_tauのanalysis.D1.complete_kernel、
及びanalysis.all_exceptional_rootsの全非零parameterのD1.complete_kernelを読む。
全350行で独立確認済みの**完全**kernelである。

|全source patch|保持するparameter|全kernel次元|rows18/19/20|
|---|---|---|---|
|BSS0 QQ(tau) generic|元selected-minor open|7|全entry 0|
|BSS0 fixed wall|−3/2,1/6,1/2,3/2,−1/2|7,7,8,7,7|全entry 0|
|BN0 QQ(tau) generic|元selected-minor open|7|全entry 0|
|BN0 fixed wall|1,−3/2,1/2,−1/2,3/2|7,9,9,8,10|全entry 0|

元generic minor／kernel分母の全因子と全根の別作者認証により、
この二generic patchと全十非零壁でtau非零全体を尽くす。
generic-rankやsample値から零境界を外挿しない。
tau零のbad三rowには実際に非零entryがあり、本表へ混ぜない。
その全源は別のrank八束証明でrank二以下へ送る。

同reportのD2 nullityはBN0 tau一で二、tau半で一に増えるが、
その全free P2もbase pairsだけなので以下の因子証明には影響しない。
全非斉次P2を保存するため、D2 inverse又はkernelによる消去すら不要である。

## 5. 元actual五次の因子

任意の全Poisson solutionではP1は上記完全kernelに入る。
従って全P=P0+P1+P2の全w-rowがwで割れ、
H=(w=0)はactual Poisson invariant hyperplaneである。
Hのprojective tangent dimensionは三なので、
その上のbivector squareは零、すなわちactual退化五次はH上恒等零となる。
従ってwは元のsigned Euler-Pfaffianを割る。

直接にも、wを含むcomplementary raw Pfaffianはw-rowによってw倍、
wを含まない唯一のcomplementには外側Euler係数wが掛かる。
よって全自由higherを含む全actual五次がw倍である。
F零ならrank四でなく、F非零なら一次因子を持ち指定既約性に反する。
Euler補正、全actual normalizer逆と全PGL frameを戻してもこの結論は不変。

これはwhole necessary linear source自体に働く因子であり、
その内部の未計算nonlinear compatibilityを捨てた結果ではない。

## 6. 高multiplicity co-Higgs rank二以上の集約

本稿のtau非零二whole族を、
NDSS0／NDN0／BSS0(0)／BN0(0)の全rank二証明と合わせれば、
M(p)非零の高multiplicity全四源は空。
M(p)零の高multiplicity全源は旧独立reviewの共通U／二次因子／三row矛盾で既に空。
rank M三の全高multiplicityは
[tau零全二源定理](PLANE_COHIGGS_FULLRANK_HIGH_M_WHOLE_RANK_TWO.md)でrank二以下。

従ってこれら新producer純証明を独立受理した後には、
O(1) co-Higgs、rank M≥2、m_D≥3の全枝が指定rank四・既約quinticで空となる。
m二のM(p)非零の残finite壁、rank M一以下、O(0)、scalar及びzero-symbolは別義務。
本稿から全rank零plane又は全Eの完成を結論しない。
