# 全 R19\N27 は滑らかで、その全ての点で交叉は通常二重構造である

2026年11月3日。これはR\Nの全射影軌道分類を用いて、一般点の結果を補集合全体へ拡張する定理である。
R∩Nの軌道・局所環・primary decompositionを尽くす主張ではない。

## 定理

1. R\Nは滑らかな有理19次元多様体である。
   12次元のbig paired-flag orbit上のA7-bundleに同型である。
2. CA,CBのproper旗模型は、R\Nを含むtargetの開近傍上で同型となる。
   従って両成分はR\Nの全点で滑らかな28次元多様体である。
3. 全ての閉点p∈R\Nで、実交叉D=CA∩CBの完成局所環は正確に

       C[[u1,…,u19,epsilon]]/(epsilon²)

   である。この部分のDはCohen–Macaulayであり、embedded associated primesはない。

全Poisson完成局所環は、同じ点でも既存LPB39とCA/CB28の存在によって非CMである。
二成分交叉と全Poisson環を同一視しない。

## 1. 残る六軌道の全実計算

AF≠0の一変数族Pκの全点は既存のconstant minorと二つのcokernel covectorsで証明済みである。
全R\N軌道分類により、残るのは次の六代表だけである。

    Q=h²+hz, hz+z², h², z², hz, 0.

それぞれで、同じ実際のcanonical150係数空間において

    rank J_CA=29,
    rank J_CB=29,
    rank [J_CA J_CB]=37,
    rank J_R=20,
    rank [J_CA J_CB Q2diff]=38

をQQ上で正確に検査した。
前二者は滑らかなproper sourceの全affine cone微分であり、R sourceもbig caseの全八係数と群作用からの微分である。
各29次・20次非零minor、全cokernel covector、二次差への非零値を保存した。

二次差は全A,B,Fを記号のまま使っても

    PA(epsilon)−PB(epsilon)=−18epsilon² z²∂x∧∂y

となる。先行の有限source curveを全係数で再使用し、六点ごとの標本だけから式を推定していない。
各代表で二次差を零化しないcokernel covectorが存在するので、全六軌道で真の二次obstructionが非零である。

## 2. 六代表でも全幾何旗fiberは一つ

Q≠0なら、全invariant hyperplanesはh,x,zだけである。
CAのtrace条件 `f_k=5m` とCBの `f_k=(5/2)m`、共通の `f=m*k−C(k)` により、
CA旗は(L,W,H)=(C∂y,span(h,x),Ch)、CB旗は(C∂y,span(x,z),Cx)に一意に定まる。
この計算はAF=0を除外していない。

Q=0ではwもinvariant hyperplaneとなる。
この唯一の追加候補ell=wはm=2x/5を持ち、k=xならf=−2x²、f_k=−4xとなる。
CAの必要値5m=2xにも、CBの必要値(5/2)m=xにも一致しない。
ほかのh,x,z候補は先行と同じ計算であり、新たな旗は増えない。
Q=0でもintrinsic factorは一意、定covector Casimirはないことを正確に確認済みである。

従って全六代表でも全geometric fiberは一つである。
rank29によりprojective sourceの28方向にkernelがなく、fiberは長さ1となる。
proper有限化とNakayamaの論証で、両proper sourceは近傍のCA,CBと同型となる。
群同変性と全軌道分類により、これは全R\Nを覆う。

## 3. big incidenceは補集合全体でproperな滑らかなsourceとなる

全properpaired旗incidenceをR\Nへ制限する。
66表のほかの65 relative casesは全てNへ入るので、その全supportはbig orbit上にある。
big orbitが閉となる旗積の開Uを選ぶと、incidenceの全supportはU内にあり、
big orbitへのbase changeは閉部分スキームとなる。
これが一定rank8 linear intersectionからなる滑らかなP7-bundleの開部分であり、R\Nへの射はproperである。
全incidenceがschemeとしてbig orbitを因子に持つとsupportだけから推論する必要はない。
閉部分スキームとなる滑らかなsourceを使えば十分である。

sourceでt=0ならf_CA=0なのでNへ入り、逆に全t≠0標準形はNに入らない。
したがってこのsourceは厳密にt≠0の部分、すなわちP7内の一つのhyperplaneを除いたA7である。
共通安定化群はt=0 hyperplaneを保つため、12次元paired-flag orbit上のA7-bundleとなる。

全点の幾何fiberは一つであり、rank20はprojective19次元sourceの全接方向がinjectiveであることを与える。
再びfiber長さ1・properness・Nakayamaにより、sourceは全R\Nと同型になる。
従ってR\Nは全て滑らかである。
底はBruhatの議論で有理であり、affine bundleはZariski局所自明なので全体も有理となる。

## 4. 全点での通常二重構造

全代表でCA,CBは滑らかで、実交叉tangent dimensionは29+29−37−1=20である。
CAの28方向内でCB方程式の独立な八つの一次項をimplicitに除き、正則20次元環Sを得る。
被約台Rは滑らかな19次元なので、S内でRのidealは一つのregular parameter(t)である。

既存のgeneric double定理により、残るintersection ideal Iの(t)でのlocalizationは(t²)である。
SはUFDなので、I⊂(t²)となる。
各代表での非零二次obstructionは、Iにt²*unitがあることを意味する。
したがってI=(t²)で、他のembedded relationは残らない。
AF≠0族と六軌道の全てでこれが成立し、全R\Nを尽くす。

この点ごとの完成局所環表示を、全域のnil line bundleが自明、あるいは全域で直和分裂するという主張へは広げない。
局所の二重構造とCM性・embedded primeの不存在がここでの結論である。

## 証明書と未完了部分

* `code/ambient_R19_entire_complement_double_20261103_v63.py`
* `reports/AMBIENT_R19_ENTIRE_COMPLEMENT_DOUBLE_20261103_v63/report.json`

全六代表の29/29/37/20/38 ranks、source minor、CB十九係数、全150-covectorと非零二次値を保存した。
R∩Nの全幾何・全交叉のglobal primary decomposition・全Poisson分類は引き続き別問題である。
