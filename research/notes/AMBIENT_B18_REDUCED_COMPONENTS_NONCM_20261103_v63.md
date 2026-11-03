# B18全体で被約CA/CB/Nは非CM・非正規・特異である

2026-11-03。これは**個々の被約target variety**についての結果である。
source ramification scheme自身が非CMであるという前稿の結論とは区別する。
実CA/CB交叉Dの完成型を決定するものでもない。

## 1. generic B18点での有限moduleの議論

固定代表P0で、受理済みの実imageは三族それぞれについて

    A=R⊕I*t ⊂ B=R[t]/(t²−nu)

と書ける。Rは実target座標により作られたregular complete local ringであり、
CA/CBではcone次元29、Nでは**別に構成した**cone次元28である。
Aは個々の被約targetの完備局所環で、R上有限、dim A=dim Rである。
Iは非零idealでR-rank1をもち、先の独立監査によりmu(I)はCA/CBで10以上、Nで8以上である。

仮にAがCohen--Macaulayであるとする。
R上有限かつlocalなのでm_R*Aはm_A-primaryであり、Rのregular parameter列はAのparameter列でもある。
従ってCM性からこの列はA-regularとなり、depth_R A=dim Rである。
regular local ring上のAuslander--Buchsbaum公式によりAはR上自由となる。

しかしA=R⊕I*tはR-moduleの直和である。
よってIは自由R-moduleのdirect summandで、local R上projective、従ってfreeである。
R-rank1なのでIはprincipalとなるが、mu(I)≥10/8に矛盾する。
したがってAはCMでない。

この証明はnuがR/I上の非零因子であることを仮定しない。
Iの像modulo nuとI/nu Iの同一視、normalizer同士のtensorによる交叉の代用も行わない。
局所環のCM性はcompletionと同値であり、非零coneのscalar smooth因子を除いても同値なので、
projective被約CA/CB/N自身もP0で非CMとなる。

## 2. 非正規性と全B18への閉包

P0の全proper normalizer fiberは一つの台を持つdual number長2である。
P0近傍でnormalizerは有限birationalであり、targetがnormalなら有限normalizationは同型になるはずである。
その場合のfiberはreduced長1なので矛盾し、三targetはP0で非正規である。
regular点はnormalであるから特異点でもある。

P0の全150係数GL微分のrankは19であり、scalar方向もその像に入ることを別に確認する。
従ってprojective PGL軌道次元は18である。
既証のB18は18次元既約G不変閉集合でP0を含むので、この軌道はB18に稠密である。
三targetとCM性・normality・regularityはG不変なので、全軌道で同じ失敗が成立する。

有限型複素varietyのCM locus・normal locus・regular locusはいずれも開である。
従って各失敗の閉集合は稠密軌道の閉包B18全体を含む。
結論として

    B18 ⊂ NonCM(CA) ∩ NonCM(CB) ∩ NonCM(N),
    B18 ⊂ Nonnormal(CA) ∩ Nonnormal(CB) ∩ Nonnormal(N),
    B18 ⊂ Sing(CA) ∩ Sing(CB) ∩ Sing(N).

これは**包含**であり、これらのbad lociがB18に等しいという主張ではない。
他の非CM点、例えば既知の非零Qa上の点を排除していない。
各環のexact depthや、actual交叉Dのnilradical・embedded primes・primary型は別の未解決課題である。

* code: `code/ambient_B18_reduced_components_noncm_20261103_v63.py`
* receipt: `reports/AMBIENT_B18_REDUCED_COMPONENTS_NONCM_20261103_v63.json`

独立受理済みfinite-image/odd-ideal入力をハッシュ固定し、GL軌道の全係数rankとscalar包含を追加検証する。
