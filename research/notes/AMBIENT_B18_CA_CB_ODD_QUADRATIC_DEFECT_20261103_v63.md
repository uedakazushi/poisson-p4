# actual CA/CB odd2jetの欠けるUV class：nuはIに入らない

2026-11-03。generic B18での各actual image A=R⊕I*t⊂B=R[t]/(t²−nu)を使う。
対象は個々の被約CA/CB target完成環であり、actual交叉Dの全完成型ではない。
新しく決定したNのI=I+(nu)を、CA/CBへそのまま延長できないことを示す。

## 全sourceと二次vertical微分

全sourceの順序はfixedfiber20＋全旗9である。
CAのfixedfiberはlinear、CBではq*fを含む全非linear sourceを使い、
fiber微分もF_i−qG_i、−Σf_iG_iを保つ。
CBではflag列にも二次source項があるので、Neumann計算にJacobianの二次項を含める。

前稿での三変数critical tangentを(W,U,V)とする。
CAではs3=(3U+V)/10,s7=(7U−V)/20,s0=(V−U)/10,s2=W、他のsは0。
CBではs3=(U+V)/2,s7=(U−V)/2,s0=−7U/9+V,s2=W、他0である。
これは全source critical schemeのmonic消去の**一次部分**に沿う実source族である。

source旗E_xh列（全列23）をfreeとし、他28列をB0+B1+B2、free列をb0+b1+b2とする。
共通28 actual target Y行に対して、次数2までその値を保つvector fieldを

    c0=−B0_Y^(-1)b0_Y,
    c1=−B0_Y^(-1)(b1_Y+B1_Y*c0),
    c2=−B0_Y^(-1)(b2_Y+B1_Y*c1+B2_Y*c0)

で求める。
全150 target微分の次数0,1は零、二次はS2=b2+B0c2+B1c1+B2c0である。
Y全28行とbranch row61のS2も零である。

Nと同じrank2 traceの議論で、source tracefree tのcritical tangentへの制限は次数3以上となる。
Dt(0)=1、全odd b_iの一次source制限は零なので、S2はactual b_iの二次制限そのものである。
source ramification idealの生成元をactual odd係数と呼び替えたものではない。

## all actual generatorの二次classを計算する

既存のactual R29 odd一次行列のrankは7である。
選んだ七rowの基底に対する全150rowの線形従属係数を取り、
各actual b_iからその七元の定数結合をsubtractする。
こうして**R全firstjetが零**となる全odd generatorを作る。

その二次制限の張る空間は、両族で正確に

    span(WU,WV)

という二次元空間である。
全150rowを検査し、UV,U²,V²,W²の係数は全て零である。
独立なrow46,95の小行列式も保存する。

一方、actual generatorの定数結合でfirstjetが正確にdnuとなるchosen lift i_nuを作ると、

    CA: i_nu,2 = −UV/75,
    CB: i_nu,2 = −16UV/9

となる。これはchosen liftのclassであり、一次0の別元を加える自由はある。
しかしその自由の二次制限はWU,WVだけなので、UV成分を消すことはできない。

## nu∉Iの直接証明

仮にnu∈Iなら、i_nu−nuはI内でR-firstjetが零の元となる。
ところがsourceのtはcritical tangentで次数3以上なので、nu=t²の制限は次数6以上である。
従ってi_nu−nuの二次制限は上の非零UV成分を持つ。

任意のfirstjet零のIの元をΣf_i(R)b_iと書く。
その二次制限でf_iの非定数部分は寄与しない。なぜならb_i自身の一次source制限が零だからである。
定数部分はactual一次行列の関係であり、既に全150rowの基底引算で全て計算済みである。
従ってその二次制限はWU,WVのspanに入らなければならず、UVは入らない。
これは矛盾である。

結論は**nu∉I_CA、nu∉I_CB**であり、tracefree t³はそれぞれのactual imageへ入らない。
またJ=I+(nu)に対するNと同じfull Nakayama surjectionは成立しない。
これは以前の撤回済みR19-conductor推測とは独立な、actual係数の厳密な有限証明である。

なおR\N上では個々のCA/CBは滑らかでnormalizer旗fiberも一つである。
CAの旗とCBの旗が異なることを、個々のnormalizationでの二preimageと取り違えてはならない。
R\Nを個々のconductorの第三枝として使う議論は本稿では一切使わない。

## 未決の高次部分

全Iはまだ決定していない。
leading dataはI=(s6,nu+cUV,WU,WV,U³)のような候補と整合するが、
actual cubic class、全生成数の上界、formal同値の全次数性は未検証である。
nuがR/Iでnilpotentか非零因子か、exact depth、normalizationのsubintegralityも今回結論しない。
actualDの全primary型も別の問題である。

* code: `code/ambient_B18_CA_CB_odd_quadratic_defect_20261103_v63.py`
* evidence: `reports/AMBIENT_B18_CA_CB_ODD_QUADRATIC_DEFECT_20261103_v63/`

全physical Jacobian、variable fiber列、Neumann012項、全150 actual二次制限、
全zero-firstjet生成元、UVを含むchosen nu liftを固定する。
