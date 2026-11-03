# actual CA/CB odd generatorsのU³ class

2026-11-03。前稿の二次欠損証明を変更せず、全critical monic消去に沿うU-axisでactual三次classを計算する。
full Iの生成数上界やCA/CBの全完成型を決定した結果ではない。

## 全monic sectionを使う

CAではs3=3U/10,s7=7U/20、CBではs3=s7=U/2とし、W=V=0とする。
七つのmonic critical消去式を**二次・三次項も含めて全て**代入する。
CAではe5=6U²/25、CBではe5=5U²/27となる。
source fの第三係数やHの三次補正を落としていない。
両族ともratio gはこのsectionで0だが、CBのratio微分列−f*Gも可変のまま保つ。

全source JacobianをJ(U)とし、前稿と同じ共通actual28target Yを保つvertical vectorをNeumann三次まで求める。
各次数dについて

    c_d=−B0_Y^(-1)*(b_d,Y+Σ_(i=1..d) B_i,Y*c_(d−i))

である。全150 output derivativeの次数0,1,2が厳密に零、さらにbranch row61の三次も零となる。
従ってtracefree tはこのsectionで次数4以上となり、even成分は次数3に寄与しない。
全150三次係数はactual odd b_iのU³制限そのものである。

## 以前のclassを全てsubtractする

まずactual R全firstjet行列の七basis rowの定数結合をsubtractし、zero-firstjet actual generatorsを作る。
次いで二次WU,WV classの二basis row46,95に対する全係数をsubtractする。
その残るactual三次classは、raw row63のみに非零となり、

    CA: (12/125)*U³,
    CB: (−4/243)*U³

である。

このsectionでは全J=I+(nu)のgeneratorは次数3以上に制限されるので、
m_R Jの元の制限は次数4以上である。
従って上の非零U³係数はJ/m_R Jの最後のcubic minimal classを検出する。
actual Iの像は七つのlinear class、二つのquadratic class、一本のcubic classの計10次元を持つ。
既知J/mJは11次元で、欠けるものは前稿で証明したUV quadratic classである。

I自体のminimal生成数は少なくとも10となるが、I∩mJをmIと同一視してはならない。
したがってactual Iが十生成、あるいはIの全idealが表示した十元で生成されるとは結論しない。
nu∉Iの結論は前稿のまま有効である。

候補I=(s6,nu+cUV,WU,WV,U³)はleading dataと整合するが、
高次で追加generatorが入る可能性、nuのnilpotent order、exact depthは未解決である。
R19\Nは個々のCA/CBで滑らかなため、個々のconductorの追加枝として使わない。

* code: `code/ambient_B18_CA_CB_odd_cubic_class_20261103_v63.py`
* evidence: `reports/AMBIENT_B18_CA_CB_ODD_CUBIC_CLASS_20261103_v63/`

全source section、variable Jacobian、Neumann四項、全三次classと以前のclassのsubtract係数を保存する。
