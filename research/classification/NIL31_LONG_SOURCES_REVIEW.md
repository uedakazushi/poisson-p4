# 冪零(3,1)長い鎖top：二つの全五源の独立確認

2026-11-07。主課題 `E_multiple_hyperplane`。
対象は `MODULAR_NIL31_LONG_SUPPORT.md` のE源(r=2/3)とF源(r=3/2)である。
未決だったr=−1の二源は今回再検査していない。

## 入力と有限検証

入口は `CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json` のcase1、
全11変数の六括弧と六残Jacobi条件である。
全90係数からこの線形源への還元は既存依存とし、今回のreviewでは再計算しない。
全源の代入は `code/data/classification_nil31_long_sources.json` に従い、
producerの保存式とsource/spec/hashを照合した。

独立checker `code/classification_nil31_long_sources_review.py` は、
標準ライブラリのFractionによる疎多項式演算を用いる。
producer、SymPy、NumPy、legacyのJacobi・階数関数はimportしない。
既存の独立reviewで用いた算術基本関数を再使用している。

有限job `jobs/nil31_long_sources_review.json` をgate後に実行し、1.44秒で全項目PASS。
結果、standalone checker全bytes、依存hash、非零minorは
`reports/CLASSIFICATION_NIL31_LONG_SOURCES_REVIEW_20261107.json` に保存した。

## 全パラメータの恒等式

Eは自由変数(n2,n3,n7,n8,n9)、rank四開n8≠0。
Fは自由変数(n2,n3,n4,n7,n9)、rank四開n4≠0。
代入後の全六残条件と全10 strict Jacobiが恒等的に零になることを確認した。
非斉次Pfaffianは、(x,y,z,w)=(x1,x2,x3,x4)として厳密に

    E: −n8 [x−yz+z³/3+n2+(3n3/5)(w−n9)],
    F:  n4 [y−z²/2+n7].

それぞれx,yの係数から、記載したrank四開は必要十分であり、いずれも既約である。
これらは五自由係数を標本に固定した恒等式ではない。

分母を除く倍率はEで15、Fで2。
全源をhで斉次化し、全係数の保存式・Pfaffianを再現した。
非斉次modular場は両者とも

    N=x2∂x1+x3∂x2+∂x3

であり、斉次化すると `N_h=x2∂x1+x3∂x2+h∂x3` となる。
倍率sを含むPに対し、div P=sN_h、
`Π=P−(s/5)E∧N_h` のdivergence零と全射影Jacobiも確認した。
affine源自体をdivergence零と仮定した計算ではない。

## 全接空間と成分性

各自由係数を全て1にした同じ有限代表を用いた。
全150二次bivector係数、全350射影Jacobi微分D、
Euler25列・GL5方向25列・scalar一列・実際の五源微分を全て保持した。
全56列Fに対してDF=0をFractionで厳密確認した。

| 源 | rank_QQ D | Euler実効階数 | 軌道・Euler・scalar | 全実源込み |
|---|---:|---:|---:|---:|
| E(r=2/3) | 103 | 24 | 47 | 47 |
| F(r=3/2) | 103 | 24 | 47 | 47 |

modular非零minorはQQ上の下界にだけ用いた。
103+47=150および厳密DF=0が上界も与えるので、表示値は正確なQQ階数である。
Eulerの上界24は、その25列の対角和が零になる恒等式も保持している。

従って射影全接空間上界と実族下界はともに22である。
PGLと各既約五源の像閉包は22次元の全成分となり、代表点は滑らかである。
軌道だけでも22次元なので、各成分には稠密開軌道がある。
これは特定の五パラメータ値だけの存在証明に留まらず、
その点を含む既約全源の像閉包全体の成分性を与える。

既存成分との因子型・canonical cofactorによる相違、
他のnilpotent源を含む全支持の被覆は別の幾何証明に依存する。
本reviewはこれらを再証明したとも、全Poisson空間を分類したとも主張しない。
