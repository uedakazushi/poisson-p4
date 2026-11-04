# 固定G115／target19のD6：独立QQ source・積戻し監査

2026-11-04。対象は一つの固定六minorのみ。producerをimport／再実行せず、別monomial順
DEGREVLEXの通常QQ多項式積を用いる。旧120-permutationの大きなparameter式監査を反復しない。

全1325 source演算のatom payload、既知QQ/a分母、Laurent offset、固定source列を照合する。
行名G20/G115/G85はそれぞれ**current15の2/11/7**というprefix後の行の略称であり、裸の
actual135 generator番号20/115/85を意味しない。full57への正確なsource持上げは既受理の
演算子が全prefix補正を保持する。四行P4についても同じ略称の規約を使う。

bottom行は保存actual行から文字通り
Ksmall=g*current11/c10−h*current7/c6を構成し、literal current2をその上へ置く。
選択列は元32の[22,11,13,25,15,19]。top第二・第三行だけからaを一つずつ除くため、
raw D6=a²*normalized D6である。Bの未知除算やp/R/h22の反転を追加しない。

十五のtop4-column選択、各top4det、bottom2det、Laplace符号を自由matrix-entry多項式で
独立確認する。最初のP4detだけは既受理の厳密な保存P4係数証明へ戻し、二つのa行factorも
含める。十五の共通QQ分母scaleが整数であることを個別確認する。nativeの整数格子加算は
使わず、全sourceを最後の使用で解放しながら十五のQQ積をin-place加算して、保存9047係数と
全て比較する。保存結果はsource値を全解放した後だけ読み、live上界147444を守る。

結論はこの指定D6の全QQ恒等式と非零性である。D5は展開しない。新しいtyped companionが
元4×11の一成分へ結ぶまでは、それ以外のpivot、rank4、global57、quartic239を結論しない。
別の行・列・点・素数・GB・rank searchを計算していない。

再現: `python code/s61_nonsquare_prePOWER_G115_t19_literal_pair_entry_review_20261104_v63.py`。
