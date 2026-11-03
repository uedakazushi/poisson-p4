# 固定12一次切断の有限証拠と失敗scope：独立監査

2026年11月3日。既存GBは再実行せず、known-plane、低次数native projection/rank、shadow assessmentの
軽い三つだけを再実行して照合した。nativeのSym2/Sym3 product cache、monomial multiplicities、
全74入力とFF rank実装もsourceで監査した。計算は標数32003の正確な有限行列でありQQ rank等号へ移していない。

known13平面に制限した12×13整数行列とprimitive kernelを別に検査し、p単元minor46494と
f(kernel)modp=11541を確認した。既存の全cubic restrictionと合わせ、この **既知平面内の** 切断は
V[t]/t³で正確に長さ3である。全23変数sectionがArtinianという意味ではない。

GFでrankQ57、linearQ1210、全degree3が1227なので、sectionのHF≤3は(1,23,219,1073)。
linear syzygy空間は101次元である。module項順序で101個の異なるleading linear module monomialsを選び、
各componentのr_i変数を掛けて得るdegree2 leading shadowはg(r)=23r−r(r−1)/2以上。
concavityにより最小配置は23,23,23,23,9、合計1275である。
これらはdegree4 Macaulay kernelへ独立に入り、h4≥14950−16123+1275=102となる。
有限なら長さ≥1418。この範囲で全degree4空間を張るMacaulay証明は不可能だが、
低次数GBによるArtinian証明自体は排除しない。新しいquartic initial generatorsも未計算である。

単独GBのsignal9は最初のcheckpoint前で、cofactor証拠がなく、停止原因も特定されていない。
観測RSSだけから3GiB不足と断定せず、未完GBから正次元/零次元/上界12を結論しない。
上のGF h4下界をQQや全local initial idealへ移さない。既存の全R35上界22は変更しない。
code/q_linear_section_review_20261103_v63.py と同名 REVIEW report が有限証拠とscopeを保存する。
