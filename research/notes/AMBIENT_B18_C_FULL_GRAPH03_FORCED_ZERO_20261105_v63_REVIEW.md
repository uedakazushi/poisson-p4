# 全 C の三次次数五又は六と五 graph 係数の零：独立監査

独立報告 `b86b737f...` は主報告 `2f6ee029...` の有限 D 加群による障害を確認した。
同じ full C、D=QQ[[ell]]、m と、その商 C3=C/ell³C を区別する。
既知の厳密な m³C=ell m²C と C3 の実 HF(3,6,6,4,0) を使う。

F=m³C とおく。全 C の HF3 が四なら、mF=ellF、rank_D F=4 より F は D-free 四。
ellF=m⁴C⊂ell³C⊂F で、F/ellF と F/ell³C はともに長さ四である。
従って ell³C=ellF=m⁴C となる。
任意の x に ell³x=ell y、y∈F を選ぶと t=ell²x-y は ell-torsion である。
u∈m に対し ut∈F も torsion なので零となり、m ell²C⊂ell³C を得る。

一方、完全155 source・alpha 四値零・CUT1 により確定した実 C3 の三次基底は
A r²ell、A r ell²、B r²ell、B r ell² である。
独立 checker は no-cut NF から cut NF を組み直し、r ell²A と r ell²B の列が
文字通り異なる二つの単位ベクトルであることを検査した。
従って上の包含は矛盾し、全 C の HF3 は五又は六である。

いずれかの graph d03 が非零なら、既に独立確認した二列 minor d03² により
全 HF3 は四以下となる。従って五行63、64、68、73、74の d03 は全て source だけで零。
graph03 数値 worker の出力 directory が存在しないことも確認した。
準備済み SOURCE・code・custody は未実行のまま保持し、不要な数値監査は行わない。

長い指数は依然 e>=3 であり、e=3 を排除していない。
Tor_D(m³C) と m³Tor_D(C) の等式は仮定しない。
これは comparison C の結果で、元 H0 又は元 nil ideal の三乗を決めるものではない。
再現は `code/ambient_B18_C_full_graph03_forced_zero_review_20261105_v63.py`。
