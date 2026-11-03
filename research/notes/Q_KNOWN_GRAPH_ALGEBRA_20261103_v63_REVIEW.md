# 正確なcolon・全模型非平坦性・二既知素idealの剥離：独立監査

2026年11月3日。T39/T40両方のQQ九変数商と全節点差商を独立に検査した後、その純代数的帰結を監査した。
二つのfrozen proof-receipt producerも再実行して完全一致した。新たなGBや巨大行列式展開は不要である。

P/Nが整域、Cbar≠0、CN⊂I⊂Nなら、fC^k∈IをP/Nへ移すことでf∈N。
逆包含もCN⊂Iから従うので全k≥1でI:C^k=Nとなる。
従って×CによるP/N→P/Iはwell-definedな単射であり、正確な加群列が得られる。
これは単位的な環の埋込みや、R全体の直積を意味しない。

N39とN40はそれぞれa−2/9とa+8/9を含み、対応するCを全R内で零化する。
Cのgraph像は非零多項式で、各Qaのgraph原点の完成でも非零である。
そのことは正確なc39*epsilon^44とc40*lambda^275の非零curve restrictionからも従う。
従って **全Rの局所環・完成環** はa方向でtorsionを持ち非平坦となる。
閉じたsubsliceの非平坦性を全環へ移した議論ではない。

二つの一次基底因子の差は10/9で単元だから、C39*C40=0、(C39)∩(C40)=0。
従って和idealは加群としてP/N39⊕P/N40、C39+C40のannihilatorは(N39∩N40)/Iとなる。
これを反転したときの二つのidempotent C_i/(C39+C40)が、二つの既知graph局所化への正確なscheme直積を与える。
全Rの直積分解ではない。

さらにf∈N∩(I+C)をf=i+Cgと書けば、非零Cbarとdomain性からg∈N、従ってCg∈Iとなる。
よって正確にI=N∩(I+C)。またC40∈N39,C39∈N40であり、順に同じ論証を使うと

    I=N39∩N40∩K,   K=I+(C39,C40)

となる。Kは各Ciを含みCi∉Niなので、既知素成分全体を再び含むものではないが、その部分境界を残し得る。
残余Kの被約成分・primary structureの分類は依然未解決。
Qstarではa+8/9が単元なのでC40∈I,N40=Pとなり、局所化・平坦な完備化後も
Ihat=N39hat∩(Ihat+C39)が成立する。

再現コードは code/q_known_graph_algebra_review_20261103_v63.py、同名 REVIEW report が両側の係数監査をhashで固定する。
