# 一次文献によるSp4不変関数体の純有理性：独立監査

2026-11-04。複素数体上だけの結論である。一次PDFをmanaged web toolで開き、
印刷頁18・21（PDFの0-based頁6・9）の画像も直接確認した。

N. I. Shepherd-Barron, *Invariant theory for S5 and the rationality of M6*,
Compositio Mathematica **70** (1989), 13–25。
[一次PDF](https://www.numdam.org/article/CM_1989__70_1_13_0.pdf)

**文献の適用範囲。** Proposition 11（印刷頁18、証明は21頁まで）は、滑らかな
quadric Q 上の `P(V(1,2))×P(V(2,1))` を全 Aut(Q) で割った商の有理性を
証明する。連結成分だけでなく、二つのfactorを交換する成分を含む。
Theorem 15（21頁）は五次元既約S5表現の射影商についても有理性を与える。
character twistに関するLemma 8（16頁）とも整合する。ここでは前者の直接の
二形式同定を主たる適用に使い、M6の有理性だけから別の商へ類推しない。

**今回の帰結。** 独立受理済みの五点・vector descent packetでは、まさにこの
projective-pair quotientの体をBとし、元W16のSp4不変関数体について
`F=B(s1,s2)` を証明した。Bは次元4で、二倍率を落としていない。
上の既存定理によりBはC上純有理、従って

    C(W16)^Sp4 ≅ C(z1,z2,z3,z4,z5,z6)

という **純有理性・超越次数6の存在定理**が従う。
これは1989年の定理への適用であり、新しい有理性定理として扱わない。
従来のstable level≤2は正しいままだが、この帰結によりpure rationalityが得られる。

本レビューは一次定理の読解と仮定の照合であり、原論文の証明全体をCASで
再証明したものではない。六つの明示生成元を今回構成したとは言わない。
QQや任意標数への移送、NEW25の方程式・次数上界、NEW sourceとtargetの有限指数、
全多項式不変式環・全orbit分類も結論しない。新しい点・係数・GB計算はない。

直接引用は用いず、一次文献の記述は該当命題の短い要約だけに留めた。
