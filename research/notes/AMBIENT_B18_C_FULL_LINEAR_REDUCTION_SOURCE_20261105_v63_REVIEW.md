# 全比較加群 C の ell 還元：独立 SOURCE 監査

同じ実 regular 局所環 R0、極大 ideal m、D=QQ[[ell]] 上の比較加群 C を扱う。
主証拠 `27575d10...` を author 実装を import せず監査し、独立報告
`fe2efb42...` を固定した。新しい graph 係数・physical jet・gauge は読んでいない。

まず有限 D 加群であることを source から確認した。二つの実 branch は、
既知 unit graph change 後の r 三次係数が QQ 単元なので、それぞれ D-free 三である。
active beta 二次の固定 minor は -5/192 であり、IA+IB は m-primary となる。
T=(IA∩IB)/(IA IB) は IA+IB で消える有限長加群である。
基本完全列と二 branch への埋込みから R0/(IA IB)、従って C が有限 D 加群と分かる。
これは K/JK の Tor 補正を省略する議論ではない。

ell 三乗で切らない三次 target 全495単項を独立に列挙した。
z を含む483項、H の平面四項、A/B の r 三乗二項を既知 initial で消すと、
A/B の r²ell、r ell²、ell³ が各三項残る。ell による二次からの写像は I6 である。
従って gr3 C=ell gr2 C、m³C=ell m²C+m⁴C となる。
Q=m³C/(ell m²C) に Nakayama を適用し、厳密に m³C=ell m²C を得る。
C 上で ell を消去してはいない。

五 graph source について、実 mixed delta 二次十四係数を各行で確認し、
tauA/tauB の平面三次がそれぞれ d03 ell³A、d03 ell³B になることを監査した。
別の普遍 minor companion は、六次元 NF の自由 ell³ 座標から minor d03² を再構成する。
非零 d03 があれば全 HF3=4、m³C が D-free 四となる、という条件付き含意自体は正しい。
その前提は後続の `C_FULL_GRAPH03_FORCED_ZERO..._REVIEW` が純粋に排除した。
条件付き torsion 型1123や full HF366444…を実際の結論として採用しない。

元 H0、元 nil 三乗、長い torsion 指数 e の値はここでは決定しない。
再現は `code/ambient_B18_C_full_linear_reduction_source_review_20261105_v63.py`
と `code/ambient_B18_C_full_graph03_minor_review_20261105_v63.py`。
