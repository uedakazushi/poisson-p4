# 22最大minor多項式のexact GCDとlocalizer支持因子の除去

GBを再試行せず、既存22多項式のGCDだけをnative FLINTで計算する。全22個の
literal quotientと再構成を保存し、別のFraction積和で確認する。greatest性は
逆順の別実装SymPy exact multivariate GCDで照合する。点・有限体の結論を使わない。

gからgcd(g,U)を繰り返して除き、各除去因子がUをexactに割ることも保存する。
最終hとUのcoprimalityを独立GCDで確認する。a,d,F,B,C,Jの既約性・相異性を
仮定せず、未知因子のfactorizationやprimary decompositionを行わない。

hが非定数なら、全minorに共通の非定数因子がD(U)上にも残り、B-rankbadの
実際のdivisorを与える。hが定数なら、この開のheight1共通因子が無いことだけが
分かり、rankbad supportは空又はcodimension>=2（dimension<=1）。unit idealや
FULLU全域rank3はそこから推論しない。fullquartic239も別の課題である。
