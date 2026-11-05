# 元H0のell消滅を判定する三つのBockstein証人（条件付き計算法）

2026年11月6日（日本時間）。これは既証DVR構造からの一般的帰結と別ルートの候補であり、実際のBockstein像を計算した結果ではない。

D=k[[ell]]とし、**元の**比較対象を表すD-flatな複体Fを使う。H1(F)=KのDVR型はD自由四個とD/(ell^ei)三個である。Fbar=F/ell F上のBocksteinは、cycle xbarを持ち上げたxについてdx=ell*yと書いて、beta[xbar]=[y mod ell]と定義する。これはbeta:H2(Fbar)→H1(Fbar)=Toldを与える。

同じ複体のuniversal coefficient完全列は

    0 → H2(F)/ell H2(F) → H2(Fbar) → K[ell] → 0,
    0 → K/ell K → Told → C[ell] → 0

である。定義からbetaは最初の右端K[ell]を経由し、その像をK/ell Kへ自然に入れてからToldへ入れる。D/(ell^e)のell-kernelの生成元ell^(e−1)は、K/ell Kにおいてe=1のときだけ非零となる。従って

    rank_k beta = #{i : ei=1} = 3-j,
    length T2 = 17+j = 20-rank_k beta.

特に、正しいFbarのcycleを三つ構成し、それらのBockstein像が**実Told**で独立と確認できれば、H2(Fbar)全体の基底や全617 upper像を計算しなくてもrank beta=3となる。既証の三torsion summandが上界三を与えるからである。そのとき全ei=1、元H0の長さ三とell H0=0を得る。既証のN³⊂ell H0から元nil三乗も零となる。二証人だけならT2長≤18、一証人だけなら≤19という片側上界に限る。

このルートには重要な前提がある。元二項対象を単に途中で打ち切った自由表示のH2を、K[ell]と同一視してはならない。実際の全複体を保つか、同じ導来対象を表すことを証明したD-flat置換が必要である。正しいFでもH2(F)/ell H2(F)が残り得るが、その部分のbetaは零なので、上の像の議論には影響しない。有限conductor商で先にkernelを取ったものも代用できない。

現行の完全C2比較・有限factored homotopyはsourceとして存在するが、三つの適切なcycleとRegDivell後の像が既存の低次cacheだけで得られるかは未確認である。高次Phiやhomotopy係数の新照会、任意のtruncated kernel計算をこのノートから許可したことにはならない。現行case554と接続写像の低次計算を先に進め、失敗した場合のSOURCE候補として保持する。

依存する既証記録：`AMBIENT_B18_ELL2_HONEST_T2_DVR_TYPE_20261105_v63`、`AMBIENT_B18_ELL2_HONEST_T2_LENGTH_BOUNDS_20261105_v63`。元比較の打切り注意は`AMBIENT_B18_COMPLETE_NORMALIZER_COMPARISON_PREFLIGHT_20261104_v63`にある。新しいrank、T2長、H0長、nil指数はこのノートでは確定していない。
