# 全 C の四次上界と一変数の高次 HF：独立純粋監査

主 SOURCE `deed64d5...` に対する独立報告 `25856b35...` を固定した。
実 full C の HF 初四項(3,6,6,6)、m³C=ell m²C、D-rank 四を使い、
新しい係数・jet・gauge は計算していない。

全 no-cut 四次 NF を curve の再帰還元で組み直した。
その行零と行三は、旧 C3 の Ar²ell²、Br²ell² の射影に文字通り一致する。
実 HB0=ell P0-r P1 は、ell³ cut より前に全三次が相殺する source であり、
行零の既知値 -3/5 は全 C でも非零の四次 initial となる。
従って六次元 cover から少なくとも一関係が落ち、HF4<=5。
rank と ell 還元の下界から HF4 は四又は五である。
未知の r ell³、ell⁴ 座標を零とはしていない。

別の独立証明は HB 値を必要としない。M=m²C、S=Tor_D M とおく。
M と ell M はともに rank 四・最小六生成なので、S は指数がともに二以上の二巡回成分を持つ。
ell²S は、比較 C の型(1,1,2,e)から得られる巡回加群 ell²Tor_D C に埋め込まれる。
従って二成分のうち一方の指数は正確に二であり、
M≅D⁴⊕D/(ell²)⊕D/(ell^b)、b>=2 となる。ただし D 加群としての非標準的同型だけである。

厳密な ell 還元から、j>=2 で

    h_j(C)=4+1_{j<4}+1_{j<b+2}

となる。Hilbert 級数は (3+3t-t⁴-t^(b+2))/(1-t) である。
特に HF4=4 は b=2、HF4=5 は b>=3 と同値である。

さらに ell²Tor_D C⊂S から b>=e-2。
実 C3 の m⁴=0 により ell²M⊂ell³C である。
有限 D 加群には Tor_D(ell^nC)=ell^nTor_D C が成立するので、
ell²S⊂ell³Tor_D C となり、長さの比較で b-2<=e-3 を得る。
従って e は b+1 又は b+2 である。
この等式は D の構造定理を使うだけで、R0 加群や環の分裂を仮定しない。

b、e の実値と全四次 initial の完全性は未決である。
元 H0、nil 三乗、全分類への結論はない。
再現は `code/ambient_B18_C_full_quartic_bound_source_review_20261105_v63.py`。
