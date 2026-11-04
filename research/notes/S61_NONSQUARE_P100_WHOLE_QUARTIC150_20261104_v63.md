# 許可された P100 全点での四次階数150

本 packet は二つの独立レビューをゲートとする論理的結合である。行列式・係数・normal form・階数・点の新しい計算は行わない。元の有限計算、五つの選択小行列式、及び全ての凍結資料を変更しない。

T=T100=QQ[z,1/(z(z+20)(z+60)(203z+2700)(z+12))] とする。M は全実319×239四次写像である。受理済み NF117、H2 の単元一方向、H1 の全23方向の section を用いた有限演算は T 上で大域的に定義され、coker(M)=coker(Gamma319×98) を与える。全 pure tail と元319源を保持する。

既存の五つの実9小行列式のイデアルは厳密に(H)であり、H はモニック十五次である。ここでは元小行列式と保存された primitive 多項式の違いが T の既知単元だけであるという source binding も用いる。この主張は Conn の `P100_K7W_FINITE15_RESCUE_REVIEW` において、旧五つの H による除算と cofactor の既存素数での Bezout1 を含め独立受理された。

別の実有限代数計算は Hbar の全七因子で Gamma の階数98を証明する。元の固定選択行列、全必要分母の逆元積、chi と H1 section の恒等式を含め Mark の独立レビュー受理を必要とする。モニック有限自由 Z_(32003)[z]/H と有限 Nakayama により、標数零の H=0 の全有限候補で階数239となる。特殊ファイバーの一部や一つの点での推測ではない。

## イデアルとしての結合

源・標的の冗長な表示を含め Fitting イデアルは cokernel のみに依存するので

    Fitt_89(coker M) = I_150(M) = I_9(Gamma).

五つの実小行列式により H は I_9(Gamma) に入る。一方、Gamma modH の全射から

    I_98(Gamma)+(H)=T

を得る。大きい小行列式のイデアルは小さい方に含まれるから I_98(Gamma)⊂I_9(Gamma) であり、従って I_9(Gamma)=T、すなわち I_150(M)=T である。

ゆえに許可された P100 曲線の**全複素点で実四次階数150以上**。H=0 ではより強い239が成立する。単一の固定150小行列式や新しい展開済み Bezout 係数を計算したわけではない。

**全 P100 で239、全 FULLU で三次57・四次239、一般アフィン尽くし、全成分数上界はいずれも主張しない。**

## 保存物と再現

結合コードは `code/s61_nonsquare_P100_whole_quartic150_20261104_v63.py`。同名 report の phase は両レビューを待つ状態から先に保存し、両方の FINAL SHA256 と all_passed が揃ってから結合証明書を凍結する。報告書は旧一次報告、新一次報告、両独立レビュー、全源の型付き演算を固定する。

    .venv/bin/python code/s61_nonsquare_P100_whole_quartic150_20261104_v63.py \
      --output reports/S61_NONSQUARE_P100_WHOLE_QUARTIC150_20261104_v63 \
      --finite-peer PATH_TO_FINAL_MARK_REPORT --finite-peer-sha FINAL_SHA256

この再現は既存の保存値・SHA と論理ゲートを読むだけで、数値計算を再実行しない。
