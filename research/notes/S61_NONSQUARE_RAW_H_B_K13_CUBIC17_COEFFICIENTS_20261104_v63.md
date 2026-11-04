# Ka13上のactual cubic17×24全係数入力

2026-11-04。既受理のsource/cost preflightが指すPOWER-free17×24行列を、actual Bbadの全Ka13上で評価する。新しい点・素数を選ばず、Gaussian・rank・minorも計算しない。

普遍的なflat sourceは17×24、numerator84253項である。Ka13上の1113個のmonomialを、既知の認証済みnode高々三つの積として計算する。その全値はmemory cacheだけに保持し、旧証明書を複写しない。共通整数分母Dを使って各source和をQQ[a]次数12以下の整数線形結合で作る。

全分母c*u（u=1,a,B）に対し、C_ij=(cLCM/c)*(aB/u)*p_ijを計算する。ここでaB/uは小さい既知Ka13値として保持する。C=cLCM*aB*Mなので、既受理のFULLU逆元により元の行列と同じrank/全射性を持つ。既存8段階の単元minorも全Ka13上で単元である。

保存する各source恒等式はc*Cnum*den(P)*den(mult)−cLCM*den(C)*num(P)*num(mult)=g13*Qであり、整数除法の余りが零であることを確認する。Qを重複展開して保存せず、その一意な長除法recipe・degree・canonical係数hashを保存する。Pの全演算子はpinned flat sourceと既存monomial factor recipeで具体的に決まり、未知のcofactorやsolveはない。独立checkerは各行の全sourceから同じ有限係数恒等式を再構成できる。

上限は一回60秒・2GB・出力全体32MiB・明示150000項。input/cacheをartifactとして重複保存しない。matrix値を一成分ずつ保存し、guardに達した場合はprefixだけを保持する。rank17ならfull cubic57を含意するが、このpacketはrankを調べない。rankが不足してもfull57否定へ進まない。

再現: `.venv/bin/python code/s61_nonsquare_raw_H_B_K13_cubic17_coefficients_20261104_v63.py --output work/K13-cubic-coeff-fresh`。新しい係数packetの独立レビューは別途必要である。

実行は4.603秒、最大RSS52916KiBで全408値を保存した。354成分が非零、QQ係数4602個、最大分子2076bit・分母1916bit。全408 source恒等式の整数長除法の余りは零で、商の次数は高々11だった。商環整数演算1358回、live明示項33120以下、input/cache複写なしで**出力全体2952077 bytes**だった。

次の独立したpositive-only判定案は、C:Ka13^24→Ka13^17のrestriction of scalarsである。1,a,...,a^12を基底とし、各C_ijの乗法13×13blockを並べてQQ上の221×312行列を作る。許可された既存素数32003で、g13のleading coefficientと全保存分母が単元であることを先に確認し、221の非零係数minorを得られれば、同じQQ minorが非零なのでQQ-linear写像は全射、従ってKa13-linear写像も全射である。元のfull cubic57への含意はQQ上の既受理単元prefixを通して行うため、そのprefixのGF単元性を新しく仮定しない。modular rankが不足した場合は結論できず、別の点・素数・QQ巨大Gaussianへ自動移行しない。今回はこの判定を未実行である。
