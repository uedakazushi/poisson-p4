# 固定二tailの有限体Kronecker gcd：標数零の共通因子は未決

独立に整数掛け戻しを済ませた N0（4500項）と N2（2317項）を対象に、a=t^646、e=t^19、d=t とし、固定素数32003で一回だけ univariate gcd を求めた。e次数33以下、d次数18以下の全箱で指数符号化は単射となる。入力の重み次数10544と7975は還元後も保たれた。

二本の uint32 vector と4096の metadata だけを保持する Euclidean worker を新規に用いた。係数積は30 bit以下、全保持量25000、45秒 soft/60秒 hard・2GB・RAW4MBに固定した。全67641091積・4077 Euclidean stepで完了した。

結果は次数3899・94項の monic gcd である。gcd が一なら Gauss と leading degree 保存から標数零での非定数共通因子を排除できるが、今回はその仮定を満たさない。従って coprimality、空補集合、全FULLU57/239、実bad Fittingの等号は結論しない。これは数値的なbeta/source点の標本化ではなく固定有限体への係数還元と単射単項式符号化であるが、正の gcd を標数零の因子へ読み替えてはならない。

次の候補は94項の有限体支持を参考にした bounded CRT lift と元二整数多項式への全係数除算・掛け戻しである。候補が得られても正確な恒等式なしには採録しない。gcd の leading coefficient は元二式の leading coefficient の整数gcdで scale し、最後に content を除く必要がある。単一の追加素数を無計画に繰返さず、固定有限個のprime・全保持量・整数増大・prefix保存を先に証明する。まだ lift code や次の値 GO はない。

`reports/S61_TWO_TAIL_FIXED_KRONECKER_SOURCE_20261105_v63` と `reports/S61_TWO_TAIL_FIXED_KRONECKER_ENTRY_20261105_v63` に code、compile/help、算術前正予算、monic gcd、原入力SHAを保存した。この診断の独立再計算はサービス容量制約により未実施。
