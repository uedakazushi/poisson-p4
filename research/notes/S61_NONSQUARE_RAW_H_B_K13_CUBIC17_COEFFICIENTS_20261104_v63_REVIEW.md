# Ka13 cubic17×24 全408係数の独立 source 監査

元の25×32 unit checkpoint から、残る17 source rows と24 target columnsを直接読み、
全84253 numerator coefficients と denominator=c*(1,a,B)を独立照合する。そこから作る
flat input は旧 P80 source17x24-input.txt と全byte一致し、P100の同名fileとも同一である。
これは (a,e,d) の普遍入力であり、P80/P100曲線の zeta 特殊化済み行列ではない。

独立受理済み actual Ka13 source nodes から、全1113 monomial factor recipes の exponent
と値を別実装の fixed-table QQ 演算で作る。全408の source 和、単元 clearing、保存値と、
保存された整数 g13 quotient の一意な係数hashを長除法で確認する。作者の producer や
Ring、逆元探索、rank、点、素数、GB は実行しない。

行列は正確に C=cLCM*a*B*M で、cLCMは非零QQ定数、a,Bは受理済み FULLU inverse に
より全Ka13で単元。既存8prefixのdeterminant−a⁴FB²C²/4608も単元なので、この実行列の
rank17は actual full cubic57 の十分条件となる。ただし本 source gate だけでは rank を
判断せず、後段の有限線形代数証明を別に監査する。結果と入力hashはreportに保存する。
