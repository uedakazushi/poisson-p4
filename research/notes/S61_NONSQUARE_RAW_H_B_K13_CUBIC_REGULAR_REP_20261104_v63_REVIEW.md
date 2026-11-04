# Actual Ka13 bad locus 上の full cubic57: source408と固定minorの独立監査

先に独立検証した全408 source receipt と、固定素数32003の regular representation を結ぶ。
全分母408個と g13 leading coefficient の p-unit 性を自分で確認し、保存 modular input の
全5304係数を QQ source から再構成する。新しい点・素数・rank search は使わない。

別実装は monic g13 による a^0,...,a^24 の表から、指定221列の係数行列を直接作る。
列(j,k)は C_ij*a^k の係数、行は(i,l)であり、aを根に評価していない。指定 row permutation
に対する全48841の P*A=L*U を別計算し、対角積を確認する。保存 LU を解放後、異なる
FLINT fixed determinant algorithm で221×221 minor の値30542≠0も確認する。

V=Z_(32003) 上で A_V=V[a]/monic(g13) は free13、全 C_ij は p-integral。非零 modular
minor は同じ QQ minor の非零性を証明する。従って C:Ka13^24→Ka13^17 は QQ-linearにも
Ka13-linearにも全射である。任意の residue field への base changeも全射なので、全成分で
rank17となる。体性、既約性、reducedness、étale性は仮定しない。

characteristic0の source bridge は、普遍的 pre-curve17×24入力と actual25×32 checkpoint の
全byte・全84253係数の一致、C=cLCM*aB*M、全FULLU既知単元8prefix、product41、実 Cdelta
source quotient criterion を保持する。したがって actual Ka13 bad locus 全体で full cubic57
が成り立つ。元の GF cubic/prefix が単元という新しい仮定は置かない。

この結論を全FULLU cubic57、quartic239、P100の149下界、generic Poisson exhaustionへ
広げない。証拠・範囲・入力hashは同stem reportとindependent-proofに保存する。
