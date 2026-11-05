# ell³ product base39 の独立 SOURCE 検証

2026-11-05。PRIMARY report
`AMBIENT_B18_ELL3_RAW39_SOURCE_PREFLIGHT_20261105_v63/report.json`
SHA256 `03ee57279bd4143098018d1273d1d2161b2bdb12b1b48f0ad2741330268d13c4`
を対象とする。独立 REVIEW report は
`AMBIENT_B18_ELL3_RAW39_SOURCE_PREFLIGHT_20261105_v63_REVIEW/report.json`
SHA256 `f300322e219f1dfad4c7d2996886dedbdec400432717dd4e814670321aeffb33`。
author program の import、商の新しい rank・乗法表・physical source 値の計算はない。


## 同じ実際の DVR 上の証明

ℓ は実際の同一関数 ℓ=5(sA-rA) であり、各枝は既受理の Weierstrass 証明により
(D=\mathbb Q[[\ell]]) 上 free rank3 である。(r=r_A) は両枝の同じ transverse
coordinate として使用する。(J=I_A+I_B=(z_0,\ldots,z_6,(r,\ell)^2)) も同じ
実際の source からの等式である。今回、beta2 の固定3行 `[0,4,5]` を
((r^2,r\ell,\ell^2)) の順に読んだ小行列式も独立に確認し、値は (-5/192\ne0)
となった。これは商の rank search ではない。完全環での Nakayama により
三つの plane 二次単項式が (J) に含まれる。

(U=R/(I_A\cap I_B)) は二枝の direct sum の部分加群なので D-free rank6である。
fiber product の共通商は (k[r,\ell]/(r,\ell)^2)。枝の基底を (1,r,r^2) とすると、
定数係数の差は ℓ²、r係数の差は ℓ で割り切れるので、次が U の実際の D-基底となる。

\[
e_0=(1,1),\quad e_1=(r,r),\quad u=(r^2,0),\quad v=(0,r^2),\quad
p=(\ell^2,0),\quad q=(\ell r,0).
\]

この順の重みを (0,1,2,2,2,2) とする。両枝の Weierstrass 関係の係数は

\[
r^3=A_2r^2+A_1r+A_0,\qquad
\operatorname{ord}_\ell(A_2,A_1,A_0)\ge(1,2,3)
\]

を満たす。例えば (ru=A_2u+(A_1/\ell)q+(A_0/\ell^2)p)、
(rp=\ell q,rq=\ell u) であり、すべて重みを1以上上げる。
B枝側の式では (A_1e_1,A_0e_0) などが追加されるが、同じ下界を保つ。
各 (z_i) は A枝で0、B枝で order2以上の (-\beta_i(r,\ell)) なので、
重みを2以上上げる。これらの全六列の形式的 valuation 下界を独立コードで確認した。
最大基底重みが2であるため

\[
\mathfrak m^5U\subset\ell^3U.
\]

(B=R/(I_AI_B))、(T=(I_A\cap I_B)/(I_AI_B)) とすると (JT=0)、従って
ℓ²T=0である。U が D-free なので、(0\to T\to B\to U\to0) を
(D/(\ell^n)) と tensor しても左端は単射である。n≥2では、B₃→B₂ の制限は
同じ T への恒等写像となる。よって ℓ³商の ᵐ⁵ は U₃ で0、従って T に入り、
B₂での像も0ならもとの元も0である。

既存 B₂=raw33 の9個の座標作用を使用し、次数5の可換単項式1287個すべてが
単位元を0へ送ることを独立に確認した。既受理の可換・巡回作用から ᵐ⁵B₂=0であり、
従って **実際の ᵐ⁵B₃=0** が従う。新しい quotient table は構成していない。
既存の全作用の分母 LCD と整数化した列 L1 から、この固定検査の有理数高さは
187 bit以下とも独立に評価した。`degree5-height-and-metric-scope.json` に記録した。
凍結 evidence の `source_arithmetic_observed_bits=15` は beta/g/plane 配列だけの値であり、
全作用積の高さと取り違えない。

## Exact cubic truncation と有限入力

三次 Taylor pair も同じ J、同じ二枝の D-free3、既受理の同じ B₂ を持つので、
同じ証明から ᵐ⁵=0である。generator の省略誤差は order4以上、product 誤差は
order5以上であり、actual・truncated の両方の ℓ³商で0となる。
両方向の ideal 包含により

\[
(I_AI_B,\ell^3)=(I_{A,\le3}I_{B,\le3},\ell^3)
\]

が厳密に成り立つ。z依存 beta3 を既受理の (I+M_2) unit generator changeで除く際に
生じる plane order4項も、この product の等式では許される。

独立な Fraction 二変数合成で、s=r+ℓ/5、CBの実際の一次座標変換、beta/gの全項を
再構成した。75 ambient monomial の別列挙から、21 pair と14g の全35関係、
必要な plane 倍62非零行・204項のすべてが PRIMARY と一致した。
49 monic (z_i(z_j+\beta_j)) で多重z項を消し、pair関係で対称性を保つ。
各関係への z作用はこの plane module の関係内に入り、gAgBは次数6なので0である。

長さ20・33と free rank6、ℓ²T=0から、Tの指数1・2の個数はそれぞれ7。
従って B₃の長さは (18+21=39)、raw cokernelは二つの枝の9+9を加えて57である。
これは新しい行列 rank の計算結果ではなく、実際の source と既受理の module 型からの結論である。

個々の branch quotient には四次係数が寄与し得るため、枝自身を cubicで切ってよいとは
主張しない。Cprime・genuine torsion・比較kernel・nil idealへの追加結論もない。
