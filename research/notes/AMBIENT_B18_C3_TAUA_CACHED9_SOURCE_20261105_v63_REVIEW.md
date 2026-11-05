# 実 tauA 九射影の独立 SOURCE・prepared worker 監査

2026-11-05。SOURCE PRIMARY `AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63/report.json`
SHA256 `8234e80934cba6d761eb72666d3d4337fcfa07a1c5d12ca9b0776f85d9339bd0`
に対する REVIEW report は
`reports/AMBIENT_B18_C3_TAUA_CACHED9_SOURCE_20261105_v63_REVIEW/report.json`
SHA256 `b7cd179090ff8a22e20de1fa007ffbe3d54f25bd7168bd374dd6d14f8bbf594b`。
SOURCE検査では九個の三次dotをpairしていない。

同じ実際のCとR0、受理済みfirst21/quad-source frameを使用する。
(h_{li}=(b_{B,l})_{z_i}) とすると、実際の source

\[
q_l=F_le_A+\sum_i h_{li}R_{H,i}
\]

は一次項が0、eA成分が厳密にdelta_l、plane二次成分が−bB_l,2 eHである。
3×12の bB2 map、固定3×3逆行列の両積、全27 kernel恒等式、
free-tail I9を持つ12×9行列Kを独立に確認した。
固定minor5/768がunitなので、pivot成分は−P逆*Bfreeで一意であり、これは全九次元kernelである。

Kの各列でqを組み合わせると二次plane部分が消える。
残るzを含む二次成分は実際のfirst21で消せる。
**eAの消去には厳密な列zi eAしか使用しない**ため、高次eA carryはない。
RBのbeta2やRHのdelta/bAの高次部分は残るが、eB/eHにしか入らず、
今回のeA射影はそれらを消す。このためcorrected三次のeA成分は正確に

\[
\sum_lK_l\delta_{l,3}
\]

である。未知の他成分を零に置いたわけではない。

元gA3をs=r+ℓ/5で別実装の二項展開により合成し、係数

\[
g_{A,3}=\frac{125}{24}r^3+\frac53r^2\ell+\frac16r\ell^2+\frac2{375}\ell^3
\]

を確認した。(gA3,ℓ³)を法とするeAのr²ℓ係数functionalは

\[
\varphi_A(\delta)=\delta_{(2,1)}-\frac8{25}\delta_{(3,0)}.
\]

raw33の全74 normal formについて、この二成分の148個の値がliteral Taylor coefficient
そのものであることを確認した。元多重z単項式のmonic reductionはorderを上げつつzを残し、
filtered normal formは高次項からpureplane三次を作らない。
従って受理済みC2 cacheのindex30と29をそのまま読める。
未保存のrℓ²係数やgauge3に対するこのfunctionalの係数は0であり、値を0と仮定していない。

九個のunpaired adjoint全18項をKと8/25から再構成して一致を確認した。
24 leafの高さだけを読み、weight最大5 bit、leaf最大18 bit、各dot最大2項から48 bit上界を得た。
prepared codeは `5891d81746e4eb9ad9df0e5d93879c0d8246cadfc496e03a17e79ccf42249888`。
parse、積、和の各Fraction演算の前に2048-bitの上界を検査する。
入力＋全出力のowned envelope18000/hard30000、RAW400000 bytes、60秒、2GiBを確認した。
最初のoverflow fixture harnessが自分のassertionを同じexceptで捕捉した履歴は、
別companion `f8eba09d…` に保存されている。worker自体は変更していない。
レビュー側は元のguard式と演算前順序も直接確認した。

非零値が出た場合は、plain length21 modelの外に実際の三次関係が一つ出る。
さらにℓを掛けると、eAのCI(3,3)商のr²ℓ²係数が同じ非零値になる。
r*gA3、ℓ*gA3はr⁴,r³ℓ上のtriangular unit blockを持つので、r²ℓ²は残る。
従って次数の異なる二関係からC3長さの上界19が従う。
これだけで上界18や残り二torsion指数両方2とはしない。
全九値が零なら、この九射影だけのnegative結果である。
