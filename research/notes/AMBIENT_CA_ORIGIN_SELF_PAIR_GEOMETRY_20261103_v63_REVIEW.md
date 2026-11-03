# CA実self-pairの全二branchと非CM性：独立監査

2026-11-03。対象は個別の **reduced CA** であり、full Poisson scheme、CA/CBの交わり、原点のfull odd idealではない。

著者の係数Lieルーチンを使用せず、raw CA20族からdivergenceを除き、実5×5テンソルの微分

    L_E P = E P + P E^T − Σ (EX)_i ∂_i P

で全GL25方向と両flag9方向を再構成した。QQ(a)で両sourceはrank29であり、共通fixed20を除いたminorは
−11943936 a^9/15625、864 a^14/15625。両方とも全a≠0でimmersionである。
第二sourceの **全29方向** を実行列g⁻¹で同じP(a)にpushbackし、全58列のrank37を独立に確認した。
著者の共通fixed20を除く150×29行列も全係数で一致し、minorは4076863488 a^25/3125だった。
GL軌道rank19、scalar追加rank19、dP/da追加rank20も厳密に一致する。
従って射影軌道は18次元、この族のG-saturationは射影19次元である。

## 全fiberから実completionへ

既に独立監査された全proper CA flag-sourceと、a=0での**全flag fiberのlength2**を依存として使用する。
properかつそのfiberが有限なので、a=0を含むtarget openで全source mapはfiniteになる。
有限moduleのfiber lengthの上半連続性により、さらに縮めて全fiber length≤2となる。
多項式arcへの逆像はa=0を含むopenであり、そのa≠0の点には、先行監査で全次数を確認した二つの異なる旗がある。
従ってそのopen arcでは全fiberが正確に二つの被約点である。
この議論は全a≠0への一様なfiber主張ではない。

有限sourceの代数は、実reduced image CAの局所環を含む。completionはflatなので、その単射は
二つのsource completionの積B₁×B₂への単射として保たれる。
sourceはsmooth29であり、上のrank29によりambient形式環T→Bᵢはsurjective、Bᵢ≅C[[29]]となる。
二つの接空間の和の次元37>29によりkernel I₁,I₂は異なる。
全fiberの二点を用いたので

    Â_CA = T/(I₁∩I₂),   T/I₁≅T/I₂≅C[[29]]

であり、第三のformal branchはない。global CAの既約性は変わらない。
射影ではsmooth28の二つのformal branchになる。

## meeting schemeとdepth

実scheme intersection C=T/(I₁+I₂)の接空間はT₁∩T₂なので、そのembedding dimensionは正確に21である。
これは**whole CAのembedding dimensionが37に等しい**という主張ではない。
全GL×parameterのpaired-sourceが両branchに入る実20次元族を与えるので、20≤dim C≤21。
clean intersection、meeting ideal、conductorの式はいずれも推論しない。

完全列

    0 → Â_CA → B₁⊕B₂ → C → 0

においてdepth(B₁⊕B₂)=29、depth C≤edim C=21である。
depth lemmaからdepth Â_CA=depth C+1≤22<29。
従って、このopen arcとそのG-saturation上で個別のreduced CAはnonnormal・singular・非CMである。
射影のdepth上限は21、次元28である。
この結果を原点completion、全bad locus、R19/Nへのmembership、full Poisson schemeへ拡張しない。

## 再現と保存

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CA_origin_self_pair_geometry_review_20261103_v63.py --output reports/AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63_REVIEW.json --replay work/v62-review/CA-origin-self-pair-geometry-replay

著者のfresh certificateと保存certificateのbyte一致、全係数行列、QQ(a) rankと厳密minorを検証する。
全propernessとa=0全fiber length2は既存の独立監査を明示的に再使用する。
GBは再実行しない。公開証拠を変更せず、code/note/reportを別保存する。
