# B 次数だけによる G₃₂ edge 候補の除外

2026-11-04。保存済みの低次数 source 係数だけを読み、G₃₂ が NEW25 上で恒等零でないことを証明した。**source32多項式の展開、source polynomial 同士の乗算、点の評価は一切行わない。**

分母多項式を H=147456Y−(3072W−1)²、spectral invariant を D₀=p₀ と書き分ける。低次数 source map の B データは次であり、全係数ファイルから再確認した。

|量|numerator B次数|分母B指数|分母H指数|
|---|---:|---:|---:|
|S=F₄/20|0|2|4|
|D₀=p₀|0|3|8|
|p₂,p₃|1|4|8|
|N₁₂|1|5|12|

各式はさらに非零有理整数定数を分母に持つが、係数体 QQ(a,Y,W) 上で B 次数を考える際には単元である。

\[
 K_{16}=D_0(p_2-S^2+p_3-D_0)+(p_3+S^2)^2.
\]

これを B⁸H¹⁶で払った numerator の B 次数は高々2。各項 D₀p₂,D₀S²,D₀p₃,D₀²,p₃²,p₃S²,S⁴ の次数を個別に数えて確認する。

\[
 G_{32}=D_0N_{12}^2-4SK_{16}N_{12}+4K_{16}^2
\]

を B¹⁶H³²で払うと、三つの項の B 次数は順に

\[
 5,\quad\le4,\quad\le4.
\]

最初の項は元の分母 B¹³H³²から追加の B³を持ち、N₁₂ numerator が正確に B 一次なので、B⁵の係数は

\[
 \frac{\operatorname{num}(p_0)
       \operatorname{lead}_B(\operatorname{num}N_{12})^2}
      {d_{p_0}d_N^2}\ne0.
\]

保存した p₀ numerator と N₁₂ の B¹項の係数列は両方非零。係数体は整域だから、この積が非零であることには積の展開が不要である。他の二項は B⁵を打ち消せない。従って G₃₂ は source 上で非零で、**NEW25 方程式の候補として除外**される。

普遍的な ordered-chart 恒等式

\[
 G_{32}=64p_4D_0^2(gj)(hi)
\]

と、NEW25 の spectral base への優越性を合わせる。NEW25 は既約で、その一般点では D₀p₄≠0かつ G₃₂≠0。よって既知の b=e=1 chart の一般 NEW25 点では

\[
 ghij\ne0.
\]

従って R=g/h,T=gh を使う Laurent chart は一般 NEW25 に本当に交わる。ここで gj,hi は cross block の積で、元の Poisson family の source u,v とは別である。

これによって、八次表示の別の分母 L が NEW25 上で単元であることは分からない。L 単元の従来の一般係数体の証明には m の base 上の超越性を使い、NEW25 へ同じ議論を移せない。NEW25 multisection の次数も、定義方程式の次数下界も、この結果から新たに決めない。

保存：`code/q_new_family_edge_boundary_valuation_20261104_v63.py`、
`reports/Q_NEW_FAMILY_EDGE_BOUNDARY_VALUATION_20261104_v63/{report.json,certificate.json.gz}`。

証明書は全B次数・分母指数、K₁₆各項の上界、唯一の最高項の二つの保存済み因子を保持する。N₁₂側は元ファイルの B¹項のリストとして保存し、その B を一つ除いたものが leading coefficient である。低次数 map の独立レビューと普遍 edge 恒等式の hash も固定した。実行約0.014秒。
