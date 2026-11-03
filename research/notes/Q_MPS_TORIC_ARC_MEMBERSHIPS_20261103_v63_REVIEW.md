# 全 Qa の actual T02/T03 帰属の独立監査

2026-11-03。新規 producer と全保存資料は時間以外一致した。
別の SymPy 行列・線形場 wedge から generic diag(0,1,2,x,y) arc の全150係数を再構築した。
Qa 極限、Euler-wedge 差、単位下三角 eigenframe を厳密に検証した。
可換 tracefree 行列冪からなるため全パラメータで actual canonical Poisson arc である。

4×4 skew 行列の線形 Pfaffian dual は A D(A)=−Pf(A)I、D²=I を満たす。
保存 MPS の −Pf(B)B⁻¹ の規約と正確に一致する。
全120置換を六 basis で別に構築し、新 chart の bivector dual と置換 biresidue の符号を照合した。
bivector係数へ直接 period 条件を課す取り違えはない。

両曲線について x に関する primitive 性を確認し、次数を保つ法5特殊化の既約性を
producer の Rabin とは別の Berlekamp 分解でも確認した。
Gauss と preserved x degree により QQ[x,y] 既約性が従う。
実 period basis と保存 period parameters の全六成分差の分子は f で割り切れた。
distinct eigenvalue、rho の分母・零・±1、接微分、rho²+1 に対応する各剰余は非零である。
QQ 曲線環が整域なので積も非零であり、複素化した少なくとも一つの幾何成分上に
同時に admissible で δa≠0 となる非空開がある。絶対既約性を仮定していない。

その各点で t≠0 arc は actual T02/T03 の toric boundary period 族に写る。
toric period 族自体は Pfaffian dual による多項式族なので、別の smoothing denominator が
落ちる period も閉包として同じ成分に入る。
t=0 の Qa と、非定数 a の稠密像に対する閉性から、全有限複素 a の帰属が従う。
有限体標本から複素 membership を外挿した証明ではなく、mod5 は既約性だけに使う。

固定・一固有値 scan はこの方法の正の証人の記録に限り、他成分への非所属ではない。
旧 Qa1∈T37 は元の P(a)+tτ(a) の actual MPS identity と τ(1)=Qa1 の射影極限に由来すると確認した。
旧大計算の再実行、新成分数、全通過成分、局所枝の尽くしは主張しない。
監査結果: `reports/Q_MPS_TORIC_ARC_MEMBERSHIPS_20261103_v63_REVIEW.json`。
