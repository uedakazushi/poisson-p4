# 選択したθ/H二wordのsource展開：整数scaleとvector recursionの事前検査

2026-11-04。全Inv₂₄基底に必要となったθ⁵Hとθ⁷HはHに関して一次であるため、16×16のHやその高い冪をsource上で密に展開する必要はない。このpacketは既存定数行列と保存済み形式平面係数の検査だけを行い、新しいsource不変式の係数は作らない。

旧整数source座標をb=1080BHβとする。Ωと全ρ(Aⱼ)は整数、2Gram⁻¹も整数。既存の全136 quadratic coefficientについてμᵀ=−C+2Dと2μの整数性が確定している。μ=JSならSは対称で、その成分Sᵢⱼが`sp_basis()`の座標だから、Θ=2ρ(μ)=2θも整数quadratic行列である。

\[
 Z=(\rho(A_j)b)_j,\quad V=\Omega Z,\quad G=2\operatorname{Gram}^{-1}
\]

とすると、k=5,7について

\[
 \operatorname{tr}(\theta^kH)
 =\frac{\operatorname{tr}(\Theta^{k+1})+
 \operatorname{tr}_{10}(G V^T\Theta^k Z)}{2^{k+1}}.
\]

従ってZ₀=Z、Zₖ₊₁=ΘZₖを16×10の二つの生存levelだけで計算すればよい。7levelはscalar polynomial product高々17920回。dense16×16 H word cacheを持たない。

純trace項はμ̂=2μのs₁=tr(μ̂²)/2、s₂=det μ̂により

\[
 \operatorname{tr}\rho(\widehat\mu)^6=264s_1^3+408s_1s_2,
\]
\[
 \operatorname{tr}\rho(\widehat\mu)^8
 =1032s_1^4+3488s_1^2s_2-4208s_2^2.
\]

これらはWの16個のweightの6乗・8乗の和をQQ[a,b]で展開して確認した。regular semisimple元の稠密性により普遍的恒等式である。さらに同じ保存済み形式平面のF₄とpᵢを用い、s₁=F₄/2、s₂=p₀−2p₁+4p₂−8p₃+16p₄から、既存のθのみのtrace6・trace8の**全係数**にmod32003で一致することも確認した。新しい点は使っていない。

選択した二つのC/D wordと二つの旧wordについて、scale2による精密化後のraw JSON出力の厳密な上界は次の通り。

|word|項数上界|係数bit上界|raw bytes上界|
|---|---:|---:|---:|
|C/D 000101|27418|468|4962658|
|C/D 000111|27418|467|4962658|
|θ/H 000001|27418|480|5072330|
|θ/H 00000001|80125|639|18669125|

全上界33,666,771 bytesは依然32MiBを112,339 bytes超える。従って後続の一回のsource実行には60秒/2GB/出力100k項/中間150k項/全raw32MiBのcapを設け、完成した出力だけを逐次保存する方針とした。これは実時間・memory・保存量の完走保証ではない。

保存は`code/q_new_family_hessian_source_vector_preflight_20261104_v63.py`と
`reports/Q_NEW_FAMILY_HESSIAN_SOURCE_VECTOR_PREFLIGHT_20261104_v63/report.json`。既存係数・定数行列・character weight・source numeratorのhash、全検査、上界を含む。約0.62秒。別予算による実際のsource展開は`Q_NEW_FAMILY_HESSIAN_LOW_INVARIANT_EXPANSION_20261104_v63`に保存する。
