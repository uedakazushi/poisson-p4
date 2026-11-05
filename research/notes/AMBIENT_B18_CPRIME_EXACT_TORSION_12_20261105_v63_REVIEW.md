# Cprime 指数(1,2)・全 Hilbert 関数の独立純代数監査

2026-11-05。source/value を新たに計算せず、受理済みの同じ実 R0-module Cprime を M、D=QQ[[ell]] と書く。ell=5(sA−rA)で、cut変数rAとは異なる。D上有限で rank2、torsion指数(1,2)、実HF0..2=(3,3,3)、`m³(M/ell³M)=0` を使用する。

まず `m³M⊂ell³M`、また ell∈mから逆包含もあり、`m³M=ell³M`。従って

`ell³M ⊂ ell*m²M ⊂ m³M = ell³M`

となり、`m³M=ell*m²M` を得る。ell はM上の非零因子とは仮定せず、一度も消去・除算しない。これにmの冪を掛けるだけで、j>=2に対し `m^{j+1}M=ell*m^jM`。最初の等号を基底とする帰納法でj>=3なら `m^jM=ell^jM` である。

非標準のD加群型 `D²⊕D/ell⊕D/ell²` により、j>=2でell^jMはD-free rank2となる。ゆえに実際のm-adic全HFは `(3,3,3,2,2,...)`、seriesは `(3−t³)/(1−t)`、multiplicityは2である。M全体でm³が零になるという主張ではない。D torsion長3はCprimeに限るもので、元H0(D0)長やnilidealの指数ではない。

標準写像 `M[ell] -> ker(ell:M/ell³M)/im(ell²)` はR0-linearである。D構造定理で調べると、free二本の有限商kernelはちょうどimell²として除かれ、指数1,2のtorsionはell³Mと交わらずそのまま残る。従ってこの標準写像は二次元の真正targetへの同型であり、選んだD分解には依存しない。実基底・作用行列・connecting値は別作業である。

追加の純粋な帰結として、k>=3なら `length(M/m^kM)=2k+3`。また `ell²M⊂m²M` と既知の長さ7および6から、`length(m²M/ell²M)=1`。この商はm³M⊂ell²Mによりmに消されるのでkに同型である。これもMだけの結論である。
