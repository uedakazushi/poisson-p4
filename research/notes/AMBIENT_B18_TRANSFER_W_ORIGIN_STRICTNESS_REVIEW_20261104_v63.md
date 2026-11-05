# transfer 原点像と真の商の独立検討

2026-11-04。新しい512 source の再評価ではなく、元の式の次数から原点像を検討する。
同じR0原点のC4は k[tauA,tauB]/(tauA²,tauB²) である。proper座標はその極大idealへ
写る。正規化nil二生成元 dh、dcl はproper次数二・三なので、それらの四basis倍pの
tauB単独成分は原点で零となる。任意A0について even_C(A0 tauA)=0 も従う。

transfer44の原式では全項に beta、b、bA、delta のいずれかが入り、これらは
R0原点で零である。従って全512のz44原点は零。保持される三つのwedge定数成分に
関わるのは m0 の行7だけで、これは(5/2)e20である。原式

    z20=even_C[beta*A_j5*(p_B+tauA*p_AB)+indicator*beta*b_j*e_gamma]

の一次係数も、前段のp_Bとeven_C(A0 tauA)の零性から零になる。他の微分項は
beta0またはb0を因子に持つ。Koszul cycleの一次恒等式により対応する反対側のr係数も
零であり、列ごとの元形式函数Euler持上げの保持wedge定数は零となる。top二成分は
z44原点から零。section変更がその選んだliftを固定することを確認した原証明と、
既受理のminimal projectionを合成すると、全512のminimal-origin像は零である。
元函数を有限jetへ先に落としてEulerを適用してはいない。

U=Tnil、A=N/ellN、B=K/ellK、C=image(W→A) とする。正確な列
0→A→U→Cprime[ell]→0 と右商の長さ二、mu(U)=5から、U/mU内のAの像は
少なくとも三次元。独立検証済みの自然Told→Uのorigin rankは二で、Bの像は
その像に含まれる。Cの原点像は零である。従ってAの
U/(mU+image(Told))への像は非零であり、B+C+mAはその核へ入る。
これによりP/(psi(W)+ellP)は非零、psi(W)はPの真の部分module、K+WはNの
真の部分moduleとなる。未知のK/ellK基底や、P/ellPの正規化fiberへの単射を仮定しない。

これはW=NやN=K+Wを用いた三乗零の経路を排除するが、N³の零・非零は判定しない。
全次数512写像そのもののscalar係数再構成、P/psi(W)の全長・階数、全nil生成元も
この検討の対象外である。
