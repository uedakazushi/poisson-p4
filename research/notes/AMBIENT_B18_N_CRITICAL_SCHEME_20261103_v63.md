# generic B18でのN正規化射：全critical schemeはline＋squarezero junction

2026-11-03。CAと同じB18代表

    P0=∂y ∧ (x²∂h+h*x∂z+(h²+z²)∂w)

における、Nのnormalizer source上の全ramification Fitting schemeを計算する。
target Nや実CA/CB交叉Dの完成環を計算したという意味ではない。
proper旗fiberと有限正規化近傍には受理済みB18 quadratic-image reviewを用いる。

固定旗19次元空間はCAの20空間からbinary fの第三係数を除いたものである。
同じ16次元parabolicによる全19基底の安定性を全150係数で確認する。
P0でのparabolic fiber orbitはrank10であり、19係数のうち10個への非零minorを保存する。
従って旗9＋parabolic10の計19個のsmooth source因子をetaleに除き、9変数断面を得る。

固定fiber19列を一定の単元blockで除いた全150×9旗matrixをQとする。
行5,15,22,23,31,35,50,65、列0,1,2,4,5,6,7,8のblockは断面上で恒等的にdet1である。
残る列3の全Schur成分が、最大minorのFitting idealと正確に同じidealを生成する。

断面座標をt0,...,t8とすると、全64非零Schur式のidealは

    (t0−5t6/2, t2−11t6/2, t3,t4,t5,t7,t8, t1*t6,t6²)

に等しい。両ideal包含をQQ Gröbner消去で検査し、別のMacaulay2計算でideal等号も確認する。
w=t1, v=t6とおくと、全cone critical完成環は

    C[[g1,...,g19,w,v]]/(wv,v²)

である。projective版のsmooth変数は18個となる。
CA critical schemeを単にbinary係数0で切って得た結論ではなく、N自身の全fiber19微分から再計算した。

primary分解は

    (wv,v²)=(v)∩(w,v²).

従って被約critical locusは滑らかな一本のlineであり、junctionに非零squarezero nilpotentが一つある。
関連素は極小素とembedded素の正確に二つである。
cone critical環はdim20/edim21/depth19、projective版はdim19/edim20/depth18、重複度は1となる。
前稿のN-source ramification cotangent rank7=28−21と整合する。

これはnormalizer sourceの相対微分のFitting schemeであり、
targetの全odd ideal、実Dのnilradical・associated primes・depthの結論は含まない。

* code: `code/ambient_B18_N_critical_scheme_20261103_v63.py`
* evidence: `reports/AMBIENT_B18_N_CRITICAL_SCHEME_20261103_v63/`
