# ell-colon filtration と実cubic11×16からの小四次criterion

2026-11-04。FULLU=D(adFBCJ)で同じTn frameを使う。これは全点四次全射の証明ではない。新しい入口はquarticを展開せずに作れる、実cubic11×16のmapである。

## 1. 16と48のexact sequence: CIやunmixednessを仮定しない

S=R[x0,…,x4,Z,W]、J=(q0,q1,q3,q5)、I=J+(ell Z,ell W)とする。ell_x1がFULLU unitなのでellを一つのlight coordinateにできる。raw四quad mod ellのquadratic-cofactor84×126 mapがrank78となる開U78を取る。そのkernelは六Koszulだけである。この事実からdegree3だけのcolon

    (J:ell)_3=J_3

が従う。実際、ell*f∈J4のcofactorをmod ellに取り、六Koszulを元四quad間のKoszulで引くと全cofactorがellで割れ、cancelできる。base上でもKoszul6は四quadの既存split dualからsubbundle、rank78 kernelもrank6 bundleなのでNakayamaで同じkernelとなる。generic equalityをtorsionなしと仮定したものではない。

従って(I:ell)_3=J3+(Z,W)S2。degree4の標準colon exact sequenceは

    0→V16=(S/(J,Z,W))3 --ell→B64=(S/I)4
                                      →V48=(S/(I,ell))4→0。

V16のdimension16は **実際の五変数35−19** である。元four quadsの20 light productsは全aC-openでrank19、一つのliteral linear syzygyを持つ。naive four-quad CIのHF3=15を使わない。V48のdimension48は126−78。左右はconstant-rank quotient bundlesであり、middleもrank64 bundleとなる。ellがS/J全次数で非零因子、CI、unmixed、CMなどを仮定していない。U78の補集合が既存証明で有限ということと、全有限点がこの開に入ることも区別する。

## 2. 122 sourceにあるcanonical11の像

heavy7開でf_i=H*x_i−p_i、pure r1,…,r9を既存と同じ可逆変換で取る。fullquarticは122 pure rows→B64で表され、全32 genuine quartic residualを残す。

そのsource内でnine combinations ell*r_jとtwo combinations H*(ell Z),H*(ell W)はV48へのprojectionで零である。V16への像はexactに

    r1|ZW,…,r9|ZW,p_Z|ZW,p_W|ZW。

後二式はH Z=p_Z、H W=p_Wという既存heavy rewriteから従い、未知式を新たに零にしない。sourceの11 vectors自体は、ell_x1 unitとsixquad prime-basis変換det−80からsplit独立である。**そのtarget像が独立という主張とは違う。**

別の構成では、16 genuine cubicをZ=W=0へ制限し、残るH*五light projectionのkernel11を取れば同じmapが得られる。product quadricの19方向を35 pure cubicから割れば実11×16となる。これは32 quarticsや282 matrixを使わず、保存57×91 cubic sourceだけから定義される。

## 3. 59×16と、追加rank11 chart上の48×5

122→V48が全射の開では、そのkernelはrank74。既知15 fullkernel方向はこの中でsplit独立なので、quotient source59→V16が元fullquartic cokernelを保つ。これはD43のdual presentationを16 target・59 sourceへ戻すexact criterionである。

さらに上のcubic11→V16がrank11なら、source/targetからその11方向を可逆消去できる。既知15 fullkernel方向とは像がzero/independentなので交わり零。残りはsource96→target53で、V48へのprojectionが全射なら最終mapはsource48→target5となる。generic kernel rank43という数も48−5と一致する。

この5×48の全射は追加chart上でfull239と同値である。rank11が落ちたlocusをfullquartic defectとみなさず、逆に除外もしない。122→V48 projectionが不全射のlocusも別に残す。heavy7欠損、U78外の有限点、これら二gateを全て保つ必要があり、単一generic minorから全FULLUへは移さない。

## 4. 保存old QQ matrixによる一回のgeneric-rank preflight

新しいcore pointを取らず、既に凍結・source-boundのexact QQ57×91 matrix at(a,e,d)=(1,2,3)を読む。元qfourの同じ点でのZ=W制限と16 genuine cubicsから36×40を作る。20quad productsはrank19、genuineのheavy5はrank5。二つをexact QQ Schurで除き11×16を得る。original heavy7のrank/非零minorも確認する。

oldpointでrank11なら、保存したrational-function Schur minor circuitがその点で非零なので、QQ(a,e,d)上のrankは11である。上界11はsourceの本数そのものである。これは全点へのpromotionではない。rank不足ならgeneric不足とは解釈せず、別点やsolverへ再試行しない。選んだquad19/heavy5/11minorは追加chart factorsであり、FULLUの既知unitには加えない。

不足43方向やfull239を解決したとは主張しない。新しい実計算対象はこの小cubic mapのglobal rank locus、projected V48 gate、または最終5-target mapのsource-bound coefficient certificateである。
