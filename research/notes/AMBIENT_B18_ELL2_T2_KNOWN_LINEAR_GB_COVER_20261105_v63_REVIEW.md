# 全432零reduction・global linear GB・coarse cover停止の独立検証

全270 NF nodeを独立に代入し、全48真の一次式×九変数の432列が零になることを確認する。rankは零であり、既知線形coverの二次層は六次元である。新しいrank探索やproducerの再実行は行わない。

元多項式moduleに48 normalized linear式と六ell²e_targetを加える。全171同一target線形S-pair、および48個のell²とのcritical pairが尽くされ、全て零reductionを持つ。別target間にmodule S-pairはない。既証の乗法的well-orderとBuchberger判定によりglobal GBとなる。

standard monomialはdegree0の六generatorと、正次数ではbottom三targetのr^d/r^(d-1)ellだけなのでHilbert seriesは6/(1-t)、先頭四層は(6,6,6,6)で長さ24になる。I²⊂m4だけを追加しても次数三までは変えられない。実T2の既証上界20を改善できないため、このcoarse cover routeを止める。実T2には未知の真の二次初項があり得るため、実HF2/3の下界・等号・全長は主張しない。

完了evidence後のreport path shadow失敗は元code・hash付きで保存され、metadataだけ修復されたことをsourceとして束縛する。数学evidenceを独立に全検査し、producer arithmeticは再実行しない。
