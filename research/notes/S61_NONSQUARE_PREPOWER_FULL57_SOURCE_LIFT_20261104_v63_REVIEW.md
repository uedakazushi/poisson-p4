# projected prePOWER から actual full57 への source lift：独立監査

2026-11-04。全181 typed node、元63の正規化、全91 physical target、eta補正を監査する。
既存係数は外部参照のままにし、新しいsource値・jet・rank・POWER207展開を計算しない。

元matrixは63×91、六依存行を除いたprimitive matrixは57×91である。保存parserは
normalized row × normalization = original row の方向を持つ。全63倍率が非零QQ定数で
あることを個別確認し、old25 quad sourceのQ0倍数は9/20、scaled Q5倍数は
1600×(−5/4)=−2000となる。Lは保存された同じ五係数をphysical(v0,v1,v2,v4,v7)に
掛ける。restriction上だけ同じになる別のliftへ交換しない。16 genuine倍率muも全て
nativeのC除去指数と一致する既知C冪である。relative6を除いた15 sourceのactual135番号、
full57番号、元63番号を全て照合し、元63へ戻す逆倍率も保持した。

Cdeltaはdelta*(full57row47−full57row42)、delta=−64000C²である。正のsourceはactual135=61、
負は20。二行のgenuine clearing muは等しくdeltaで、pi/psiの差は厳密に零。
これはfull physical source自体が零という意味ではない。physical91の全順序と
pure84/Hlight7を保持し、H係数行は実constant7changeの転置を右から掛けてoldXへ移す。

全physical91→ZW/L31→psi追加32→八unitのtarget削除→GL6→GL3→current22を照合した。
Lのx1係数は−40CというFULLU既知単元。九monomial削除はL²/Lx0のproductであり、
後の最初のgenuine pivot original32=26は **H*L** である。これをpure pivotと数えない。
後続二つのtarget変換はcurrentの四H成分とpsiを変えないので、κおよび残余四行の
currentheavy5零性とprefixのHL零性からoldXの五H成分が零となる。

残るhZ,hWについてc=(cZ hZ+cW hW)/etaとし、sourceをw−c Cdeltaへ補正する。
固定小多項式の全QQ係数を独立に照合して

    cZ*zdelta+cW*wdelta = eta = 192*a*F*B²*C³/125

を確認した。別の六形式変数の恒等式

    eta*hZ−(cZ*hZ+cW*hW)*zdelta = −cW*psi,
    eta*hW−(cZ*hZ+cW*hW)*wdelta =  cZ*psi

も独立確認した。従ってηのみの正則除算でactual raw H7が全て零となる。
HZ/HWの個別係数や未知因子を反転しない。Cdeltaのpi/psiは零なので補正後のprojected
imageは同じだが、元pure84 imageが補正前後同じだとは仮定しない。

この有限source operatorはκ、残余四行、新十一pivot行をactual full57と元63へ明示的に
戻す。κはFULLU全体、残余四行と十一pivotはFULLU∩D(Delta407*Saux*T)で定義される。
既受理のproduct41とrestricted product19のprojection theorem、Cdeltaのprimitive H単元を
合わせると、指定chartで残余4×11のrank4が得られれば元full cubic rank57が従う。
これでsource liftのgapは閉じる。rank4自体、全FULLU57、全quartic239は未計算である。
全八Q9 tailの効果とchart補集合を引き続き保持する。

再現: `python code/s61_nonsquare_prePOWER_full57_source_lift_review_20261104_v63.py`。
