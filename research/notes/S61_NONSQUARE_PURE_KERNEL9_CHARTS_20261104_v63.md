# 残るpure cubic kernelのsource付き七chart

2026-11-04。H7開集合とBgood上に残り得るfull cubic defectを調べるため、既存sourceの最小のcharted detectorを構成する。新しいminor・rank・点・素数・GB・parameter多項式展開は行わない。

canonical rank10の範囲を先に固定する。最初の保存10minorはgeneric lower10だけを与えた。その後のLFRAME_BORDERS_REVIEWは、10Qと5heavyをSchur消去したcanonical11×16でG6行が実際に零であることを全係数確認し、generic exact10を与えた。これは全16genuineの固定pure projectionへのupper boundではない。現在の14×21には10pure Q pivotと一つのmixed H*L genuine pivotがあり、11pure pivotとは数えない。

保存済みcompleted-R11x16の零行4（genuine relative6）を落とすと、実際にmaterializedな10×16行列Rがある。nonzero parameter係数は95376項である。元のsourceはLframe正規化済みG16で、既存ZW25_AUGMENTED_INPUTのpsi16は全く同じ正規化を使う。Q10のpsiは零で、G6とG1のpi・psiは等しい。

既存transcriptの最初の206 ELIMINATE・231 NORMALIZE（733511 bytes、parameter項30417）をpsi、HZ、HWの三列へ同じ順で作用させる。最初のcanonical pure pivot(row14,col25)の直前で止める。G6行のpsiは構造的に零になり、残る10行に追加した10×17=(R,psi)を得る。各row operationはparameter係数のみの線形演算なので、初期恒等式psi=-wdelta*HZ+zdelta*HWも保存される。全26×26 source変換と逆演算recipeも保存する。

全136個の非定数u/g逆因子の出現と、さらに15個の実pivot係数をlocalizerとして保持する。単にu/gを可逆にするだけではpivotを削除できるとは限らない。known FULLU unitとの照合は保守的なliteral testだけで、UNKNOWN分類は真のbad divisorだとの主張ではない。このchartの補集合がBbad13や既知Hdrop曲線で尽くされるとは仮定しない。

残る10 sourceのgenuine relative labelは[0,3,4,5,7,9,12,13,14,15]である。psiが恒等零の3行はG18・G39・G68、relative[0,3,7]。従ってpsiの非零候補は七つである。H7かつ保存H5 chart上では、Cdeltaのprimitive lineを除いた残りのheavy quotientがrank6なので、この七つのpsiは同時には零にならない。各psi_q≠0 chartで、三つのpure sourceと六つのpsi_q e_j-psi_j e_qがkernel9の基底を与える。非pivot9列のminorはpsi_q^6という対角形である。

**元の七heavy係数の補正を省略しない。** piとpsiを殺したsourceでも、そのraw HZ/HWがCdelta lineに残ることがある。既受理のcZ*zdelta+cW*wdelta=eta=192*a*F*B²*C³/125はFULLU全体の単元恒等式である。各kernel vectorの残るHZ/HWを(hZ,hW)とし、(cZ*hZ+cW*hW)/eta倍の元source(G6-G1)を引く。恒等式

    eta*hZ-(cZ*hZ+cW*hW)*zdelta=-cW*psi
    eta*hW-(cZ*hZ+cW*hW)*wdelta= cZ*psi

により七heavy係数がすべて零になる。Cdeltaのpiは零なので、R上の16純係数像は変わらない。保存するのは初期26 source cofactor、元57 sourceのgenuine16 cofactor、補正量、kernel9×16の全rootである。

このchart上ではrank(aug25×32)=15+rank(10×17)=16+rank(kernel9×16)で、既存17×24との関係はrank(17×24)=8+rank(kernel9×16)である。従ってkernel9×16 rank9はfull cubic57の十分条件である。実際のfull cubicの同値条件は、raw H7 kernel9からpure43への写像の単射性であり、ここで計算する16方向へのprojectionのrank不足は反例を与えない。

再現: `.venv/bin/python code/s61_nonsquare_pure_kernel9_charts_20261104_v63.py --output work/pure-kernel9-fresh`。一回60秒・2GB・32MiB、active parameter係数150000項を上限とし、新しいparameter productを展開しない。新しいtyped source packetの独立レビューは別途必要である。

構成は1.290秒、scalar node10586、active parameter係数140820項、全artifact約0.44MBで完了した。三pure行はactual135のG18,G39,G68、残り七候補はG40,G43,G86,G115,G120,G129,G134である。全26 source変換、逆順の逆演算、七chartそれぞれの9source cofactorと元genuine16 cofactorを保存した。Cdeltaの二つのheavy-zero恒等式は独立な六形式変数の整数多項式として全係数確認し、七つのsource basis minorの対角係数も確認した。

## 同値chartと、より広い十分条件chart

元の136 inverse factor出現と15 pivotは、rank同値・source変換の可逆性にはすべて保持する。一方、実際のforward kernel cofactorとprojected imageのancestorを静的に調べると、逆元nodeは28個だけで、27個はNORMALIZEのscale、残りは既知FULLU単元etaだった。ELIMuの逆元や15pivotの逆元はforward式には新たに必要ない。非定数NORMALIZE_gは38出現（保守的literal testで未知26・既知単元12）である。

従ってforward式のsource恒等式とheavy-zero恒等式は、FULLU∩D(product38g)全体で評価できる。uやpivotが零でも、これらの前向き多項式・局所化恒等式自体は消えない。この広い領域のH7上で9個のprojected pure imageが独立なら、9個のactual heavy-kernel sourceが独立となりfull57を与える。ここでは先にsource変換が可逆だとかkernel basisだとは仮定しない。rank同値の元chartと、この広い**十分条件だけ**のchartを区別する。38gの零集合を除去したとの主張はない。
