# full135 Jacobi–Bianchi source operator の独立監査

2026-11-04。元の134式と追加G6を保った135式、nuisance零のQQ[a,e,d]上で監査する。
producerや旧full universal verifierは再実行せず、immutable v02のencoded coefficientsから
別の疎有理多項式演算で現在の全core sourceを再構築する。
C122、四次Macaulay行列、rank/pivot、点・体・GBの計算は行わない。

## 全sourceへの束縛

61変数の元150係数を、固定J、61個の自由raw coefficient、immutable eliminated-raw式から直接再構築する。
betaを(a,e,d)coreへ代入した後、qについてaffine-linearであることを検査する。
affine chartのbivectorを

    P_ij−y_i P_0j+y_j P_0i

という元の座標式から作り、qに依存しない全項をt倍してhPとする。
これを保存されたhPの全spatial/q/t/core係数と比較する。
従ってq-leading二次部分だけでなく、t*J、t*L、t*Q_beta²、t*C_beta³も保持される。

別に元134個のencoded関係とuniversal G6を同じcoreへ代入し、q/tについて二次に斉次化する。
保存された135個の全係数との一致を確認する。G6列を消したり、tを単元として扱ったりしない。

旧80/140 raw Jacobi mapsの**全61変数universal**な正しさは、既受理の
UNIFORM_CUBIC独立監査とその元コード・報告hashを入力とする。
このgeneric61検査は今回重複しない。
代わりに今回の全core hPから、Jacobi trivectorを直接微分・乗算して新たに計算する。
旧154式→134式の20 redundancyも、constant部分とbeta-linear部分を別に再構築する。
保存された80個の空間三次、140個の空間四次の係数を、再構築したFh135へ戻した値が
直接Jacobiの全係数に一致することを確認する。
空間低次の60係数と四次を超える係数は零である。
これはfull raw mapの現在の全core特殊化について、単なるq-leading照合より強い独立計算である。

## 全Schouten行と135恒等式

raw Jacobiをfree generator symbolで記したTと、上の全hPとのSchouten bracketを
alternating/unshuffle添字から別に作る。
保存された空間次数4/5/6の35/56/84行、全175行の全係数と一致する。
空間次数2/3のoperator係数が134へのredundancy後に零であることも確認する。

Tへ実際のFhを入れたものは、前節で直接Jac(hP)と確認した。
従ってgraded Jacobi恒等式

    [hP,[hP,hP]]=0

により、今回再構築した175行はすべて元135式の**文字どおりの**syzygyとなる。
producerと同じ175個の三次多項式再展開を繰り返す代わりに、
独立に確認したJacobi correspondenceとSchouten構成からこの結論を証明する。
報告では、この証明と「175三次多項式を再展開した」という別の計算を区別する。

空間次数4のq-partは、旧canonical35ベクトルと全係数で照合する。
追加の非Bianchi関係は、旧125項cofactorから独立に

    t*G6h−ΣK_i(q)*F4h_i

を組み立て、この一つの三次多項式は全core係数を実際に展開して零と検査する。
これにより134→135の扱いとtの係数・符号も束縛する。

## 実際の重み部分と小block

全175行と追加G6行を、q_fixed=1へのdehomogenizationの後、元のlocal coordinate/equation weightsで分ける。
全weight piecesのsource係数とそのindexを再構築して保存されたものと比較する。
実際に非零なspatial operator row数は

    degree4:35, degree5:56, degree6:0。

84個のdegree6枠は全て文字どおり零である。
weight4 piecesはdegree4から5個、degree5から10個、計15個に限る。
35個のleading canonical関係をfull operatorと同一視したり、個数から80個の独立関係を主張したりしない。

weight4のgenuine quartic33列へは、constant multiplierの係数だけを抜く。
そこで現れる非constant weight-zero multiplierはlocal座標19,23だけである。
旧正確なweighted chartで、T34_normalとU34のこれら二行が**全要素零**であることを確認する。
従って固定normal graph上で二座標は全次数で正確に零であり、blockから該当132項を除いてよい。
全operatorと全weight pieces自体からこれらを削除してはいない。
この制限は対象normal graphでのものであり、元の45変数chart全体で二変数が零という主張でもない。

得られるconstant blockは15×33、非零行11、非零core項88であり、全係数を比較した。
まだT34/Tnの全normal transport、implicit graphへの全係数代入、319行への35削除のliftは行っていない。
このblockのsourceから後の四次kernel方向を採用するには、それらのexact liftが別途必要である。
今回の監査はblockのpivot/rank、80方向、全239列全射を証明しない。

## frozen sourceの保全

authorのsource-verified-checkpointと最終operatorの共通する全source/135式/175行/追加G6行が
完全一致することも確認する。後段のweight-zero block処理がfull sourceを変えていないことを区別して記録する。

専用コード：`code/s61_nonsquare_full_bianchi_operator_review_20261104_v63.py`。
受領書：`reports/S61_NONSQUARE_FULL_BIANCHI_OPERATOR_20261104_v63_REVIEW/report.json`。
上限60秒・2GB。途中で上限に達した場合は進行位置と未完了を保存し、自動的な大きいretryは行わない。
