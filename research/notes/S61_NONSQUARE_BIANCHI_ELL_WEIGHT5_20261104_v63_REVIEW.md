# full135 weight5 Bianchi／ell除算の独立否定証明

2026-11-04。凍結 `S61_NONSQUARE_BIANCHI_ELL_WEIGHT5_20261104_v63` の独立監査。
新しいkernel方向やrank探索ではなく、指定された20候補空間に対する否定証明である。
producerを呼ばず、受理済みfull175+tG6 operator、weighted chart、Tn/T34 adapter、実ellから
restriction行列を再構築する。core評価は**既受理の同一点(a,e,d)=(1,2,3)だけ**。

## sourceと全変数を保持するframe

元operatorの真のlocal weight5 piecesは20。各cofactorの全135 source indexとweight条件を
保持する。normal8を(x0,x1,x2,x3,x4,Z,W,H)とし、旧normal変数の
`v3=W,v6=Z+x1−3x0` を座標式から独立にJacobian化した。
T34の追加base列1と全34 U-normalizerを残し、43独立変数を使う。
零にするのは検証済みnuisance rows19,23だけで、Tnormal/Tbase1/Uの該当行は恒等零。
Tbase列2,5はそこに対応する二つのunit座標である。

同一点における全weight blockのdeterminantを、別のinteger-cleared Bareiss法で計算し、
非零と確認した。よってgeneric corefield上のこのfull43 linear frameは可逆。
非線形normal graph、source319の削除、任意のnormalizerの零化は使わない。

実ellの全8係数を元anchorから計算し、ell_x1=−16CJ/243≠0、H係数零を確認した。
作者が「cofactorを写してからellを消去」したのとは逆に、本監査はまず全raw-coordinate
linear functionalをell=0へ制限し、それから全135 cofactorを作用させる。
これで得る5805×20行列の非零1166行を、保存された全行とexact QQで照合した。

F5 constant cofactorは全て保持する。weight5 source indexは29個あり、そのうち同一点で
27 constant行が非零。残る零行もsource layoutから削除していない。
F6とsource index134のG6hも保持する。source index134のlocal equation weightは3であり、
20個のweight5 piecesでそのcofactorが零であることを実データから確認する。
一方、**追加operator row175** のtG6h identityは真のlocal weight6なので、この20候補には
入らない。この二つのindexを区別する。135ではなく134以下で検査したという意味ではない。

## 独立minorと結論

指定された20行を独立再構築した行列から取り出し、各行の分母を払って整数Bareiss消去を行った。
全exact divisionを確認し、作者の有理Gaussian消去とは別に

    det = 2987708695128866016 / 32890625 ≠ 0

を得た。分子62bit／分母25bit。各entryは同一点でregularであり、非零specializationを持つ
minorはQQ(a,e,d)上でも非零。従って20候補のcorefield係数による非零結合で、全135 cofactorが
同時にこのellで割れるものはない。ellはlinearでx1係数がunitなので、ell除算とこのrestrictionの
kernel条件は一致する。

結論は**指定された20個の真のweight5 rows**、43local変数について定数かつcoreについて有理な
結合、そして指定の実ellだけを扱う。高weight、一般polynomial cofactors、他の延長、全saturation、
新しいuniform kernelやfull239を否定・証明しない。generic rankの否定から全core点での否定を
主張することもない。新しい点・prime・GB・source再生成・作者producer再実行は行っていない。

再現：`PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/s61_nonsquare_bianchi_ell_weight5_review_20261104_v63.py`
（60秒/2GB上限）。
