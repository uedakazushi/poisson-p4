# DB=0の補助rank5：実係数による独立監査

2026-11-03。author packet `S61_NONSQUARE_FIVE_RANK_DB_ZERO_20261103_v63`を受理した。
結論はa*P*J≠0のDB=0全体で、補助五二次式(Q0,Q1,Q2prime,Q3,D)のrankが正確に5となること。
F,B,Cの非零条件は必要とせず、DB≠0やpair gateの結論を含まない。

## 実sourceを再構成した計算

archiveと全input SHAを照合し、804個の実quadric係数から四式を再構成した。
Q2primeとDの定義はliteral QQ polynomial identityとして再確認し、Wをinvertしていない。
W=0のQ1=A*B/6075、Q2prime=40*x0*B/6075、3Sc−7Vc=20Pも再確認した。
Aのx1,x3係数は−3P*Sc、−6P*Vcなので、P≠0ならAとx0は独立、Bも非零である。

split inverse frameのdetは1/(2a)。
4Q0=uR−5vtを直接確認したので、非退化四変数部分のrank4が実可逆frameで成り立つ。
authorのHessian conjugationとは別に、a²Q0,a²Q3,a²Dを実frameへ代入し、各二次monomialの係数を抜き出して
全75 matrix entryを再構成した。保存された三行列と一致した。

## DB=0上の二つのcase

DB=20P²+3Pe+20P+180a、L=20P+3e+20、K=−130P+168e−5とする。
DB=0ならa=−PL/180なので、a,Pが単元である限りLも既知単元。
Bのr係数は−20a*DBで零となり、Q0はrを含まない。
独立に再構成したH3,Hdについて

    H3_rr=0,
    Hd_rr=−P^5*L^4*J²*K/35872267500

をQQ[P,e]上で確認した。
さらに既存のglobal proportional exclusionの五つのcross-minorとliteral QQ Bezout identityを、
今回の実H3とBから再び確認した。よって(b,h3)のrankはDB=0全体で2である。

K≠0ではHd_rrが非零なので、Q0,Q3,Dの制限の定数係数関係のr²項がDの係数を零にする。
次いでr-linear項とrank(b,h3)=2がQ3の係数を零にする。
Q0はrank4 split formの一つのhyperplaneへの制限なのでrankは少なくとも2であり、零ではない。
従ってQ0の係数も零となる。

K=0は捨てず、e=(130P+5)/168として別に扱った。
このとき

    L=125*(10P+9)/56, J=(10P−117)/28.

実(b,h3,hd)の全3×4行列を再構成し、四つの3-minorを独自のFraction polynomial permutation和で計算した。
authorのdet/gcdex routineは使用していない。
各minorと保存primitive polynomialの間の因子は、記録されたP,L,Jの非負冪と非零QQ定数だけであることを確認した。
保存された四つのQQ cofactorを独立のsparse polynomial積で掛け、和が正確に1になることを確認した。
従ってrank(b,h3,hd)=3で、Q3,Dのr-linear classはBを法として独立。
このcaseでもQ0,Q3,Dの制限は独立である。

Kをglobal unitとみなしていない。
K≠0とK=0の二caseが全DB=0を覆う。

## 五式への復元とscope

五式の任意の定数係数関係をB=W=0へ制限すると、上の結果でQ0,Q3,Dの係数が零になる。
残る関係はW=0でB*(c1*A+40*c2*x0)=0。
B≠0とA,x0の独立性からc1=c2=0となり、五式は全DB=0,a*P*J≠0で独立である。

これは補助rank gateのDB=0部分の解決。
DB≠0、pair gate、元raw24/product41、quartic、ambient removalへ拡張しない。
新しいGB・素数・数値標本は使わず、凍結sourceを変更していない。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/s61_nonsquare_five_rank_DB_zero_review_20261103_v63.py --output <new-output-directory>

全再構成行列とcase証明書:
`reports/S61_NONSQUARE_FIVE_RANK_DB_ZERO_20261103_v63_REVIEW/`。
実行約5秒、制限120秒・2.2GB。
