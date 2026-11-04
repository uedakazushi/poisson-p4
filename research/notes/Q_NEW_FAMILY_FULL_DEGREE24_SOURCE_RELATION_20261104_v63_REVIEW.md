# 全Inv24のsource制限は単射：独立67係数minor

2026-11-04。受理済みの普遍67基底と、実際の10個の低次数source不変式mapを
入力とする。source点・beta点は評価しない。普遍70候補のsemantic定義から
10生成元の指数を**新たに**作り、62,68,69を除く67列と照合する。

全列をB12H24へ払う際の整数constant denominatorが32003単元であることを
確認する。独立nmod engineは変数順(B,W,Y,a)、LEXを使い、著者の整数・
(a,Y,W,B)・DEGLEX workerを再実行しない。保存された67 coefficient行の
全4489成分を再構築し、別のPython eliminationでminor9857 mod32003を得る。

従って整数係数の非零minorがあり、QQ上でも、さらにC上でもrank67。
入力67本は**全Inv24空間の基底**なので、以前の41部分空間だけの除外とは
異なり、全degree24不変式のsource制限が単射となる。

実際のNEW beta像は既受理の既約15次元斉次hypersurfaceである。Sp4は連結
半単純なので、その斉次prime equationは不変式に取れる。受理済みの次数
下界24と今回の全Inv24除外により次数は24より大きい。中心−IはW16上で
beta→−betaとして作用するので不変式の次数は偶数。従って**定義式次数は
少なくとも26**。source quotientと実際のquotientの有限指数1はこの推論に
不要で、必要な支配性は既受理の有限source gaugeとrank5による。

degree26以上の計算や方程式の構成、未知Poisson成分の除外は行わない。
一つの非零係数minorから得る範囲に限定する。

```sh
. /workspace/cloud-setup/activate.sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_full24_source_review_20261104_v63.py --output work/full24-source-independent --work work/full24-source-independent-native
```

60秒・2GB、一回の独立係数worker。自動retry・別prime・別次数への拡張は無い。
