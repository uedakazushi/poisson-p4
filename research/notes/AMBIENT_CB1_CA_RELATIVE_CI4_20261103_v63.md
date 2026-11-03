# Q_bのfirst CB branchとCA：相対CI4と実double-S quotient diagram

2026-11-03。pointのclosed transverse CI4と、新しいQQ(b)全150係数計算を組み合わせる。
actual Dの全環への非flat basechangeはここでは行わない。

## 完了したgeneric計算の保存

一回のbounded workerは98.21秒、2.2GB・120秒cap内で完了した。
全58source frame、36次inverse、全150×3 raw quadraticとimplicit補正後の全150×3を保存した。
公開workerは実行済みsourceとbyte同一であり、codeに別保存する。
通常のpublisherはこの完了済み証拠を再利用し、全150行のQQ(b) rank・厳密式・point specializationを検証する。
98秒の物理計算を再実行したとは記録しない。必要時のみ`--recompute`で全workerを再実行できる。

generic obstruction rankは正確に2、基底rowは63,68。
二式の一次multipleが全三次monomialを張るdetは

    1024 b⁴ / [25 (b²−3)⁴]。

全frame・inverse・raw150×3・corrected150×3がb=1で先行point証拠と一致する。
b²−3はこの選択frameの分母であり、幾何的例外の分類を主張するものではない。

## 相対全環：jetからの飛躍を避ける

smooth20 paired-image Sをbaseとし、その二normal idealをm=(r,z)とする。
全S上のactual intersection tangentは22なので、Sをsectionとしてstraightenした全residualはm²に入る。
GL×bはSをdominantに覆うため、generic QQ(b)上のintrinsic quadratic rank2はS上での3次minor恒等零を与える。
b=1の2次単元minorにより、完成base R=C[[s20]]上で全二次行は二つの選択二次式のR-linear spanとなる。

degree3 spanはb=1で単元なので、選択actual式F_i=Q_i+O(m³)についてm³⊂(F1,F2)+m⁴。
有限moduleへのNakayamaからm³⊂(F1,F2)、従って全higher residualも既に二式に入る。
全相対idealは**正確に二つのhomogeneous二次式**であり、R-free rank4、相対HF(1,2,1)である。
closed fiberだけのlength4をgeneric lengthへ外挿していない。
非退化なGorenstein multiplication pairingはb=1で単元であり、complete char0 base上のformal GL2/Hensel操作により

    Z = C[[s20,u,v]]/(u²,v²)

となる。

## 実CB meetingへのquotient

二CB branchの実meetingは既にdouble smooth20、C2=R[rho]/rho²と分かっている。
CB1∩CAの二normal(r,z)における、その実rho tangentの座標はQQ(b)で

    (c0,c1)=(−b³+b,−2b²)。

これは全source58frameを使った実quotientから計算したものであり、抽象的なlength2 quotientを任意に選んでいない。
kernelのlinear generatorは

    eta=c1*r−c0*z=−2b²r+(b³−b)z。

保存されたQ1,Q2に対し、厳密に

    eta² = [9b²(b²−3)²/8] Q1 + [15b²(b²−3)²/4] Q2,
    (b*r−z)² = (b²−3)² Q1/24

である。二つのisotropic coordinateの変換detも、後者のrho tangentへの値も−b²(b²−3)で、b=1で単元。
従って相対完成環と実quotientのdiagramを同時に

    R[u,v]/(u²,v²) → R[v]/v²,  u↦0

へ正規化できる。kernel(u)はR-free rank2でsquare-zero。
これはeta²の一点零だけからの推論ではなく、QQ(b)恒等式とdominant GL×bによる相対identityである。

このdiagramだけではactual Dの全ideal、nil module、Tor、primary decompositionは決まらない。
旧P_aのfull-D modelをラベル交換することは依然許されない。

再現（完了証拠の短いQQ(b)監査）:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CB1_CA_relative_CI4_20261103_v63.py --output <copied-evidence-directory>

evidence directoryにはgeneric-quadratic.json.gzとquotient-map.jsonを先に置く。
全物理workerの明示的再実行は同じcommandへ`--recompute`を付ける。
新しいGBや無制限なhigher-jet runは不要。
