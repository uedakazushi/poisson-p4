# R2L・m二 nil whole source の必要 resonance p=−2/7

2026-11-08、PRIMARY G。全 p の源完全性／Ore 還元は別担当。
本 note はただ一つの必要値 p=−2/7 の **全 fiber** を扱う。
元三十四源 report の generic complete K1 をこの値へ正確に評価し、
新しい M1 profile、minor/jet 列又は標本値を作らない。
p=−2/7 は保存済み generic complete-kernel/minor の分母・分子の零点ではない。

## 1. 元 P0 と source

元全 p family は

```
{a,c}=−a², {b,c}=2ab, {c,u}=(p−3)au,
{b,v}=au, {c,v}=bu+pav,
```

残 brackets 零、actual G2 は `(p−4)a³u²`。
全 K1 十一と **FULL P2 九**を保持する。
P2 を homogeneous kernel 三へ制限しない。
`T=coefficient(P1_av,u²)`、`m=coefficient(P1_cu,v²)` を元係数として定義する。
必要 Euler gauge `E∧(∂cPbc ∂b)` は両 functional を変えない。
T 零及び m 零の全分岐は別の whole-p 証明へ戻し、本計算では捨てたままにしない。

## 2. Tm 非零の全実 normalizer

actual push を

`(a',b',c',u',v')=(A a,B b,C c,U u,V v)`

とし、全 tensor scale A/C を施す。
P0 を正確に保つ条件は `C=B²/A, V=BU/A`。
この群は p を動かさず、元 source functional を

`m'=A³m/(B²U), T'=A²T/(BU)`

へ写す。従って Tm 非零の全 fiber は

```
A=1, B=m/T, C=m²/T², U=T²/m, V=T
```

で **m'=T'=1** へ入る。
全変換は source の unit 上で定義された actual PGL map と tensor normalization、
inverse は非零 diagonal の逆である。
これは二つの任意係数を標本的に一と置いた切片ではない。
全 other K1 coefficients と P2 九も同時に push し、全 normal weight source 内に戻る。
この actual identity は別作者 quartic_review にも係数で確認された。

## 3. 有限な original-support 判定

保存済み三十四源 report の full30 basis／complete generic eleven kernel を使い、
p=−2/7 へ評価した全 source を組み立てる。
m=T=1 による actual normalizer を適用後も、残九 P1 と全九 P2 の十八変数を保持する。
元全三百五十 Jacobi と actual Pfaffian をこの全 source で一回計算する。

使用を許す有限 elimination は、この元二次 equation の QQ-rowspace に現れた
linear equations を正確に解く Gaussian elimination だけである。
毎回少なくとも一つの元変数を消すので回数は十八以下、次数・jet・座標源を増やさない。
各 linear row の元三百五十式への lift と substitution を保存し、
元 zero-locus と substitute 後の zero-locus の全単射を確認する。
constant equation が出れば全 unit branch は空。
actual D の polynomial proper factor が全 residue source 上で得られれば同様に空。

全 filtered group fiber の下界十は三十四源 report で証明済み。
m,T の二 normalization を課した後の fiber は dimension theorem で少なくとも八。
残 source が高々八変数へ exact Gaussian elimination されれば、非線形式を全部解かなくても
`24−8+8=24` が全 projective image の上界となる。
それ以上の変数が残る場合はこの粗上界で閉じたと呼ばない。

本 source の whole-p 認識、p五及び他 operator-wall の global coefficient bridge、
T 零／m 零／p−二は別担当の proof と必ず合わせる。
この一 fiber を処理しただけで全 rank二 plane 又は全研究の完了とはしない。
