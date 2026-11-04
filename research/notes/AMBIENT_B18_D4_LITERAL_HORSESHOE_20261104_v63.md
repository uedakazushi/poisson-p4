# Actual同R0 D4のliteral 9・44・112 horseshoe

2026-11-04。独立受理されたactual N coneと、凍結済みnil-delta extractorの36 residueを
結合する有限DAG構築。extractorそのものをこのconstructorが独立再監査したとはしない。
新しいsource solve、jet、point、prime、GBは無い。
T3は別のstructural packetで零と証明されたが、本packetはその証明に必要ではない。
T0,T1,T2と元全Dのdepthはまだ判定しない。

## 入力と同じ係数環

R0=R8[[v]]、R8 regular coordinatesはz1,...,z8。
actual N module generatorはn0=1と、physical rows18,19,28,29,45,46,60から
実even partを引いた七つの関数である。
nil generator Fは旧actual target function `2*aa+v0+z25`。

nil extractorはZ=z1をparameterとするpointed curve77-unit系と、linear-dual77-unit系を
全次数で指定する。各N関係のD4 liftへの値は正符号でc_alpha(Z)*Fであり、
Fはnilcurve上free basisである。
曲線は同R0のregular ideal

```
P=(z2-psi2(Z),...,z8-psi8(Z),v-psiv(Z)), Z=z1
```

を持つ。今回、extractorの全pure-unit scalar DAGを別名でcloneし、parameter Zだけを
ambientのactual z1 nodeへ結合した。curve側N56とambient N56の未知数を同一視しない。
全77+77式、未知数、単位determinant、c36、psi8を同じDAGへ保存した。

## d1とd2

D4 generatorを(F,n0,...,n7)の順に取る。
d1の最初の八列はP_i*F、後の36列は

```
Nrelation_alpha-c_alpha(Z)*F.
```

actual N d2の各entry S_alpha,betaはz,vの一次式であり、nuを含まない。
x=(z2,...,z8,v)として

```
h_i,beta=sum_alpha c_alpha(Z)*coefficient_(x_i)(S_alpha,beta)
```

を有限の定数線形結合で作った。
N関係間のsyzygyをactual D4へ適用すると `(sum c_alpha*S_alpha,beta)*F=0`。
ann(F)=Pなので、このsumのcurve restrictionは零。
一次性によって正確に `sum c*S=sum P_i*h_i,beta`。
従ってd2はnilcurve Koszul28列と、上段h・下段Sの84列、計112列である。

全84 curve obstruction nodeも保存する。
その原点値零だけを全次数の零性へ外挿したものではなく、上のactual identityとfree basisの
証明により零になる。普遍templateではd1*d2の残差がこれらだけであることを全係数照合済み。
新しいHadamard coefficientやsource計算は必要ない。

## 完全性

任意のD4 relationをNへ射影すると、literal Nの全36関係で生成される。
それらを引いた残りはFのみで、ann(F)=Pの八関係で尽くされる。よってd1は全kernelを生成。

次にd1間の任意のsyzygyをNの36成分へ射影する。
Nの全84 syzygyで引くと、残りはnil側八成分だけでPとの積が零。
regular8 sequence Pの全28 Koszulで尽くされる。
従ってactual同R0上で

```
R0^112 --d2--> R0^44 --d1--> R0^9 --> D4 --> 0
```

は正確である。112という数を最小Betti数と主張しない。
extractor ring mapの全inputとnilmodule同定を使うため、型や点tangentだけからの推論ではない。

## 残る比較

old honest9からこの九generatorへのunit transition B_fullと、old Qn relationの44列へのliftは
今回未計算である。旧r/g functionをratioとして零locusへ代入して再除算する操作も無い。
今後それらを完成すれば、この112列を使って低次数torsion moduleを正確に扱える。
元B18全Dのdepthや全nilradical、extensionの分裂は主張しない。

Code `code/ambient_B18_D4_literal_horseshoe_20261104_v63.py`。
Report/folder `AMBIENT_B18_D4_LITERAL_HORSESHOE_20261104_v63`。
`actual-D4-horseshoe-DAG.json`にactual coefficient functions、全unit入力、全44/112列、
generator mapと全84 curve identity nodeを保存する。
