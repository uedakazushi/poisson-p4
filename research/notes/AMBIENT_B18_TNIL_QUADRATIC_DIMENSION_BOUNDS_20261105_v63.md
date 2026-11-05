# 実Told→Tnil単射による二次層の範囲：2≤HF2≤5

2026-11-05。追加のsource係数や6×661行列の値を計算しない、既証の実module間の系である。

## 同じ係数環上の入力

`S=QQ[[zA0,...,zA6,rA]]=R0/(ell)`と置く。`rA`はDVRのparameter `ell`とは異なる。
accepted actual natural comparisonはS線形の完全列

`0→Told→Tnil→Q→0`

を与え、`length(Q)=5, mu(Q)=3`である。既証のuniserial tailと実三次層から
`length(Told)=10+s`, `HF(Told)=(5,3,2,1,...,1)`, `1≤s≤7`、かつ
`m^(s+3) Told=0`である。一方、`length(Tnil)=15+s`, `HF0=HF1=5`である。
単射のorigin rankは2であり、Told=mTnilという誤った同定は用いない。

## 上界

`h=HF2(Tnil)`とする。S線形単射から`i(m³Told)⊂m³Tnil`。
左右の長さはそれぞれ`s`と`s+5−h`なので、`h≤5`を得る。
一般に`length(Tnil/m^k Tnil)≤length(Told/m^k Told)+5`である。
もし`h=5`ならこの包含は長さが等しいため等号となり、三次以降の各graded layerはToldと一致する。

## 下界

長さ5、最小生成元数3のQにはoriginより上に高々二つの非零layerしか無い。
したがって`m³Q=0`であり、`m³Tnil⊂i(Told)`。
これとToldの上のLoewy boundから`m^(s+6)Tnil=0`を得る。

`h=0`ならTnilの長さは10となり矛盾する。`h=1`と仮定する。
`gr_m(Tnil)`は標準次数多項式環のfree module5の斉次商なので、module term orderでmonomial initial submoduleに移せる。
二次standard monomialが一つだけである。任意の三次以上のstandard monomialの全二次因子もstandardでなければならない。
二次monomialが同じtargetのpure squareなら、その変数のpure powerだけが残り得る。
二次monomialが二つの異なる変数を含む場合、三次以上に延びるmonomialは無い。
よって三次以降も各Hilbert値は高々1である。

非零であり得る次数は0から`s+5`までなので、長さは高々
`5+5+(s+4)=s+14`。実際の長さ`s+15`に矛盾する。従って`h≥2`。
コードは8変数×5 targetsの全180個の二次monomialについて三次への可能な延長を確認し、一般次数の議論は上の因子条件で証明する。

## 実source operatorへの帰結と限界

complete814一次零sourceを保持した既存6×661 operatorについて
`HF2(Tnil)=6−rank(operator)`である。従って実階数は1〜4であり、全ての値が零となることは無い。
具体的な非零entry、正確なHF2・総長さ・nil generators・nil cubeは未決のままである。
これはsource scalar値を算出したという主張ではない。
