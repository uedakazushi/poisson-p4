# 同じraw四次jetから真の初期idealを作るblock線形消去

2026年11月1日。これは係数体K上の一般的な線形代数手順であり、
特定のS61層の四次rankや局所環を単独で主張する証明書ではない。
再利用helperは `code/s61_initial_quartic_from_jets_20261101.py` である。

## 入力と出力

九変数に限らずn個の形式変数vについて、同じ式系

    E_i=A2_i+A3_i+A4_i+O(m⁵)

の三つの斉次係数行列を入力する。
単項式順は各次数の `combinations_with_replacement(range(n),degree)` とする。
低次式を変更した後のA2と、変更前のA3/A4を混ぜてはならない。

関数 `initial_quartic(A2,A3,A4,n,domain,compute_rank=False)` は
Kを体とし、二次・三次の独立先頭行と追跡行操作、四次の二つの行列
`quartic_lower`、`quartic_extra` を返す。
この二行列の行空間の和が、真の初期idealの次数4部分に等しい。
`compute_rank=True` はそのrankも計算する。
有理関数体を使った場合、特殊パラメータでの一様性は別のminor・分母監査が必要である。

## 手順と正当化

まずA2の独立行を選び、その元の式をFとする。
A2の全行をその基底で表す係数行列C2を計算し、実際の式を

    H=E−C2 F

へ替える。Hは次数3以上であり、その三次・四次尾は

    H3=A3−C2 F3,  H4=A4−C2 F4

である。以後も同じ行操作を全ての尾に適用する。

次数3の先頭行列L3を、全ての `v_j F2` と全H3を縦に並べて作る。
それと同じ順で四次尾T4を、全 `v_j F3` と全H4から作る。
L3の独立行B3を選び、全行をC3 B3と表す。
すると

    R4=T4−C3*(B3に対応するT4行)

は実際の式組合せの四次先頭項である。
これが `quartic_extra` となる。
ここには、二次式の一次syzygyが作る四次先頭項も含まれる。

残りの四次先頭項はB3の全変数倍で張られる。
これが `quartic_lower` である。
実際、次数4までの式組合せに必要な係数は次数2以下に限られる。
定数係数の二次打消しはH、一次係数の三次打消しはL3の全核で処理済みである。
二次式の二次倍とHの一次倍の先頭項は、全てB3の変数倍に入る。
したがってこの二行列が次数4部分を過不足なく張る。

追加のnonlinearな座標変更を仮定せず、入力されたraw jetの式変更だけを追跡する。
全初期idealの次数5以降や完成局所環はこの手順からは結論しない。

## observableの四次

例えばtの二次jetが零で三次jett3がL3の行空間に入る場合、
t3=Ct B3を解いた後で

    t4_corrected=t4−Ct*(B3に対応するT4行)

を使う。raw t4だけを四次先頭idealで割ってはならない。
helperはこの計算に必要なL3、T4、独立行番号を保存して返す。
observableを元の方程式idealに加えることもしていない。

## 小さい厳密検算

F1=x²+y³、F2=xyを使うと、入力のraw A4は零であるにもかかわらず

    y F1−x F2=y⁴

が四次初期関係になる。
helperで二次rank2、三次rank3、四次の下位倍rank4、真四次rank5を確認した。
これは二次一次syzygyを落とした実装を区別する検算である。

```
import sympy as s
from s61_initial_quartic_from_jets_20261101 import initial_quartic
r=initial_quartic(
    s.Matrix([[1,0,0],[0,1,0]]),
    s.Matrix([[0,0,0,1],[0,0,0,0]]),
    s.zeros(2,5),2,s.QQ,True)
assert (r['quartic_lower_rank'],r['quartic_initial_rank']) == (4,5)
```
