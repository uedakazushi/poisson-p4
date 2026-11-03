# 実際の B18 Noether 閉 fiber：一変数 source 四次による偶二次 graph 比較

2026-11-03。これは実際の共通28個の標的座標を P0 の値に固定し、残る
`z=q61−q61(P0)` と10個の奇方向を比較する計算である。全29変数の偶二次
係数は展開せず、一変数の source p 四次だけを使う。

## 1. source と精度

CA、CB それぞれで、受理済みの19 GL 方向と固定9方向を28個の source
変数に取る。kernel p の有限行列 `exp(pK)` とその固定 source の補償を
保ち、28個の実際の標的座標を固定28次逆行列で p、p²、p³、p⁴の順に解く。
CB の固定 source は `g=−ratio*f` を全て掛けたものを使い、ratio と binary
係数の混合項を省略しない。各次数の raw source、両 inverse 恒等式、全150
出力と25発散・350生 Jacobi の四次までの零を検査する。

全出力の p 一次係数は零である。fold は
`z(p)=w p²+c p³+d p⁴+O(p⁵)` で、w は CA で−4、CB で3/8。
一般の出力を `q=q0+q2 p²+q3 p³+q4 p⁴+O(p⁵)` とすれば、その trace-even
部分の z² 係数は正確に

`[q4−q2*d/w+(3/2)*(q2*c²/w²−q3*c/w)]/w²`。

実際、二次 cover の involution は
`τ(p)=−p−(c/w)p²−(c/w)²p³+O(p⁴)`。全出力に p 一次項がないため τ の
p⁴係数は不要で、source p⁵も不要となる。コードはこの involution を別に
代入し、上の trace-even 公式と trace-odd 部分を全150係数で照合する。
even の z 一次式と odd の z 一次式は既受理の実際の E、b 一次式に一致する。

## 2. 比較するのは正規化した graph

同じ10個の物理的奇生成元を両族で選ぶ。共通の実際の39接空間から、残る
各物理座標の odd 定数係数 h0 を取る。比較する係数は raw a_j,2 自身ではなく

`α_i,j,2=a_i,j,2−Σ h0_jk*a_i,min_k,2`。

共通線形奇変数 `Y=q_min−E_min R` を用いると、未知の h1 は zY にしか
寄与せず、pure z² の graph 差はこの α の CA−CB 差に等しい。最小39個の
座標では α は自動的に零であることを全て検査する。残る全111物理座標を
保存し、その中の独立な87個へ標準代表の126座標を拡張して完全性を確認する。

## 3. 12/13 判定への応用条件

独立受理済みの finite13 proof は、実際の閉 fiber で二次初項 YY と zY、
三次初項 z³を与える。従って候補 Hilbert function は `(1,11)` または
`(1,11,1)` だけである。

各個別 image の全 kernel は odd ideal の syzygies と pair products で
生成される。その二次初項には YY、zY だけがあり、pure z² はない。
共通最小39座標の外の正則 graph を両族で引いたものが残る全関係であり、
pure z² の有無は上記87個の正規化した偶 graph 差で尽くされる。

一つでも非零なら z² が消え、閉 fiber は長さ12。全て零なら二次 Hilbert
値1が残り、長さ13である。この応用には source 四次の独立検証と、全87個
の graph による網羅性の確認を受理ゲートとして残す。report の判定値は
この二つの検証を経て研究台帳へ反映する。

元の交叉 D の全環や正確な深さを決める計算ではない。p⁵、p⁶、全29変数の
source ジェット、Gröbner 基底は実行しない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_Noether_fiber_even_quadratic_20261103_v63.py --output work/B18-Noether-even-quadratic
```

60秒・2.2GB上限。全 source 段階と比較に使う raw 配列を保存する。
