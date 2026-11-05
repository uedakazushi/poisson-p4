# honest ell二乗比較の長さ範囲とH0最小長の判定条件

これは既証DVR構造からの純粋な帰結で、T2の実係数・長さを評価する計算ではない。
D=k[[ell]]のellは同じphysical parameterであり、rAとは区別する。

元の自由比較複体のH1=K,H0=Cについて、DはDVRなのでell二乗への導来base changeは
正確列 `0 -> K/ell²K -> T2 -> C[ell²] -> 0` を与える。
CのDVR rankは四、実C/ell二乗Cの長さは十四なので、torsion部分の構造定理により
`length C[ell²]=14-2*4=6`。

KはDVR rank四・torsion最小生成数三である。三つの正整数指数をe1,e2,e3と書くと

`length T2 = 8 + sum min(2,e_i) + 6 = 17 + #{i:e_i>=2}`。

従って実長は17〜20。これは全形式表示を数値評価した結果ではなく、判定対象の範囲である。
T2長が17なら、三つの指数は全て一であり、逆も成立する。

正規化E4全体でellが正則なので、D0のell-torsionは核KのD-torsionと等しい。
KはD-finiteであり、そのtorsionは有限k長を持つので、これは正確にH0_m(D0)である。
従ってH0長はsum e_i、少なくとも三。以下は同値となる。

- T2の長さ17。
- 三つのDVR torsion指数が全て一。
- ellがH0_m(D0)を零化する。
- H0_m(D0)の長さ三。

ellによる零化から全maximal idealによる零化を推論しない。H0のsocle次元、元nilの全生成元、
nil三乗の零性や正確なnil指数は決まらない。実T2長もこのノートでは未決のままである。
