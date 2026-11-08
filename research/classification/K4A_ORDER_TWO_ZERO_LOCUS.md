# K4a の actual24 親：order二零点の全大域 operator

2026-11-08、PRIMARY E。完全equivariant切片と別作者再生により、
K4aはsmooth unique24親の境界、その全非零parameter点は一dense24orbitと分かった。
親の名前はまだ未同定。以下は同じ全非零orbitの order≥2零点を全て一括して調べる別の幾何である。
新anchor・normal・profile・小行列の列へ進まない。

## Projective tensor の intrinsic な二次零点条件

実Piを任意の斉次二次代表 P、Euler場をEとし、`T=E∧P` の十cubic entryを取る。
point q が Pi の order≥2零点であることは、全 ambient first derivatives

```
∂_l T_ijk(q)=0, (i,j,k)⊂{0,...,4}, l=0,...,4
```

の50個の二次式が零になることと同値である。
Eulerの公式によりこの条件自身がT(q)=0を含み、radial方向も失わない。
zero値でのintrinsic firstjetなので、Euler代表変更とPGLに不変。

これら50 quadraticの全15 monomial係数をAとして、
`A*(q_i*q_j)_{i≤j}=0` の50×15QQ行列にする。
kernel零ならpoint無し。kernel一dimなら5×5 symmetric monomial-value行列のrank一を調べ、
rank一ならunique projective q、他ならpoint無し。
kernelが高次元なら全kernelとrank一Veronese条件を保存し、q chartの列へ探索しない。
これは全global schemeの一つのoperatorであり、任意のPoisson点を選び直す実験ではない。

## T36/T38 全成分への必要な閉条件

accepted completed T36 と T38 の全source代表は x0 を全く含まないquadratic coefficientsを持つ。
従ってq=[1:0:0:0:0]で全T firstderivativesが零になる。
全source parameter/PGL飽和について「order≥2零点が存在する」という条件が成り立つ。
point flagはP4でproperなので、この存在条件は全ambient tensor空間の閉条件であり、
T36/T38の**全 component closure**でも成立する。
period denominator零をgeneric sourceから勝手に除外した結論ではない。

よってK4a dense24orbitにorder二零点が無ければ、親はT36/T38のどちらでもない。
T37はこの必要条件を上記source表示から得ていないので、陰性だけでT37を除かない。
陽性なら、unique q を実PGL/Eulerで標準位置へ送る次の actual source 認識が可能になる。

計算は既literal graphの `t=174`（整数化parameter e=1）を用いる。
全非零parameterは既証明の一orbitであり、parameterの選び直しや新rank計算は不要。
全50×15、全QQkernel、非零minorと全quadratic式を保存する有限一回のjobとする。
