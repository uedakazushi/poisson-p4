# m三squarefree四candidate接続のglobal quadratic lift

2026-11-07、主義務G。producer手展開。全係数と全Jacobiは次の有限batchで再構成する。
`PLANE_DOUBLE_LINE_M3_SQUAREFREE_FINITE_CONNECTIONS.md`§§4–5の四c-dependent候補に限定する。
flat meromorphic connectionだけからglobal quadratic liftの存在を推論しないための追加整合条件である。

base(a,b,c)、normal(z,w)、f=a²bとする。
canonical P_fは{a,c}=−a²、{b,c}=2ab。
normal connectionβからτ=β−(dlogf/3)Iを引き、
`iE iV vol=3fτ`のtracefree linear-field inverseを取り、raw mixed bracketを−Vで戻す。
normal-normal raw bracketは未定のglobal quadratic K=k0z²+k1zw+k2w²として保持する。

## C=a²、T=a^g b^(4−g)、g=2/3/4

mixed bracketsの全非零項は

`{c,z}=(4−3g/2)az`、

`{a,w}=(3−g/2)az`、`{b,w}=−(g/2)bz`、

`{c,w}=(g−3)cz−2aw`。

それ以外のa/b/z mixed bracketsは零。
tracefree fieldsはdiagonal entriesでVzz=(3g/2−4)a∂c、Vww=2a∂c、
off-diagonalでVwz=diag(g/2−3,g/2,3−g)であり、表示したβを全global componentsで戻す。

chart b=1、t=a/b、r=c/bでの全affine bracketsは

`{t,r}=−3t²`、`{t,z}=0`、`{t,w}=3tz`、

`{r,z}=(3/2)(4−g)tz`、`{r,w}=(3/2)(g−2)rz`、

`{z,w}=K−(g/2)z²`。

従って全normal-quadratic Pfaffianは

`−3t²[K−(6−g)z²]`。

位数三の必要条件D2零は一意にK=(6−g)z²を与える。
g=2,3,4それぞれのK=4z²,3z²,2z²は正当なglobal quadraticである。
この一意liftとfull350 projective Jacobiを有限batchで検証する。

## C=b²、T=a²b²

mixed全非零項は

`{c,z}=az`、`{b,w}=−3az`、`{c,w}=4aw`。

global βのoff-diagonal inverseはVwz=3a∂b。
chart b=1では

`{t,r}=−3t²`、`{t,z}=0`、`{t,w}=3t²z`、

`{r,z}=3tz`、`{r,w}=6tw+3trz`、`{z,w}=K−3tz²`。

全normal-quadratic Pfaffianは

`−3t²K+18t³z²`。

従ってD2零にはK=6tz²=6(a/b)z²が必要。
これはglobal K∈Sym²<z,w>に入らない。
homogeneous coefficientで見るとKの項はa²b、定数残差はa³z²であり、全三K係数で消せない。
この候補は全first-symbol sourceが不整合として空となる。
有限batchでは全30normal-quadratic coefficientsの3未知数affine systemに、rank三対augmented rank四の証明書を保存する。

この追加obstructionはconnectionのglobal pole-order判定とは別であり、両方を維持する。
