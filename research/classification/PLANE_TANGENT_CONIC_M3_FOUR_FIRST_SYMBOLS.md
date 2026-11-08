# 接線conic+line・m三：全四first symbolsのliteral Poisson代表

2026-11-08、G、producer plane_collision_review。純手証明、独立review前。
`PLANE_CONIC_LINE_M3_ROOT_COVERS_AND_EXCLUSION.md`の全root-cover還元を仮定し、その四profilesのglobal P0を具体化する。
full P1/P2の計算はまだ行わず、全G上界・E帰属を主張しない。

## 全共通座標と十brackets

base=(a,b,c)、normals=(u,v)、L=c、Q=ac−b²、f=cQ。
全profileのbase bracketsは

`{a,b}=2ac−b², {a,c}=2bc, {b,c}=c²`。

V=2b∂a+c∂b、H=−a∂a+c∂cとする。
mixed columnsを`{x,u}=uX`, `{x,v}=uZ+vY`と書く。

| profile | X | Y | Z | {u,v} |
|---|---|---|---|---|
| S | V | V | 0 | 0 |
| N | V | V | −3V | 0 |
| D | −V/2 | V | 0 | 0 |
| M | V | −V/2 | 3H/2 | 3u²/2 |

ordinary connectionとraw mixed fieldの係数は
`βii=(dlogf−α_(diagonal field))/3`、off-diagonalは`−α_Z/3`。
したがってS/N/Dの表示はそれぞれscalar、scalar+nonzero nilpotent、half-logQ/logLの四profile表と一致する。
Mでは

`i_E i_(3H/2)vol=3Q db−(3b/2)dQ`

だからoff-diagonalは`(b/(2c))dlogQ−db/c`である。
D2のu²係数は`3f K20−Z(f)−det(E,X,Z)`。
Xのcharacter p=dlogcより、この零条件はK20=p(Z)=3/2を与える。
他のnormal coefficientsは零で、全normal quadratic consistencyを戻している。

## P0全体のPoisson性をchosen numerical testなしで確認する

S/N/Dでは全linear fieldsがV又はVと可換なconstant normal fieldsの組で、V(f)=0、trace V=0。
従ってP_f及びmixed partsの全Schouten bracketが零、raw strict Poissonである。
actual Pfaffianもnormal quadratic零なのでrank≤2である。

Mについてc=1 chartを取る。uはfull first-symbol PoissonのCasimirで、(a,b,v)の三bracketsは

```
{a,b}=3(a−b²),
{a,v}=−3(au+bv),
{b,v}=−(3/2)(bu+v).
```

三vector J=({b,v},{v,a},{a,b})は
`J=(-(3/2)(bu+v),3(au+bv),3(a−b²))`。
curl J=(-9b,−9/2,(9/2)u)、従ってJ·curl J=0が全a,b,u,vのliteral恒等式。
ゆえにMもprojective Poissonであり、u Casimirからrank≤2。
dense chart上のprojective Jacobi零は元homogeneous全350式へ延びる。
これはraw strict Jacobiとprojective Jacobiを混同した証明ではない。

## 全regular leading spaces

全G3は以下の全parallel spacesに入る。

* S：`Q Sym³(u,v)`、dimension四。
* N：`C*Q u³`（zeroをsupersourceへ含めるならC line）、dimension一。
* D：`span(c²u³,Q uv²)`、dimension二。
* M：`C*(bu+cv)(a u²+2buv+cv²)`、dimension一。

NではQ-meridianのmonodromyはnontrivial unipotentである。
rational parallel tensorはそのinvariantsに入り、Sym³上のnilpotent kernelはu³だけ。
scalar equationを戻すとcoefficientはQのconstant倍となる。
Dは各monomial u^k v^(3−k)のcoefficient equationから、k=3でc²、k=1でQ、他kで非rational square rootを要求する。
Mはframe z=bu+cv,w=uで(az³+bQzw²)/cとなり、polynomial性がa=bを強制する。
全leading amplitudesは固定せず、Dの各zero-coefficient境界も保持する。

## 有限whole-source計算を実行するときの正しいdimension budget

tangent pairs (L,Q)の固定plane内moduliはdimension六、anchor amplitude一を戻して七。
Sはnormal GL2-conjugacy dimension零、Nのnonzero nilpotent coneは二、Dのdistinct diagonal pairは二、Mは高々三。
従って全first-symbol sourceの安全なupper Bは順に

`S:7, N:9, D:9, M:10`。

Mのlinear coefficient identificationはfixedQのisometry frameであり、constant normal GL2が全てを覆う。
Rを固定したまま別のbase Aut自由度をさらに加える必要は無い。
全higher P1三十、P2九を戻し、P1にweight一Jacobiと全G3-space membershipを課す。
任意P1上のP2 fiberは空又はfixed weight二operator kernelのtranslate。
その両kernel dimension和をdとすれば、全plane六とtensor projectivization一を戻す像は

`B+d+5` 以下。

よって必要な最弱targetsはd≤12,10,10,9の四つ。
これは単一P0 witnessのdimensionを全familyへ推測した値でなく、全root-cover正常化のactual inverseと全first-symbol source countから来る。
次に計算する場合、四全profilesを一回のbounded batchで扱い、全parallel spaceとactualD3 imageを保持する。
profileが未閉なら全kernelと幾何障害を保存し、同じminorやjetを追加しない。
共有runnerは所有者調整前には変更していない。
