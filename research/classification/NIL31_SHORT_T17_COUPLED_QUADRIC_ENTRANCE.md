# Short nil31 r=2/3：T17の二factor入口と一つの正則coupled衝突の限界

2026-11-08、PRIMARY E。純手展開、新しいoperator、point、rank、minor、jet又はCASは無い。
全targetは `MODULAR_NIL31_SHORT_SUPPORT.md` §7のtd非零一軌道と、そのd零閉境界である。
以下は全T17非所属又はtarget帰属ではなく、具体的な一実族の正確な診断。

## 1. Whole T17の正しい親と二つのfactor assignment

既独立受理の `TWO_QUADRICS_TAN_SIMPLE_SNC_SEPARATE_REVIEW.md` §2,4から、全親は

```
F=pAB+uH², G=qCD+vF, b=p-a,
P=F ∂A wedge ∂B+G ∂C wedge ∂D
 +(A∂A-B∂B) wedge(aC∂C+bD∂D).
```

全五parameter p,q,a,u,vを保つ真正completed T17 sourceで、G中のvuを独立係数へ緩めない。
actual DはHFG、Fはrank三以下、全Jacobiはb=p-aにより恒等零。
named sourceのp零もpolynomial閉性で含める。

target Dは `h²*Z*Q`、`Q=hY-ZW` はC上既約rank四。
proper factor incidenceで親のH,F,Gを標識して延長すると、outer Hのspecialはh又はZ。
Fのrank三以下は閉なのでFはQへ移れない。
従って全割当は次の二つだけである。

```
outer H -> h: F -> hZ, G -> Q;
outer H -> Z: F -> h², G -> Q.
```

factorに任意の非零scalarを許す。これは正しいfinite幾何入口だが、tensorの実族は別途必要。
旧SS12の分離companion flagsの不整合から、この二入口全体又はcoupled flagsの不存在は従わない。

## 2. 最初の割当を実現する一つのfull tensor family

親からtarget座標へのactual frameを

```
H=h, A=h+e²X, B=Z, C=h-eW, D=Z+eY
```

と置く。detは `-e^4`。
逆は `h=H, X=(A-H)/e², Z=B, W=(H-C)/e, Y=(D-B)/e`。
次の真正period/edge bindingを使い、tensorをe²倍する。

```
p=1, q=1/e, u=e²U, v=-1/e, a=e*k, b=1-e*k.
```

このとき

```
F=hZ+e²(XZ+Uh²),
G=Q-e(WY+XZ+Uh²).
```

両quadricの関係は全eで元のG=qCD+vFそのもの。
全括弧の直接Leibniz展開は次のpolynomialとなり、負e冪は一つも無い。

```
PXY=hY-khZ-e(Uh²+khY)+e²(XY-kXZ)-e³kXY,
PXZ=F,
PXW=-k(h+e²X)(h-eW),
PYZ=e(1-ek)Z(Z+eY),
PYW=Q-e(WY+XZ+Uh²+khZ)+e²kZW,
PZW=e²kZ(h-eW).
```

generic e非零ではactual T17 sourceのPGL/scalar像であり、polynomial性から全eでPoisson。
actual Euler-Pfaffianは変換則e^4/det(frame)から正確に `-hFG`。
中心Dはtargetと同じ `-h²ZQ` である。
これは単なるdivisor候補でなく、全tensorが正則な実曲線を一つ構成した結果。

ただし中心の六括弧は

`PXY=hY-khZ,PXZ=hZ,PXW=-kh²,PYW=Q,PYZ=PZW=0`。

h一chartのcolumn divergenceは
`Nlim=-2∂X+Z∂Y+∂W`、斉次化すると `-2h∂X+Z∂Y+h∂W`。
これは平方零、rank二。
shortnil targetの唯一の二重hから得るcofactorは
`Nshort=Y∂X+Z∂Y+h∂W`、Nshort²(X)=Z非零でJordan型(3,2)。
よって中心はtargetとactual PGL同値でなく、この実族は今回のwhole22軌道への正の橋ではない。
actualD一致だけで帰属を推測しない。

## 3. 同じframeのperiod又は次数を増やしても届かない

この固定frameでp,q,a,u,vを任意meromorphic functions、tensor scaleを任意s(e)とする。
親column divergenceは全五源で
`(-2p A,2p B,-q C,q D)`。
実inverseから全modular場は

```
NX=-2sp*h/e²-2sp*X,
NY=s(q-2p)*Z/e+sq*Y,
NZ=2sp*Z,
NW=sq*h/e-sq*W, Nh=0.
```

regular endpointがnilpotentなら、diagonal entriesからsp,sqは零へ行く。
残るlimitはXをhへ、YをZへ、Wをhへ送り、h,Zを零へ送るので必ずNlim²=0。
したがって全meromorphic period、任意ramification/scaleを許しても、このframeはshortnil targetに届かない。
a/u/vのLaurent次数、tail又はprofileを追加する余地では解消しない。

ここで停止する。未処理なのは、上の二factor assignmentを保った**別のcoupled moving flagsを含むactual T17曲線**、
又は別のnon-SNC親からの実族である。
任意T17閉包非所属、新component又はr2/3全E完了は主張しない。
既二次／三次local余核の零を高次計算の根拠にしない。
次は指示に従いstable二源のwhole algebroid／normalizationへ戻り、同fixed-frame計算を反復しない。
