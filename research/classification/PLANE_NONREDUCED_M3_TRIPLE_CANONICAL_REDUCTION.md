# generic triple finite正常形を九十源へ減らすcanonical obstruction

2026-11-08、producer純証明。
共通lemma `PLANE_M3_CANONICAL_TRANSVERSE_OFFDIAGONAL_OBSTRUCTION.md` は別作者手確認済み。
元の全有限入口 `PLANE_NONREDUCED_M3_TRIPLE_FINITE_NORMALFORMS.md` のbytesと全classを保持し、
以下の十三fixedだけを全source空として排除する。

## 1. R2L log-onlyでbad vectorは常にc方向

p=(dlogf)/3+αdlog(a/b)、q=2p−dlog(f/h)。
X=3αa∂c、Y=(3e−4+6α)a∂c。
qの両residue非零なら全rational coefficient比較でCbad_a=Cbad_b零。
q=dlogaなら追加b方向が残る可能性があるが、fullw³式にCb非零を代入すると
Xc係数3α=4が必要。q=dloga自体はα=5/6−e/2なのでe=0,1,2で一致しない。
q=dlogbなら追加a方向についてfullw³式は3α=−8を要求する。
q=dlogb自体はα=1/3−e/2なのでこれも一致しない。
従って全αでCbad∈C*∂cであり、同じ式からκ零が必要。

γd(b/a)非零源でもq(Cbad)のb²/ab係数から直ちに同じc方向となりκ零。
X,Yの像はc方向のまま。
δd(c/a)非零は別noteで全排除済み。

## 2. R2L log二resonanceの十二fixedを排除

共通canonical lemmaを全r八dimへ適用するとE=H零が必要。
η=1、r=dψ/xのψ=y²,y²+x²ではH非零なので全空。
η=−2、r=x²d(ψ/x²)のψ=y²,y²+xではE非零なので全空。
各eの全二＋二型、計十二fixedを排除する。
fullhigher variablesを固定せず、canonical strict weight一の必要係数だけを使う。

## 3. R3L一fixedを排除

φ零、r=d(xy)はX=Y=K零、Z=(0,3b,−3c)、P0自身divergence零。
共通lemma§6のstrict vector式からCbad零となり、全source空。
他のφ零ψ型、φ=x/x²のr型をこの排除へ含めない。

## 4. 全残sourceの正確な数

R2L γ一は各eでr0/dlogx二whole line、η一の追加一fixed。
R2L log-onlyは各eでgeneric三whole lineと
η−1四、η零六、η一六、η−2五、計二十一fixed。
R3Lは旧十からψxy一を除き九fixed。
したがって残全first symbolsは

`15 whole α sources + 75 fixed sources = 90 families`。

whole log sourceはη∈{−2,−1,0,1}を除くが、その全fiberは残二十一fixedへ入る。
γ一whole sourceは全αを保ち、η一にだけ追加offdiagonal型がある。
全normalizer及び全inverseは元有限正常形noteに明記済み。
これを新しい標本/point listとは解釈しない。

次は全九十first symbolsでglobalK/P0/全30/9/actualD3lineを一bounded batchで扱い、
D3像零、fullkernelのbad三行零（proper z factor）、又はB+k1+k2+5≤24を先に判定する。
全false源は元全supportを保存し、必要な場合だけOre/actualDへ移る。
全generic triple、全m三、G/Eの完了はまだ未証明。
