# 非被約m三：constant-root contractionによるwhole枝の事前排除

2026-11-07、producer純証明。
入力は `PLANE_M3_CONSTANT_ROOT_OBSTRUCTION.md` のrank二plane全source補題。
R被約性又はconnection logarithmic性はこの入力に要求されない。
単極residue比較のreduced-R強化は使わず、全rational expressionを直接比較する。

## 1. 全source必要式

constant normal root u、complement vでfullraw P0を

`P0^{i,u}=uX_i, P0^{i,v}=uZ_i+vY_i, P0^{u,v}=u(k20u+κv)`

と書く。fullP1のu-rowにおけるu非可除項はv²C、Cはconstant base vector。
C零なら全actual Dがu倍になり既約quinticに矛盾する。
actualD3のv³係数は3f q(C)、qはordinary complementary conormal character β22。
G3がu-rootを持つためq(C)零、かつC非零が必要である。

補助の全weight一projective equationは

`Ebase∧[(DY)C+κC]−[(2Y−X)∧C−P_{C(f)}]=0`。

これは全P1他係数/P2を落とした選択sliceではなく、原Jacobiのv³係数である。
本noteの空性では弱いscalar必要式q(C)零だけを使う。

## 2. degree(0,1)・Croot=bの全EXACT_Cを排除

`PLANE_NONREDUCED_M3_MOVING_SIMPLE_FIRST_CONNECTIONS.md` のCroot=bで、
h=a又はb、A=a又はcを取る。
EXACT_C全stratumはnormalizer前にδ非零で、

`q=ρa dloga+ρb dlogb+γd(b/a)+δd(c/a)`、

h=aでは(ρa,ρb)=(2,−1)、h=bでは(3,−2)。
Cbad=(Ca,Cb,Cc)としa²b q(Cbad)を展開すると

`Ca(ρa ab−γb²−δbc)+Cb(ρb a²+γab)+δCc ab=0`。

δ非零のbc係数でCa零、ρb非零のa²係数でCb零、ab係数でCc零。
従ってCbad零、既約quinticには矛盾する。
全A二配置・h二選択、計四whole EXTACT_C源はglobalK/operatorを計算する前に空である。
γを最初から零と仮定せず、δ非零の全stratumを直接排除した。

degree(0,1)の残finite sourceは四whole αと十一fixed。
これは全nineteen-source finite化の四源を純証明で閉じた結果で、他の源は残す。

## 3. degree(0,2)の全R2L源に対する必要方向

Croot=a²/abではqは両residues非零のclosed log formである。
Croot=b²でもq=4dloga−3dlogb+γd(b/a)、γ零/非零を保持。
同じ全polynomial比較から、いずれもCbad_a=Cbad_b零。
従ってnonzero bad vectorはC*∂cだけである。
ここからCbad零又は全G空を推論しない。
Aのc依存と残fullJacobi/actualDは別の義務である。

## 4. generic triple root・R3Lの二whole potentialを排除

`PLANE_NONREDUCED_M3_TRIPLE_ROOT_FIRST_CONNECTIONS.md` のq=dloga+2dφを使う。
φ=xy=bc/a²なら

`a³ q(Cbad)=Ca a²+2a(Cb c+Cc b)−4bc Ca`。

全coefficient比較はCa=Cb=Cc零を強制する。
φ=x²+y=(b²+ac)/a²なら

`a³ q(Cbad)=Ca a²+4ab Cb+2a² Cc−4b²Ca−2acCa`。

再び全三constant componentが零。
よってこの二whole potential枝はoffdiagonal r又はK/higherの値にかかわらず空。
normal Borelでrを零へ直す前の全r=λdφも含む。

残φ=0,x,x²の全r型（五＋三＋二＝十fixed）はまだ別義務。
特にφ=x,x²では必要方向Cbad∈C*∂cが残り、これだけで空としない。

## 5. 状態

本命題は既constant-root補題の直接適用によるpure geometryで、新operator計算はない。
source firstconnection有限化の別作者reviewはまだ受領していない。
全nonredm三又は全G/Eの完了を主張しない。
