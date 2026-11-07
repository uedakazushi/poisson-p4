# 型7 nil22長translation：全quartic支持の七正規源への還元

2026-11-07、PRIMARY E。全源有限化の厳密証明。成分帰属は次の課題であり、新成分は数えない。
cofactorは `Z=y∂x+∂y+w∂z`、全20係数源を使う。以下は複素数体上の被約支持についての議論。

## 全potential表示

補助座標 `t=y,u=x−y²/2,v=z−wy,w=w` のJacobianは−1で、逆も多項式である。
全六括弧を変換すると `P=∂t∧V+Jac(Φ)` で、

    Φ=−n17 u²−n18 uv+[n15+(n16−n13)w+(n19−n14)w²]u
       +(−3n11/4+n19/2+1/4)v²−(n6+n7w+n8w²)v
       +n3w+(n4−n0)w²/2+(n5−n1)w³/3−n2w4/4,
    Vu=(n11−2n14−1)u−(n12+3n8)v−n0−n1w−n2w²,
    Vv=−n18u+(−n11+n14)v+n9+n10w+n12w²,
    Vw=n13+n14w,  div V=−1.

全四Jacobiは `d(VΦ+Φ)=0` と同値である。従ってPoisson点では
`VΦ+Φ=k`、実Pfaffianは `F=Φ−k`、`VF=−F`。
全identityは 最終証明書の `potential_geometry` 節 で20係数のまま検証した。
この補助多項式座標は因子分解と次数を扱うためのもので、PGL同値とは混同しない。

q=n14、HをΦの(u,v)-Hessian、MをVの(u,v)-線形部とすると

    H=[−2n17,−n18; −n18,−3n11/2+n19+1/2],
    M=[n11−2q−1,−n12−3n8; −n18,−n11+q],
    MᵀH+HM+H=0, trace M=−q−1.

## Hessian rank2は全て空

`trace(adj(H)(MᵀH+HM+H))=−2q detH` なのでrank2ならq=0。
constant GL2で二次部分をUVへ直すと、conformal条件からMはdiag(α,−1−α)となる。
`F=UV+p(w)U+q1(w)V+r(w)`、Vw=s定数と書く。
UF/VFの係数を消すと、g=r−p q1は `s g'+g=0` を満たす。
非零多項式gは最高次数項によりこの式を満たせない。
従って `F=(U+q1)(V+p)` であり、全rank2は既約quarticを持たない。

## Rank1の有限resonance

constant GL2とactive変数の平方完成により
`F=U²+f(w)V+g(w)`、deg f≤2、deg g≤4とする。
active固有値は−1/2、passive固有値は−q−1/2である。

f≠0なら係数式は

    (n13+q w)f'+(1/2−q)f=0.

非零多項式fの次数d≤2から `q(d−1)+1/2=0`。
d=1は不可能であり、d=0,2はq=1/2,−1/2を与える。
q=0の定数flowでも非零多項式解はない。

f=0ならUの一次項を消す条件によりactive fieldも−U/2となり、
`(n13+q w)g'+g=0`。
g=0または偶数冪はFを平方差の積へするので、既約ならgはlinear flowの奇数冪d=1,3のみ。
従ってq=−1,−1/3。これでrank1の全quartic支持は四つのqだけへ覆われる。

## Rank0の有限resonance

`F=ℓ(w)·(u,v)+r(w)`、deg ℓ≤2、deg r≤4。
q≠0でlinear flowの根をw0とする。ℓ(w0)=0なら定数項のsemi-invariant式からr(w0)=0も従い、
Fはw−w0で割れてしまう。従って既約ならℓ(w0)≠0で、Mは固有値−1を持つ。
flowを原点に移してℓ=Σℓ_d w^dと書くと
`(Mᵀ+(1+dq)I)ℓ_d=0`。
非定数係数ℓ_dがあればd≤2、異なる固有値とtrace M=−q−1から
`q(d−1)=−1`。d=1は不可能、d=2でq=−1となる。
ℓが定数ベクトルだけなら、quartic条件はrのw4項を要求し、
定数項式の最高次数から4q+1=0、すなわちq=−1/4。
q=0では、flow非零ならℓは単純固有値−1の定数ベクトル、rは次数≤2となる。
flowも零ならFは多項式f(w)とlinear fiber式の積となり、既約quarticは出ない。

実coreにH=0とdet(M+I)=0を課すと、q=−1/4は単位理想。
q=−1は正確に `n11=0,n19=−1/2,n12=n8,n2=4n8²` のみとなる。

## 実際のaffine正規化と七core

次は補助座標変換ではなく、実際のPGL変換である。

    h'=h, x'=x+a z, y'=y+a w, z'=r z, w'=r w, r≠0.

全cofactor Nと可換、det=r²。最終証明書の `actual_core_action` 節 は全90係数から全20parameter mapを抽出し、
全六括弧の普遍pushforwardを一致させた。特にqは不変である。

rank1でn17≠0ならa=n18/(2n17)とr²=1/n17でHをdiag(−2,0)へ正規化する。
conformal式が `n11=2q+1/2,n12=−3n8,n19=3q+1/4` を強制する。
四つのqを実coreに戻すとq=−1/2だけで、n2=−3/16,n8=±1/4。
r=−1で二符号が交換されるので正符号一つで覆う。

rank1でn17=0ならdetH=0からn18=0、H22≠0。
conformal式よりn11=q+1/2。q≠−1/2ではn19=qとなり、
q=1/2,−1,−1/3の三つのcore曲線を `a=−4n8/(2q+1),r=1` でn2=n8=n12=0へ送れる。

q=−1/2ではcoreは正確に

    n11=n17=n18=0,
    (3n19+1)n8=0,
    4n12n8+(4n19+1)n2=0.

Poisson支持上で最高quarticはΦの最高次数と等しく、ここでは `w³(n8 y−n2 w/4)`。
n8=n2=0の全連続coreはquarticを持たない。
n8≠0ならn19=−1/3,n2=12n12n8であり、`a=−3n12,r=n8` でn8'=1,n12'=n2'=0。
n8=0,n2≠0ならn19=−1/4であり、n2=ρ²という全複素非零値を覆う表示を使い、
`a=−4n12,r=ρ` でn2'=1,n8'=n12'=0。
rank0の残curveは `a=−2n8,r=1` でn2=n8=n12=0へ送れる。

これらは 最終証明書の `rank_resonance` 節 に全core理想の両方向GB所属と全normalizationを保存した。
得られる七coreを(n2,n8,n11,n12,q,n17,n18,n19)で書くと

|名前|core|
|---|---|
|rank1 U|(-3/16,1/4,-1/2,-3/4,-1/2,1,0,-5/4)|
|rank1 V、q=1/2|(0,0,1,0,1/2,0,0,1/2)|
|rank1 V、q=−1|(0,0,-1/2,0,-1,0,0,-1)|
|rank1 V、q=−1/3|(0,0,1/6,0,-1/3,0,0,-1/3)|
|rank1 V、yw3型|(0,1,0,0,-1/2,0,0,-1/3)|
|rank1 V、w4型|(1,0,0,0,-1/2,0,0,-1/4)|
|rank0|(0,0,0,0,-1,0,0,-1/2)|

全qは非零なので、もう一つの実affine変換
`z'=z+s y,w'=w+s h`、s=n13/qでn13'=0へ移せる。
全coreは不変であり、この変換の全六pushforwardも 最終証明書の `all_seven_complete_families` 節 に保存した。

## 全fiberと既約quartic opens

各core+n13=0で元の全Jacobi理想を完全に解いた。
各parametric理想との両方向所属、全括弧・Pfを保存した。残る七つのparameter空間は
A4、D(n15)⊂A5、D(n3−4n6n7)⊂A5、D(n5−n1)⊂A5、A4、A4、D(n6)⊂A5である。
unitの全零補集合は平方Pfまたはw因子で排除した。
rank0の生の残compatibility n6*n10=0も保持し、n6=0全枝をw因子で排除した後にn10=0とした。
他の全parameterでPfのprimitive linear fiber又はodd-power discriminantを検証し、全Fが既約quarticである。

各正規源に逆core正規化とw-centeringのparameterを戻した固定cofactor源は次元7以下。
既済のaffine cofactor orbit12とH4−scalar1から、各全射影飽和閉包は次元22以下である。
各源は既約なので、今後一つの正規anchorが実orbit22と分かればその全飽和閉包はそのorbitclosureとなる。
成分帰属は未解決であり、この有限化や上界だけから新成分を数えない。

## 再現と独立確認

全証明書は `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_SOURCES_20261107.json`。
reusable runnerの `--nil22-long-quartic-sources` と有限契約
`jobs/hyperplane_quartic_nil22_long_quartic_sources.json` で再現する。
全phaseは元の公開nilpotent源表だけから再計算し、work内の探索結果を入力にしない。
独立Fraction/Laurent実装は全数学場合分けと全90source・potential・全actual parameter map・
core/七fiber理想の両包含・全Pf/単位零補集合・cofactor軌道12をPASSした。
公開bundleだけに入力を限定した独立V2
`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_REVIEW_V2_20261107.json`
も全PASS（0.720秒）である。成分帰属の検証とは別に扱う。
