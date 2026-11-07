# 被約退化因子のrank零曲面では、全重複度でconormal Lie代数は可換

2026-11-07。producer `root`。主義務G。
以下は新しい局所幾何証明であり、旧first-jet証明書のbytesを変更しない。
別作者 `generic_exhaustion` は証明の各段階を独立に確認した。
独立確認記録とともに採録するまでは、既存のm二限定主張を上書きしない。

**命題。** smooth complex fourfold上のgenerically symplectic Poisson tensorをpiとし、
その実Pfaffian退化因子Dが被約とする。
Sをpiのrank零locusに含まれるintegral surfaceとする。
Sのsmooth一般点で、ideal Iのconormal Lie algebra `I/I^2` は可換である。
Sに沿うDの位数二、Frobenius first Lie algebra又はordinary doubleを仮定しない。
特にSがplaneなら、この結論はDの全可能位数に適用出来る。

## 1. rank零locusにhypersurfaceは無い

smooth ambientのcodimension一primeの一般点では、そのidealは一元gで生成される。
もしrank零locusがそのhypersurfaceを含めば、piの全局所係数はgで割れ、`pi=g*rho`。
四次元では`pi wedge pi=g^2*(rho wedge rho)`なので、実Dの同hypersurfaceへの位数は二以上。
これは被約性に反する。
従ってrank零locusはdimension二以下であり、Sはその被約台の既約成分である。
Sの一般点では他の成分と交わらないsmooth openを選べる。

## 2. 局所modular vector fieldはSへtangent

局所nonvanishing volume formを選び、Hamiltonian field `X_f` に対して
`Z(f)=div(X_f)` と置く。
HamiltonianのLeibniz則とdivergenceの積公式では交差項が相殺するので、Zはvector field。
Jacobiと

`div([X_f,X_g])=X_f(div X_g)-X_g(div X_f)`

からZはPoisson derivation、すなわち`L_Z pi=0`である。
従ってZの局所holomorphic flowはpiとそのrank零被約台を保つ。
有限個の既約成分をconnected flowが置換することは出来ず、上記Sのsmooth openでZはSへtangent。
同じことは`L_Z pi=0`からzero idealを保存し、標数零でradicalにも降りることからも従う。
volumeをunit倍したときのZの差はHamiltonian fieldで、piが零となるS上では零である。
従って今回必要なS上のnormal成分はvolumeの選択によらない。
これはprojective space上の大域modular lift又はそのstabilizerを構成する議論ではない。

## 3. normal成分の全一次係数

Sのsmooth一般点のholomorphic coordinatesを`(a,b,z,w)`、`I=(z,w)`とする。
全pi係数はIに入り、特に

`{z,w}=l_1(a,b) z+l_2(a,b) w mod I^2`。

coordinate volumeで`Z^j=sum_i partial_i pi^{ij}`という符号規約を使う。
mixed係数pi^{az},pi^{bz},pi^{aw},pi^{bw}はSで恒等零なので、そのbase方向微分もSで零。
従って全normal成分は正確に

`Z^z|S=-l_2`、`Z^w|S=l_1`。

§2のtangencyにより両者とも零となる。
I内の任意の二関数のbracketについて、Leibniz則のmixed項は二つのI因子を持つ。
従って`{I,I} subset I^2`、すなわちconormal Lie代数は可換である。
初めはSのdense openでの計算だが、Sがintegralなので、正則なconormal bracketは全smooth locusで零となる。

## 4. 分類への正確な帰結

P4のrank零planeに適用すると、normal-normal quadratic coefficientのnormal次数一の項は
plane全体で零となる。
従って旧m二入口のnormal weight零35、一30、二9というEuler quotientの全sourceを、
m三以上にも使用出来る。traceless co-Higgs共通lineのO零/O一分岐と全weight零Lax式も、
その導出に使った可換conormal/Jacobiだけで同様に有効である。

これは旧m二のnormal-form表が、そのままm三以上を尽くすという結論ではない。
例えば旧表で除外したtau零、leading-Pfaffian零、rank低下又はresonant row零を
すべて復活させ、全higher operator、source dimensionと補集合を別に処理する必要がある。
全G/Eの尽くし、新成分又は位数三以上の不存在は主張しない。
