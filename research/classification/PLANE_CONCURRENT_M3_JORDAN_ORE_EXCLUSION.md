# K零 binary Ore：全 primitive Jordan 源の actual divisor 排除

2026-11-08、producer 純証明・独立 review 前。
新 CAS、matrix、minor、jet、source probe 又は数値点は用いない。
全 source 同値は `PLANE_CONCURRENT_M3_STRICT_ORE_REDUCTION.md`、
全 actual normalizer は `PLANE_CONCURRENT_M3_BINARY_ORE_SOURCES.md` §7 と
`PLANE_CONCURRENT_M3_BINARY_ORE_CENTRAL_NORMALFORMS.md` §4 を用いる。
本稿は全 Jordan source を扱い、SS／constant direction は別の全源証明へ戻す。
元 binary cubic f が reduced であることすら本排除には不要である。

## 1. 全 central coefficient を保持した linear-field 表示

`y=(a,b,u,v), D0=u∂a+v∂b`、
`C2=c0u²+c1uv+c2v²`、`N_C=−c2v∂a+(c0u+c1v)∂b` とする。
全 Jordan symbol は正確に `Q=D0∧S`、

`S=−(1+s)a∂a+(a−(1+s)b)∂b+u∂u+(u+v)∂v+N_C`。

実際 mixed rows は
`U=(u²,u²+uv), V=(uv,uv+v²)`、
`Qab=a[u+(1+s)v]−b(1+s)u+C2`、Quv 零で、元全 source に一致する。
`[D0,S]=−(s+2)D0` なので strict Q Jacobi も満たす。
元 full tensor は `Π=Q+∂c∧W`、W は全 quadratic field を含む源。
元 Jacobi の残全条件は `E_y∧L_W Q=0`。
W の元 f、normal connection、full higher を chosen slice に固定しない。

## 2. chart の vertical 項を保持する module identity

u 非零 chart で

`x=a/u, r=v/u, h=b/u−av/u²`

と置く。u一で inverse は a=x,b=h+rx,v=r、全て polynomial。
Euler は `E_y=u∂u`、Q は

`Q=π+u∂u∧Z, π=∂x∧H, Z=−∂x`。

ここで

`H=[μh+C2(1,r)]∂h+∂r, μ=−s−2`。

Z の sign は Q(u,a/u)=−u から直接決まる。
quadratic W は homogeneous weight一なので
`W=u(Y+σu∂u)` と書け、Y は x,h,r の polynomial vector field。
その horizontal module は正確に

`[Y,π]+Y∧∂x=0`。                                    (1)

実際 `[uY,u∂u]=−uY` と Z=−∂x により最後の項が出る。
twisted W を普通の Poisson vector field と扱ってこの項を落とさない。

`F=Hh Yr−Hr Yh` とすると、(1) の hr component と xh/xr components から

`∂xF=0, H(F)=(div H−1)F`。                           (2)

第一式は `Hh∂xYr−Hr∂xYh=0`。
第二式は `Hh*(1)xr−Hr*(1)xh=(div H−1)F−H(F)=0`。
これらは全 W に対する必要条件であり、低次数又は特定の W の推測ではない。

actual divisor は u一 chart で overall sign を除き F。
座標変換の Jacobian は u³、E_y と W の係数を戻すと global quintic は

`D=u⁵ F(a/u, b/u−av/u², v/u)`。                     (3)

以下ではこの全 homogeneous polynomial の必要性も保持する。

## 3. s≠−2：全 polynomial eigenfunction と global factor

actual central shift map は s≠−2 で rank三なので全 C2 を零へ戻せる。
inverse は元 source を保つ positive base shift であり、全 W source を同時に移す。
従って

`H=μh∂h+∂r, div H=μ`。

(2) は x 非依存の polynomial F に対して

`(∂r+μh∂h)F=(μ−1)F`。

μ非零なら `F=Σ h^j f_j(r)` の各 coefficient は
`f_j'=(μ−1−μj)f_j`。
非零 polynomial がこの constant-coefficient ODE を満たすには右の定数は零であり、
f_j は constant となる。μ非零なのでそのような整数 j≥0 は高々一つ。
従って `F=λh^j`、必要条件は `μ(j−1)=−1`。
もし F が零なら D も零である。

global formula (3) は

`D=λu^(5−2j)(bu−av)^j`。

bu−av は u で割れないため、λ非零の global polynomiality は j≤2 を要求する。
j一は eigencondition に不可能。
残 j零は λu⁵、j二は λu(bu−av)²で、どちらも proper factor を持つ。
この議論は全 s の値を扱い、整数 resonance を generic parameter の外へ捨てない。
μ零なら (2) は F_r=−F なので polynomial F は零。

## 4. s=−2：全 central residue も actual D 零

この rank-drop fiber の全 residue は `C2=γv²`。
γ零／非零の両方を残す。H は

`H=∂r+γr²∂h, div H=0`。

(2) は `H(F)=−F`。
polynomial change `h'=h−γr³/3` は polynomial inverse を持ち、H を ∂r にする。
ゆえに polynomial F の ∂rF=−F は零のみ。
この coordinate change は chart 内で ODE を解くための代数的変換であり、
global PGL normalizer 又は source 同値だとは主張しない。
全 critical central source の actual D は零となる。

## 5. 結論と全 source への橋

全 s、全 zero/nonzero central residue、全 quadratic W に対し、
actual D は零、λu⁵ 又は λu(bu−av)²。
従って **全 primitive Jordan K0 source の reduced irreducible quintic は空**。
global actual normalizer の inverse と full W を保持しているので、
canonical chart の一点又は generic period だけに対する結論ではない。
global 原 tensor に同じ polynomial 因子の literal式が残るとは限らないが、
actual invertible coordinate change は divisor の可約性を保つので元全 source も排除される。

constant-proportional U,V は別の essential-three 八成分へ戻す。
SS、primitive constant direction、K非零及び全 concurrent G の受理は各別証明と独立 review で統合する。
本稿だけで全 E、全 G 又は全 classification の完了とは呼ばない。
