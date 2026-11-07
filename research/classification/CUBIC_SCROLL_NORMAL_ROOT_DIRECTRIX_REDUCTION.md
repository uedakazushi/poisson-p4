# smooth cubic scroll の全 normal root と負の directrix への必要還元

2026-11-07。主義務 G。これは新しい幾何 producer の証明であり、別作者 review 前である。
全 rank 二 scroll 源の非優越又は空性をまだ主張しない。主原稿・既存証明書は変更しない。

## 1. 仮定と結論

`S=F1=S(1,2)⊂P4`、`H=C0+2f`、`L=C0+f` とする。
実退化因子 `D` が既約被約 quintic で、`S⊂Sing D`、`pi|S` が一般 rank 二である全源を扱う。
零因子 `R=div(pi|S)∈|2C0+3f|` は可約・非被約を許す。

**必要還元。** この全源では `C0⊂R` でなければならない。
従って元の全21次元 canonical linear source に、restriction の `C0` 上の二係数が零という
全域線形必要条件を加えてよい。これは一般元の open だけの条件ではない。
残る `R` の全てを排除した、又は固定 source dimension 七以下を得たとの主張ではない。

既済 `CUBIC_SCROLL_RANK_TWO_INTEGRAL_CURVE_EXCLUSION.md` の §1 により、全 nonplane `S` で
`D` の一般 normal multiplicity は二、conormal bracket は可換である。
以下ではその結果と、同 note の parallel normal quadratic の構成だけを使う。

## 2. 正規束の全域表示

`q0=x0x3−x1x2, q1=x0x4−x1x3, q2=x2x4−x3²` とする。
Hilbert–Burch の二つの syzygy は

`(x3,−x2,x0)`, `(x4,−x3,x1)`。

全 parametrization `x=(su,tu,s²v,stv,t²v)` への制限では、それぞれ
`s(tv,−sv,u)`、`t(tv,−sv,u)` となる。従って conormal の全 presentation は

`0 → O_S(−L) --(tv,−sv,u)--> O_S³ → N*(2H) → 0`。

三成分は同時に零とならない。写像 `beta:S→P2`, `beta=[tv:−sv:u]` は `C0` を
`p=[0:0:1]` へ縮約する blowdown で、`beta*O(1)=L`。
`Q` を `P2` の universal quotient bundle とすれば、全域で

`N*(2H)=beta*Q`, `N=beta*Q* ⊗ O_S(2H)`。                 (1)

`c1 N=4H−L=3C0+7f`、`c2 N=4H²−2HL+L²=9` とも一致するが、
この数値一致を bundle 同型の代りには使っていない。

## 3. primitive root の下界と a=1 の全場合

normal quadratic の height 一の各 primitive linear factor を
`ell:N→A`, `A=aC0+bf` と書く。codimension 二の零点は許し、divisorial zero は許さない。
`N|f=O(1)⊕O(2)` と `N|C0=O(2)²` への非零制限から

`a≥1`, `b−a≥2`。                                        (2)

`C0` 上で全写像が零なら `C0` が共通 divisorial zero となるため、その場合を (2) の外へ
捨ててはいない。(1) から

`Hom(N,A)=H0(P2,Q(b−4)⊗I_p^(b−a−2))`。

`H0(Q(d))=0` for `d<0` は Euler sequence から従う。従って更に

`b≥4`。                                                 (3)

`Q` の constant-vector section `p` は `p` で単純零となる。
局所座標 `(X,Y)` では `(X,Y)` であり、blowup の各 chart で exceptional factor を割ると
`(1,t)` 又は `(s,1)` である。従って `beta*Q(−C0)` に nowhere-zero section があり、dual/twist により

`0 → O(2C0+3f) → N --ell0--> O(C0+4f) → 0`。            (4)

任意の `a=1` root について、(4) の左辺から `A=C0+bf` への写像は ruling 上の degree が
`1−2=−1` なので零である。従って全写像は `ell0` と
`H0(O((b−4)f))` の積に一意に分かれる。`b>4` なら fiber の divisorial factor が必ず生じる。
よって **primitive a=1 root は `A0=C0+4f` と `ell0` の定数倍だけ**である。

## 4. ambient quadric root は全て排除される

`A=2H` の root space は `H0(Q)=C³` であり、これは全 ambient quadrics
`Q_c=c0 q0+c1 q1+c2 q2` の normal derivative に正確に一致する。
`c0` 又は `c1` が非零なら `Q_c` の rank は四である。
実際、対応する四変数 principal minor はそれぞれ `c0^4/16`, `c1^4/16`。
残る `Q_c=q2` は `C0={x2=x3=x4=0}` に沿って normal derivative が零となるため、
`N→2H` としては primitive でない。その primitive saturation が (4) の `ell0` である。

rank 四の `Q=Q_c` が normal quadratic `G` の一根なら、`D|Q` は `3S` を含む。
これは generic smooth point の二つの ambient normal coordinates `(q,z)` で直接分かる。
`D∈(q,z)²` であり、二次 normal form の `z²` 係数が零なので
`D∈(q)+(q,z)³`、従って `D|Q∈I_(S/Q)³` である。
`Q` の vertex は高々一点なので、この generic Weil-divisor multiplicity の議論から漏れない。

rank 四 quadric threefold の `Cl(Q)=Z²`, `H=(1,1)` で、integral nonplane cubic scroll の class は
`(1,2)` 又は `(2,1)` である。degree 三の他の非負候補 `(0,3),(3,0)` は ruling planes の和で
integral ではない。従って `D|Q≠0` なら residual effective class は

`5(1,1)−3(1,2)=(2,−1)` 又は `(-1,2)`

となり不可能。よって `D|Q=0`、すなわち `Q|D` となり既約性に矛盾する。
従って **primitive normal root `A=2H` は rank 一・二とも全排除**される。

## 5. rank 二 normal quadratic の全八 class と残る五 pair

`sigma=pi|S` とする。parallel 性により `disc(G)=c sigma²`。
`c≠0` なら function field 上の二根は分離しており、primitive factorization を

`G=gamma ell1 ell2`, `V=div(gamma)`, `Delta=div(ell1 wedge ell2)`

と書く。全域の divisor 等式は

`V+Delta=R`, `V=5H−A1−A2`, `Delta=A1+A2−det N`。        (5)

`Ai=ai C0+bi f`、`a=a1+a2,b=b1+b2` とすると、effective cone と (2),(3) から
`3≤a≤5, 8≤b≤10, b−a≥4`。従って全候補は

`(a,b)=(3,8),(3,9),(3,10),(4,8),(4,9),(4,10),(5,9),(5,10)`。

**注意。** `V` 又は `Delta` は `C0` を含める。従ってその `C0` との交点を非負と仮定して
候補を減らしてはいけない。この全八組には負交点の補集合も含めている。

`V=0` は `R` が reducible/nonreduced でも排除できる。
`Bl_S P4` の exceptional divisor `E=P(N)` 上で inverse form の closed residue `alpha` は、
二つの horizontal root prime だけに pole を持ち、その residues は非零定数 `lambda,−lambda`。
vertical pole は `V=0` により存在しない。
既済 note §3.1 と同じ一般 ample complete-intersection curves の residue 定理から
`A1≡A2`。`Pic(F1)=Z²` なので `A1=A2`。
しかし (5) は `A1+A2=5H` を要求し、`C0` 係数五の parity に反する。

(4)、ambient-root 排除、`V=0` 排除を全八組へ適用すると、残る unordered primitive pair は正確に

| `(a1,b1),(a2,b2)` | `V` | `Delta` |
|---|---|---|
| `(1,4),(2,5)` | `2C0+f` | `2f` |
| `(1,4),(2,6)` | `2C0` | `3f` |
| `(1,4),(3,5)` | `C0+f` | `C0+2f` |
| `(1,4),(3,6)` | `C0` | `C0+3f` |
| `(2,5),(2,5)` | `C0` | `C0+3f` |

第一行の `2C0+f` は `C0` を固定因子に持つ。他の明白な三行も同様である。
第三行で `V∈|C0+f|` が prime の場合だけを最後に扱う。
horizontal residues と唯一の vertical prime residue `mu` の Gysin class relation は

`lambda(A1−A2)+mu[V]=0`。

これは generic pole が単純であることと closed 性から従い、全点 SNC は要求しない。
`A1−A2=−2C0−f` と `[V]=C0+f` は `Pic(F1)⊗C` で独立、`lambda≠0` なので矛盾する。
従って `V` は非 prime。`|C0+f|` の非 prime member は `C0+f0`（一 ruling fiber）だけである。
よって第三行も `C0⊂V⊂R` となる。

## 6. rank 一 normal quadratic の全補集合

`c=0` なら `G=gamma ell²`、`V=5H−2A` は effective で `R` に支持される。
(2),(3) から `1≤a≤2,4≤b≤5`。
`a=1,b=5` は (4) により非 primitive、`a=2,b=4` は ambient quadric root で排除された。
従って残りは

| `A` | `V` |
|---|---|
| `C0+4f` | `3C0+2f` |
| `2C0+5f` | `C0` |

`3C0+2f` の全 effective member は `C0` を固定因子に持つ。
どちらも `C0⊂V⊂R` である。`G=0` は一般 normal multiplicity 二に反する。
これで rank 一・二の全場合を尽くし、§1 の必要還元が従う。

## 7. 次の有限検証と未完了部分

全21 source の `C0` 上で `sigma` は `O_(C0)(1)` の section である。
従って必要条件 `C0⊂R` は二つの global linear coefficient equations として、
元の全150係数・全350 Jacobi・rank 二 actual Aut patch へ戻せる。
次の計算はこの全条件の具体的復元だけを最初の目標とし、generic minor の再反復をしない。
その kernel に高 rank point が残れば、元の nonlinear whole zero/nonzero case proof 又は
実 `D` の因子を使う。二つの旧弱上界及び三枝の rank2/2/3 は成功へ読み替えない。

この note は smooth scroll に限る。twisted cubic cone、その vertex、全 cubic surfaces の G、
既知成分への E 帰属及び成分数は未解決のままである。
