# Source-bound pair・monic化・norm・有限112被覆の独立監査

対象は同日付の `RAW_H_B_COPRIME_PAIR_PREFLIGHT`、`PAIR_MONIC_COMPLETION`、`PAIR_NORM_INPUT`、`PAIR_NORM`、`FINITE_CURVE_COVER` の五packetである。B35の独立レビューへsource hashを固定する。著者worker、GB、点、prime、未知因子のfactorizationは実行しない。

## 真のLCとaだけによるモニック化

選択二多項式は、実際の22 monic minorのgen10、gen15であり、68項と59項である。旧preflightの `d_leading_coefficient` はd⁴因子を残したtop部分であり、新しい真のLCと明確に区別する。両者を元多項式から照合する。真のLCをa³、a⁴で正規化したl0,l1と `q=a+13e/50-13/75` は `l0-q*l1=-11a/9` を満たす。

この恒等式から作る `H=-9(af-qg)/(11a⁵)` を、全係数で再確認する。Hはdのmonic4次、74項、aの最小指数-5であり、反転するのはaと非零QQ定数だけである。Hが(f,g)に入るため、pair quotientはQQ[a,a⁻¹,e]上で1,d,d²,d³により有限生成となる。自由性、正確なrank4、非空性はこれから従わない。

## 実norm行列と行列式

全四つのモニック除算恒等式

```
g*d^j = H*q_j + sum_i M_ij*d^i
```

を確認する。余りがd次数3以下なので、モニック除算の一意性から全16 entryは実際の乗法行列である。Mは計1284項。固定された `s=2035768318176269531250000*a^10` でN=sMを作ると、全係数が整数、指数も非負になる。

独自の共通分母付き整数疎多項式演算で、Nの第0行を除いた四つの3×3 cofactorを実際に計算する。第0行のLaplace展開で得た361項のRを保存normと照合し、さらに全四式 `Nv=Re0` を確認する。従って、Nvだけや次数上界だけから行列式を推測したものではない。

保存されたV、Qsumと、元のf,gに対する二cofactorを再構築し、`s*g*V-s*H*Qsum=R` と `R=uf*f+vg*g` を全係数で確認する。高いa/e次数も落とさず、d係数の相殺まで整数の積和で検査する。

## 既知単位の除去と有限被覆

全31段の掛け戻しと各除去因子によるU=adFBCJのliteralな割り切れを確認する。初め四因子はaC/100、残り27因子はaであり、未知の既約因子を選んでいない。最終225項hとUのexact gcdが1であることも検算する。その結果、正確に `R=content(R)*a^31*C^4*(h/10^8)` となる。

hのe最高係数は非零QQ定数である。E=h/LCe(h)を作り、eのmonic28次式として確認する。R_base=QQ[a,a⁻¹]上で、まずEによりfree28、次にmonic Hによりfree4となるので、

```
R_base[e,d]/(E,H)
```

は実際にfinite free rank112で、基底はe^i d^j（i<28,j<4）である。これはcover代数の自由性であり、pair quotientの自由性ではない。

FULLU上のB-bad点ではf=g=0、membershipからH=R=0、aCが単位なのでh=E=0となる。従って全bad点がこの被覆に含まれる。cover上にbad点が実在すること、cover全体がbadであること、全域H7やfull cubic/quartic rankは証明していない。

このfinite free coverにはdFBCJの逆元をまだ加えていない。追加のFULLU局所化がbase上でfiniteとなることや、special fiberからのNakayama適用を自動的に主張しない。全入力hashと実測値は専用review reportへ保存する。
