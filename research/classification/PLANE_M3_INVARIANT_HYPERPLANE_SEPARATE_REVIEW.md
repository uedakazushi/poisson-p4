# Plane m三の不変 hyperplane 四補完：全源・二 critical の独立確認

2026-11-08、PRIMARY `E_reduced_nonnormal`。**独立 PASS。**
新しい source 支持を dimension 上界から推論せず、保存された全必要核と元 Jacobi の二乗式から
実五次の一次因子を示す。新成分数は零。

## 1. 全 source と元 row 順序を先に固定

座標 `(a,b,c,z,w)`。全 P1 三十係数の順序は

```
0..8: zL、9..17: wM；base pair は ab,ac,bc、各係数 a,b,c。
18..20: D_i*z² * ∂i∧∂w,
21..23: B_i*zw * ∂i∧∂z,
24..26: E_i*zw * ∂i∧∂w,
27..29: C_i*w² * ∂i∧∂z  (i=a,b,c)。
```

元 full_source1 の十括弧をこの独立 recipe と全係数で照合した。
FULL P2九は base pairs ab/ac/bc それぞれの z²,zw,w² だけで、normal column は全て零。
affine particular を含めて normal row を加えない。
homogeneous kernel 零を P2零と取り違えない。

旧独立 receipt の complete-kernel・全 complex walls・全 source/normalizer inverse を利用する。
今回、新しい M1、rank、operator、minor選択、point、profile、jetは作らない。

任意の skew quadratic P で v-row が v 倍なら、
actual signed Euler-Pfaffian `D=Σ(-1)^i*x_i*Pf(P_omit_i)` は v 倍。
v-row を残す項はその entry から、除く項は外側の x_v=v から因子を持つ。
Poisson条件や generic parameter 仮定すら不要である。

## 2. 直接 normal-row 読み戻しで閉じる全源

元 complete K1 を保存 source から読み、短い新 certificate の全 row と直接照合した。

| 対象 | 全 record | 確認した結論 |
|---|---:|---|
| C-exact e=0,1,2、generic 三＋全八壁 | 11 | C又はD三行が全列零、元P0の対応rowも可除 |
| squarefree C_DEPENDENT_A2_G2 | 1 | 元30×10 K1のC三行が全列零、P0のz-rowはz倍 |
| moving-double 弱二源＋全壁 | 10 | critical−8/3以外の全C三行が零 |
| moving-simple 初期五family＋全壁 | 17 | critical−11/3以外の全C三行が零 |

C-exact は e0/τ1/3 で w、e1/τ4/3 で z、他は zw の因子を得る。
各 generic family の全保存 complex-wall list と、今回の readback list が正確に一致することも確認。
追加 kernel のある壁を generic kernel で置き換えない。
元P0の normal-row 可除性と FULLP2九の normal零を併せ、
全 necessary supersource 上で一つの不変 normal hyperplane が得られる。
従って当該全解の D は零又は proper linear factor を持つ。

## 3. 残る二 critical：元全350と一意 inhomogeneous P2

二 critical の全 K1 はともに九次元 p0,…,p8、
唯一の bad normal coefficient は `P_cz|z0=p8*w²`。
元 complete30×9 kernel と上記三十方向 recipe から P1 を独立に組み、保存十括弧と一致した。
P0も元 exceptional record の全 source と直接一致する。

FULL P2九は r0,…,r8 のまま保持し、weight二式 B*r+C(p)=0 の
**同じ旧非零九次 minor** による一意 r_part(p) を利用する。
旧 rank九/nullity零 certificate の row、column、determinant が新 report と一致し、
元選択九式へ戻して全 polynomial が零となることを独立に確認した。
未知 source 係数の反転、P2零の仮定又は新 minor の選別は無い。

さらに、全350 weight二 compatibility、P0+P1+P2_part の元全350 Jacobi、
定義からの full signed Euler-Pfaffian を Fraction で別に再構成した。
元 row order は quad0123,0124,0134,0234,1234、各(a,b,c,z,w) degree四の降順全70 monomial。
Schouten bilinear が二倍となる規約を避け、独立実装は single cyclic を直接使う。

| 全 critical source | 元 row54 と row279 | 全 reduced point support |
|---|---|---|
| moving-double REPEATED_LINE_OTHER、τ=-8/3 | `-2*p8²` | p8=0 |
| moving-simple DEG1_HB_CA_OTHER、τ=-11/3 | `p8²` | p8=0 |

全元350 coefficient が保存式と一致し、二行だけを切り取って元系と同値と主張したものではない。
C 上の Poisson point では p8=0、その後は FULLP2 particular を含む全 z-row が z 倍。
full actual D の一致と p8=0後の z可除性も直接確認した。
leading amplitude を逆元にしないので、元の leading零補集合もこの因子結論に残る。
非被約 scheme 上の p8 の消滅又は正確な nilpotence exponent は主張しない。

## 4. 旧全源・全 module を戻した pointwise aggregate

旧独立 `NONREDUCED_M3_CONSTANT_SEPARATE_REVIEW.md` は十四全 first-symbol、全α壁、
残六 Euler/Ore fullsource、二十二 canonical module の全 affine W と actual逆を保持する。
旧未閉部分は C-exact三源だけであり、今回§2で全空となる。
残 constant-vector factor 分岐では D が c と別の constant vector 方向に不変となり essential≤3。
被約既約 quintic は essential二以下になれず、既 canonical/Saito＋全八曲線分類＋全磁気 inverse により
**constant-root(0,0)全 stratum は既八 essential-three 成分へ点ごとに入る**。
旧G次元≤23だけから帰属を出していない。

旧 moving-double 全十一 source のうち、七 source は元 leading像零、二 Ore source は
全 module の proper factor又は不整合で既独立受理済み。
残二弱源を今回の全 normal-row と critical square で閉じるので、
**moving-double(1,0)の指定全 stratum は空**。

旧 moving-simple は全三十七 first-symbol と全二十二壁、全 actual逆、
四残 Ore fullsource と十 module を既独立受理済み。
今回全十七初期 record が空となり、他初期枝の leading零／global K不整合を戻すと、
**moving-simple(0,1)/(0,2)全 stratum は HB_0/τ0 と HA_1/τ2/15 の二 canonical source に含まれる**。
両 source の actual component 帰属は別の証明に委ね、ここで空とは呼ばない。

squarefree C_DEPENDENT_A2_G2 の最後の全枝も今回§2で空に強められる。
generic triple、normal m二/四、被約R、他の source 又は全分類完了は本 receipt の結論ではない。

## 5. 有限実行と保存

実行前 gate PASS、ONE Fraction batch、**0.414秒／RSS115920KiB／1559恒等式**。
producer、shared6022 runner、CASを import せず、旧独立多項式 primitives のみ再利用した。
同じ old M2 certificate の照合であり、新 M1や rank/point/profile は無い。

`reports/CLASSIFICATION_PLANE_M3_INVARIANT_HYPERPLANE_SEPARATE_REVIEW_20261108.json` に
全 normal-row readbacks、二 critical の full-source確認、独立 inline source と全 input SHAを保存。
実行契約は `jobs/plane_m3_invariant_hyperplane_separate_review.json`。
producer/source/currentjob/主稿/global/Git は変更していない。
