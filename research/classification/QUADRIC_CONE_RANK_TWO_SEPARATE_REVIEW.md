# rank 二 quadric cone 全支持の独立監査：射影像二十四上界

2026-11-07。主義務 G。別作者による全源の幾何読解と、既存の独立 Fraction
ライブラリによる全係数再構成を行い、**全項目 PASS** を得た。producer、SymPy、
native CAS は import も実行もしていない。主原稿・global 台帳・Git は変更しない。

## 1. 結論と範囲

複素 P4 上の rank 四 Poisson 構造の実退化因子 D が既約被約 quintic であり、
その非正規曲面支持に quadric cone S が現れ、S 上の一般 restriction rank が二とする。
この**全支持**の PGL5 飽和射影像は二十四次元以下である。
S 上の零 divisor R の既約性・被約性、ordinary node、vertical divisor の不存在を
仮定しない。未知二十五〜二十七次元 G 成分の一般元から本枝を除外する結論である。
全点の E 帰属、存在・不存在、成分性、scheme 構造、全 G の尽くしは主張しない。

独立 code は `code/classification_quadric_cone_rank_two_review.py`、有限契約は
`jobs/quadric_cone_rank_two_review.json`、結果は
`reports/CLASSIFICATION_QUADRIC_CONE_RANK_TWO_REVIEW_20261107.json`。
事前 gate PASS 後の実行は 6.470 秒、最大 RSS 23,672 KiB である。

## 2. vertex を含む全 tangent 源

座標を `(u0,u1,u2,v,h)`、`q=u0u2-u1^2`、`S=(h=q=0)` とする。
S の Poisson invariance により、raw contraction は S 上で
`Pi#dh=lambda E`、`Pi#dq=mu E` となる。twist はそれぞれ OS(1)、OS(2) である。
vertex では dq が零で E は零でないので mu(vertex)=0。従って mu の quadratic
lift に v² 項はなく、q の三 gradient が u0,u1,u2 を張ることから、linear vector
field W の W(q) として全て表せる。W(h) は独立に指定でき、両 contraction を
同時に Euler gauge で消せる。vertex で割る操作は用いない。

既済の全重複度 abelian-conormal 入口により Sreg 上の bracket は零である。
S は integral CI、I/I² は vertex でも局所自由なので全 S で零となる。
global tangent lift の strict bracket は `{h,q} in I²` を満たす。
CI prime ideal の square の primary 性と飽和性により、これは global homogeneous
係数条件であり、vertex 上だけの自由な例外を残さない。

S は normal Gorenstein、omegaS=OS(-2) である。従って restriction は全
`H0(OS(2))`、次元九の F で表される。四変数の Jacobian bivector J(q,F) は
h,q を strict Casimir とし、Euler formula による restriction は非零定数倍の F
である。これを引いた ambient section は全 S 上零なので、既済の Koszul/Euler
全源から

```
Pi = J(q,F) + hL + qA
```

となる。L は全五十、A は全十、F mod q は全九係数である。Sreg で零の
ambient bundle section は integral S 上零なので、vertex 支持だけの追加 section
は無い。F=0 も上界には安全な大きい源として保持する。

全百五十 quadratic 係数に、cone parameterization `(s²,st,t²,v,0)` 上の
`E wedge Pi#dh=E wedge Pi#dq=0` を独立に課した matrix の rank は六十二。
raw 六十九方向は独立、全 Euler 二十五方向の rank は二十四、合併 rank は八十八。
従って quotient tangent 次元は `150-62-24=64=69-5` であり、表示した源は完全である。
全 conormal membership matrix 35×66 の rank 十七、完全 kernel 四十九、raw 源への
独立十一条件も再構成した。F 九方向を加えた全 conormal raw 次元は五十八である。

## 3. 全 actual normalizer と全三十三方向

h-row を `h ell_i+q a_i` と書くと `{h,q} in I²` から a は q の radical、
vertex 方向だけに残る。`a=A(h,v)=0` の全枝では H が Poisson invariant であり
実 D に H 因子があるので、既約性に反する。a 非零では
`u_i'=u_i-b_i h/a` とその逆、vertex scaling により a=1、A(ui,v)=0 へ直す。
q の変化は h linear と h² だけで全 hL に吸収され、F 九係数も可逆に動く。
反転した a の補集合は上の H 因子で処理済みである。

残る三変数 constant bivector C は SO(q) と実 scalar 正規化により、
`C=partial0 wedge partial2`、`partial0 wedge partial1`、零の全三型に尽くされる。
非零 q-dual norm と非零 isotropic の全二軌道、及び零であり、未処理 Gram chart
を捨てていない。**rank 二では C=0 を除かない。**

各型で全五十九 unknown L50+F9 の mod h original Jacobi と conormal 十一式を
独立に構成した。未知同士の quadratic 項は mod h で消える。合併 rank 二十一、
raw kernel 三十八、Euler 五方向 quotient 後三十三の完全 affine 源を得た。
保存された particular、全 kernel の独立性と全積零、全 free coordinate の inverse、
Euler gauge を検証した。さらに全百五十 tensor、全三百五十 projective Jacobi、
実 Euler-contracted D の全百二十六係数を全三十三係数の多項式として再構成した。
chosen subfamily 又は vertex 近傍だけの確認ではない。

## 4. 誤った Ore 入口を使わず、元式で Dpencil=0 を導く

全 source の chart h=1 での形は

```
P = B0(u)+v Dpencil+(Y(u)+vT(u)) wedge partial_v,
deg B0<=2,  Dpencil constant,  deg Y<=3,  deg T<=1.
```

全五 Jacobi block を独立に再構成した。

```
Jac(B0)+Y wedge Dpencil = 0,
[B0,Dpencil]+T wedge Dpencil = 0,
Lie_Y B0-Y wedge T = 0,
Lie_T B0+Lie_Y Dpencil = 0,
Lie_T Dpencil = 0.
```

全三百五十式との dense chart 恒等式も確認した。global Euler 関係により
この chart 上の零と global projective Jacobi の零は同値である。
初回 producer の v-independent B 仮定は全 mod h 必要 source では誤りであり、
sem 型の F=u0v 方向では Dpencil12=2 が残る。これは full Poisson 解の存在例では
ない。初回失敗 source/runtime と、後の正しい pencil 証明書を区別して保持する。

radial s と `F_v=ell0 u0+ell1 u1+ell2 u2` の全四係数だけが現れる元最高空間次数
rowspace を再構成した。その理想は sem 型 `(ell0,ell2)`、nil 型 `(ell0,ell1)`、
C=0 では零である。この全帰結を入れると、全三型で Dpencil は恒等零になる。
従って**入口の仮定でなく、元 Jacobi の帰結として** Ore の全三 block に戻れる。
sem/nil ではそれぞれ Fv=kappa u1/kappa u2、他の低次係数は全て残る。

## 5. 全 graded fiber と全複素共鳴

`B=B2+B1+B0, T=T1+T0, Y=-u q+Y2+Y1+Y0` とし、base を kappa,s とする。
残る全 source は weight 一十四、二十二、三三の二十九係数である。
元の全 Jacobi row は parameter weight で homogeneous であり、stage j の式は
stage j の未知に affine-linear、係数 matrix は base のみに依存する。
従って固定 base 上の全 Poisson fiber は三 matrix の nullity 和以下である。
低 weight の nonlinear compatibility や未使用 row は支持を小さくするだけで、
この必要条件による上界を増やさない。元の各 row と全 weight を照合した。

| 全型 | 三 operator | generic nullities | 全 base を戻した generic 像上界 |
|---|---|---|---:|
| C sem | 38×14、21×12、7×3 | 3,1,0 | 24 |
| C nil | 35×14、22×12、8×3 | 3,1,0 | 23 |
| C=0、Fv sem | 33×14、21×12、7×3 | 3,1,0 | 23 |
| C=0、Fv nil | 33×14、22×12、8×3 | 3,1,0 | 22 |

C=0,Fv=0 は全 actual Pfaffian の h⁰ 係数が恒等零で、実 H 因子となることを
全源上で直接確認した。C=0,Fv 非零は残る SO(q) と scalar により u1/u2 の全二軌道。
q-norm 非零/零を両方保持する。tensor scalar は normal h scaling で補償できるので、
全 F と hL を可逆に保ったまま Fv の大きさを一にできる。両 orbit 次元は三/二。

全 matrix の非零 maximal minor と、全行に対する完全 kernel を独立に検証した。
sem 型の選択 minor 積の異なる全因子は十二 line

```
2k+1,2k-1,k,k+1,k-s,k+s,k+s+1,k+s+3,
3k-s,k-s+1,k+s+2,3k-s+1
```

である。nil 型は七 line `2k-3,2k+3,k-1,s,4k-1,k,2k-1`。
二変数 QQ(k,s) 上で minor を再計算し、これらの全因子で割った残りが非零定数で
あることまで検証した。各全 line を QQ(tau) に引き戻し、そこでの全 minor 根を
正確に戻した。sem は七十二、nil は十三の root certificates、全て有理根である。
交点は重複も意図的に保持する。C=0 の全根は sem で s=-3,-1,1,3、nil で s=0。
各根でも全 kernel と非零 minor を確認し、generic 非零領域だけで終えていない。

base 次元 generic 二/line 一/point 零（C=0 では一/零）を加えると、唯一の例外
nil,k=s=0 を除く全像上界は二十四以下である。旧 graded report の
`all_rank2_quadric_cone_G24_closed=false` は、この唯一 coarse 二十五を残した正しい
未完了結果として不変保存し、次節で補う。

## 6. 唯一の全 top 五係数は元 compatibility で一つ落ちる

nil,k=s=0 での三 nullity は (5,3,0)。完全 top kernel を p0,...,p4 として、
元 weight 二の**全二十四式**の middle 十二係数 matrix を再構成した。
rank は九、全 middle kernel は三次元。全 quadratic 右辺も独立に再構成した。
保存済み pencil rows 13,59,62,75 の係数 7/8,-1/8,1/8,1/8 による組合せは
middle の全十二係数を消し、正確に

```
p1*(p0+p2/2)=0
```

を与える。全 top 五自由係数上の非零多項式なので top 支持は四次元以下。
middle fiber 三以下、bottom 零により全 normalized source は七次元以下。
free parameter は一つも反転せず、二つの線形因子の両零枝を保持している。
追加 jet、点列、native GB、平坦性仮定は不要である。

## 7. 正規化を二重に引かず、全 PGL 像へ戻す

Aut(S) は十二次元、tensor scalar を加えて十三。allowed A の sem/nil/zero
strata は七/六/四次元なので、A 固定の stabilizer 次元は六/七/九である。
C=0 でさらに Fv の全二軌道を固定すると六/七になる。全自由 frame と tensor
scalar を用いた map の fiber 作用は自由である。個別 tensor の stabilizer が零
である必要はなく、PGL frame から group 元が復元される。

重要な点として、個々の k,s の壁がこの stabilizer に不変である必要はない。
まず有限 cover **全体の最大次元**で全 normalized Poisson source X を押さえ、
その X が actual stabilizer に不変であることを用いて fiber 次元を引く。
sem 全源は d<=6、nil 全源は本 compatibility を含め d<=7。
C=0,Fv sem/nil 全源は d<=5/6 である。従って像上界はそれぞれ

```
sem:       24+6-6 = 24,
nil:       24+7-7 = 24,
C0,Fv sem: 24+5-6 = 23,
C0,Fv nil: 24+6-7 = 23.
```

ここで stabilizer は既に tensor scalar を含み、別に射影化の一をもう一度
引かない。同値に固定 S affine 全源へ A orbit を戻してから Aut(S) 十二と
tensor scalar 一を引いても同じ値になる。全 R、全共鳴、C=0、Fv=0、a=0 の
補集合を処理したので本 G 枝は閉じる。E 成分帰属又は全分類の結論は含まない。

## 8. 凍結入力と再現

入口は `inputs/quadric_cone_rank_two_manifest_20261107.json`、SHA256
`bc232e7fd705ff409f195fba55ae974b7ea8cc02364f4f2652ba89a3243b11bc`。
全二十二入力と全 transitive 独立ライブラリの hash を検証した。
四 producer report の executed source bytes の hash もそれぞれ検証し、
mutable current runner を import して既存結果に見せかけていない。

```
python scripts/classification.py gate research/classification/jobs/quadric_cone_rank_two_review.json
python code/classification_quadric_cone_rank_two_review.py
```

gate は有限契約の確認であり数学の証明そのものではない。本判定は全源幾何と
上記の全 Fraction/QQ 恒等式を合わせた独立受理である。
