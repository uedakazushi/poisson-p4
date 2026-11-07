# 最小冪零・長 block translation：全支持を三正常形へ有限化

2026-11-07。主課題 `E_reduced_nonnormal`、純冪零 cofactor 表の型10、
Jordan 分割(2,1,1)、translation top 長2。
実退化式を単純 H と既約四次 F の積とする全支持を扱う。
この型の全点は、以下の三つの実 Poisson 代表の affine/PGL 軌道の和へ入る。
三代表の全点を既存 MPS の閉和へ被覆した。C は以下の有限厳密 slice と実 Poisson 曲線を使う。
A,B の全曲線・独立確認は `HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_SEMISIMPLE_DEFORMATIONS.md`。
三代表を三つの新成分と数えてはならない。

## 1. 元の全二次係数との同値

cofactor を `Z=y∂x+∂y` に固定する。補助の多項式座標

    r=x−y²/2, s=y

では Z=∂s、体積形式も変わらない。
この補助座標自体を projective 同値として使うことはしない。
Z に不変な P を

    P=∂s∧V+π,   V=V(r,z,w), π=π(r,z,w)

と書く。`div P=Z` は `div V=−1, div π=0` を与える。
三変数の divergence-free bivector は多項式 potential の Jacobian、`π=Jac Φ` と書ける。
元の x,y,z,w における全次数二以下を保つ条件を取ると、正確に

    Φ=δr²/2+rL+G(z,w),   L=η+βw−γz,
    V=(−αr−A₂(z,w), d(z,w), e(z,w)),
    d,e は affine,   d_z+e_w=α−1,
    deg A₂≤2, deg G≤3

となる。G の定数は π に寄与しないので零にできる。
例えば `Pyz=d` に r 項や二次の z,w 項があれば、元の `Pxz=πrz+y d` に
相殺不能の y の奇数次数三以上が現れる。従って d は affine で、e も同様である。
div π=0 から `πzw=δr+η+βw−γz`、
`πrz=βr+G_w, πrw=γr−G_z` が従う。

元座標での全六括弧は

    Pxy=αr+A₂,
    Pxz=βr+G_w+y d,
    Pxw=γr−G_z+y e,
    Pyz=d, Pyw=e, Pzw=δr+L.

これは単なる一切片ではない。下記の有限証明書で、既存の全90係数源の全25自由変数と
この表示の25変数の間の明示的な可逆 affine map を復元し、90係数全てを一致させた。
元の定数項・一次項も全て保持している。

## 2. 全 Jacobi と実 Pfaffian

`Ψ=VΦ+Φ` とすると、Ga 座標(r,s,z,w)での全四 Jacobi 式は順に

    −Ψ_w, Ψ_z, 0, Ψ_r.

従って全 Jacobi は Ψ が定数 C であることと同値である。
元座標の実 Pfaffian は `F=−VΦ=Φ−C`。
δ=0 なら最高四次 Pfaffian は恒等零となり、H は重複因子になる。
対象は H 単純なので δ≠0 が必要である。

Ψ の r² 係数は `δ(1/2−α)`。よって α=1/2。
r 係数は

    L/2+Vbase(L)−δA₂

であり、`A₂=(L/2+Vbase L)/δ` は affine になる。
ここで Vbase=d∂z+e∂w、`tr(linear Vbase)=−1/2`。
完成平方を

    R=r+L/δ, H₀=G−L²/(2δ)

とすると、恒等的に

    Φ=δR²/2+H₀, V(R)=−R/2,
    VΦ+Φ=Vbase(H₀)+H₀=C

である。`H=H₀−C` と置けば、実 Pfaffian と残る全条件は

    F=δR²/2+H, Vbase(H)=−H, deg H≤3.                 (1)

potential の定数 C を差し引く操作は tensor を変えない。
R への変更は、元座標では `x'=x+L(z,w)/δ` という実際の affine 変更である。
従って非線形補助座標の同型を PGL 同値へすり替えていない。

## 3. trace 条件付き affine 場の全 weighted 分岐

Vbase の二つの線形固有値を λ,μ とし、λ+μ=−1/2。
H=0 は F が平方となるので除かれる。
固有値が重複すると両方−1/4であり、次数三以下の多項式表現の半単純固有値は
`0,−1/4,−1/2,−3/4` のみ。Jordan の冪零部分があっても固有値−1は現れず、H=0となる。

相異なる固有値では、非零固有値方向の定数場を affine 平行移動で除く。
零固有値方向の定数場 κ∂w は捨てずに残す。
この場合も、有限次元の次数三以下の表現上で半単純部分と可換冪零部分を分ければ、
H は weight−1 部分にあり、可換冪零部分にも消される。

`iλ+jμ=−1, 1≤i+j≤3` を全て列挙すると、座標交換を除いて次の三型に尽くされる。

| (λ,μ) | weight−1 の全次数三以下単項式 | H の全形 |
|---|---|---|
| (−1,1/2) | z | bz |
| (−1/3,−1/6) | z³ | bz³ |
| (−1/2,0) | z²,z²w | z²(a+bw) |

最初の二型では既約性から b≠0。
最後は `κ∂wH=κb z²` だから κb=0。
b=0 なら `δR²/2+az²` は複素数体上で二つの一次式の積になる。
既約性より b≠0、従って κ=0。零 weight の w を平行移動して a を消せる。
これで全 affine translation とその閉補集合も処理した。

残る `F=δR²/2+bz^i w^j` の非零係数 δ,b は、
z'=tz,w'=uw のみで正規化できる。
実 tensor の変換では

    δ'=δtu, b'=b t^(1−i)u^(1−j).

H=bz なら `(t,u)=(b/δ,1/b)`、H=bz³ なら `t³=b/δ,u=1/(δt)`、
H=bz²w なら `(t,u)=(b,1/(δb))` とすればよい。
複素数体では必要な非零 cube root が存在し、その有限代数選択により源を尽くす。
この変更は全て invertible affine/PGL 変更であり、追加の有限 jet は使わない。

## 4. 三つの全 Poisson 代表

元の x,y,z,w に戻し、`r=x−y²/2` と書く。
順序(xy,xz,xw,yz,yw,zw)の代表は

    A=(r/2, −yz, −1+yw/2, −z, w/2, r),
    B=(r/2, −yz/3, −3z²−yw/6, −z/3, −w/6, r),
    C=(r/2, z²−yz/2, −2zw, −z/2, 0, r).

全ての括弧は元座標の二次以下で、全四 strict Jacobi が零、
`div P=(y,1,0,0)`。
実 affine Pfaffian はそれぞれ

    r²/2+z, r²/2+z³, r²/2+z²w,

最高四次 Pfaffian は全て `y⁴/8` である。
三つの affine F は、順に変数 z に関する原始一次式、二乗と三乗の互いに素な cusp 多項式、
変数 w に関する原始一次式であり、全て既約である。
元座標への戻しは多項式自己同型なので既約性を保ち、最高次数四が非零だから斉次化に h 因子もない。

## 5. 証明書と未解決部分

実装は `code/classification_translated_minimal_quartic.py`、有限契約は
`jobs/hyperplane_quartic_translated_minimal_normal_forms.json`、証明書は
`reports/CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_NORMAL_FORMS_20261107.json`。
producer・CAS を import せず、既存独立 reviewer の Fraction と疎多項式演算を再利用した。
gate通過後、実測0.0675秒、最大RSS 14,080 KiBで以下を確認した。

- 全25変数と既存全90係数源の明示的な可逆 affine map と90係数の完全一致。
- 全次数制約、divergence、cofactor 不変性、`Pf=−VΦ`。
- 全 Jacobi と Ψ の全偏微分の一致、r²・r 係数、完成平方と V(R) の恒等式。
- 全次数三以下 weight の完全列挙、零固有値 translation の残条件。
- 三つの実代表の全 Jacobi・divergence・実 Pf・最高部 Pf。

確認したのは、この cofactor 型の対象全支持を高々三つの実 affine/PGL 軌道へ送ることである。
A,B の全点は別稿の実半単純 cofactor 曲線で既知 MPS へ被覆した。
C の全点帰属は第8–10節で閉じる。実因子の形だけから MPS 帰属を推測していない。

## 6. 全周囲の接空間・実軌道と S61 graph

同じ reusable runner の `--ambient` は、各代表について元の全150係数、
全350 projective Jacobi 式の微分 D と、実 GL5・Euler 鉛直方向・scalar の行列 G を生成する。
別作者の sparse350 engine を再利用し、QQ の完全 rank・非零 minor・完全 kernel を保存した。
鉛直像の rank は全て24で、`DG=0` も全項確認した。

| 代表 | rank D | rank G | 射影接次元 | 実 PGL 軌道次元 | ker D / G |
|---|---:|---:|---:|---:|---:|
| A linear | 97 | 47 | 28 | 22 | 6 |
| B cubic | 100 | 47 | 25 | 22 | 3 |
| C mixed | 101 | 47 | 24 | 22 | 2 |

いずれも実軌道だけが接空間を埋める smooth 判定にはならない。
明示的 rank 四基点 A,B では(0,0,1,0)、C では(0,0,1,1)を取り、
元 tensor の Pfaffian が1であることを使った。
その点を原点へ平行移動し、厳密 Darboux frame と perspective 除去を行った。
全150係数の quotient 再構成、全61変数、元の全134式は全て一致・零となった。
T39 と T40 の全45 polynomial graph 式は、三代表の全てで非零残差を持った。
全残差を保存し、どれもその二 graph に所属するとは主張しない。

T01 は退化式が五一次因子の積となる閉条件に入るが、対象には既約四次因子があるため該当しない。
この三成分の排除を、他の既知成分全てへの非所属へ拡張してはならない。
各22次元軌道がより大きい成分の境界に含まれる可能性は残っている。

計算は `jobs/hyperplane_quartic_translated_minimal_ambient.json` の gate 通過後に実行した。
全証明書は `CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_AMBIENT_20261107.json`、
実測10.02秒。新しい ambient 行列自体の独立作者による再生成は別の検証段階である。

## 7. 真の法線商上の二次初期支持

`--quadratic` では上の全 D と実 G の matrix hash を再現し、
全 kernel を G とともに張る2,3,6本の独立法線を用いた。
従って固定 H や固定 cofactor の切片だけの法線ではない。
全法線二単項式の Jacobi 係数 Q を生成し、D の非零最大 minor を使って
`Q mod im D` を厳密に計算した。
各独立な余核行について、全150列の D を消す左汎関数も保存した。

法線の順序は QUADRATIC 証明書の `normal_bivectors` と整数 scale に従う。
C の全二次初期 ideal は

    (t0*t1, t1²)

で、根基は(t1)、t0 直線が残る。
B は

    (t1*t2, t2²)

で、根基は(t2)、t0,t1 平面が残る。
A では全十三生成元を得た。

    t0²+8t1t3, t0t1, t0t2+5t3t4/6, t0t3, t0t4,
    t0t5−75t2t4/188, t1t4, t1t5−5t4²/4,
    t2t3, t2t5, t3², t4t5, t5².

その根基は(t0,t3,t4,t5)で、t1,t2 平面が残る。
実際、t3²,t5²から t3=t5=0、続いて t4²=0,t0²=0。
逆に四変数を零とすると全十三式は零である。

三者とも二次初期支持は非孤立であり、これだけから局所成分次元22を結論できない。
同時に、この初期 cone の全点が実 Poisson 変形へ持ち上がるとも主張しない。
C は次節の有限 equivariant 支持構成から実 Poisson 曲線を得た。
A,B は別稿の実半単純 cofactor 曲線で閉じた。次数を選び直す無限の jet 反復は行っていない。

finite job は `jobs/hyperplane_quartic_translated_minimal_quadratic.json`、
証明書は `CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_QUADRATIC_20261107.json`。
全三型を一回の有限実行で確認し、実測6.16秒だった。

## 8. C の全 ambient における有限厳密 slice

C の整数代表を4倍に取り、座標(h,x,y,z,w)の重みを(0,2,1,1,2)とする。
P の character は1なので、`t⁻¹ diag(1,t²,t,t,t²)_*` は P を固定する。
全150係数の相対重みは

    wi+wj−sum(e_k*w_k)−1

で、範囲は[−4,3]。実 G の全列は一つの重みのみに属し、G は重み安定である。
G の独立47列と47座標行で非零 minor を選ぶ。
その47座標を零に固定した補う103座標空間 S は、全 C150 で G の真の重み安定補空間である。
kernel 内だけの補空間を真の slice と仮定していない。

S 内の101座標と全350式のうち101式の微分に非零 minor がある。
残る二自由座標を u,v とすると相対重みは−1,−2。
元の全350式も各座標 monomial が重みを持つため、選んだ101式の空間も重み安定である。
形式陰関数定理の解は一意であり、従って C* 同変となる。
u,v の重みが負なので、相対重み−kの implicit 座標は weighted degree k の多項式しか持てない。
全 C150 の重み下限−4により、全 implicit 級数の weighted degree は先験的に4以下である。
これは法線だけの重みから全高次を推定する議論ではない。

その一意な有限多項式を全て構成し、選択101式が恒等零、
group cut の47座標も恒等的に固定、二自由座標が正確にu,vであることを確認した。
全350式を文字通り代入した残余 ideal は正確に

    (u*v,v²).                                                   (2)

実際、残余行111は12v²、行158は3uv/2。
他の全残余 monomial は uv 又は v² で割れるため、ideal の両向き包含が明らかである。
したがって真の形式被約 slice は v=0 の直線であり、C の22次元軌道だけでは成分にならない。
この点を通る被約成分は一つで、射影次元23である。
正確な nil 指数の追求は、この支持・帰属の結論の前提にしない。

全座標 cut、全 minor、全補正、全350残余を
`CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_WEIGHTED_MIXED_20261107.json` に保存した。
有限契約は `jobs/hyperplane_quartic_translated_minimal_weighted_mixed.json`。
実測1.50秒であり、任意に jet 次数を増やした反復ではない。

## 9. 全 u に対する実 Poisson 曲線と半単純入口

上の厳密 slice の v=0 は全 u に対する実多項式 Poisson 曲線であり、全六括弧は

    Pxy=2hx+u²x²/16−y²,
    Pxz=−2yz+4z², Pxw=−8wz,
    Pyz=−2hz−u²xz/8−uz²,
    Pyw=2uwz,
    Pzw=4hx+u²x²/8+uwz−2y².

P の h 行は零である。u=0 は元の4倍Cで、全350式は全 u に対して恒等零。
トーラス作用は u を u/t へ送るため、u≠0 の全点は一つの実 PGL・scalar 軌道になる。
u=1 では全150方向と350式の QQ rank は `rank D=102,rank G=48`。
射影接次元と実軌道次元がともに23なので、唯一の23次元全空間成分の滑らかな点である。
この計算だけから新成分だとは結論しない。

実 divergence は

    Z_u=4y∂x+(4h+u²x/4)∂y−uz∂z+uw∂w.

u≠0、h=1 で modular zero は

    (x,y,z,w)=(−16/u²,0,0,0).

そこで `Pxy=−16/u²,Pzw=−32/u²`、他の四括弧は零。
定数 Pfaffian は512/u⁴であり、定数交代形式は非退化である。
ξ=x+16/u² とし、負固有座標と正固有座標を

    (a,b)=(uξ/4−y,z), (c,d)=(uξ/4+y,w)

と取る。線形部分の行列式は−u/2で、cofactor 固有値は(−u,−u,u,u)。
定数 crossing 行列は `diag(−8/u,−32/u²)` なので全 u≠0 で rank 二である。
tensor を u で割れば、これは既済の半単純型20、定数 crossing rank 二の全源定理の仮定を満たす。
個々の座標係数の標本を適用条件の代わりにはしていない。

`HYPERPLANE_QUARTIC_LARGE_SEMISIMPLE_SUPPORT.md` の全支持 MPS 定理により、
u≠0 の全点は既知40成分の有限閉和に入る。
この和の閉性と実曲線から u=0 も同じ閉和へ入る。
PGL と非零 scalar を戻せば、型10の C に帰着する全点も被覆される。
唯一の23次元成分は既知 MPS の一つであり、新しい成分数を追加しない。
T39/T40 の graph が不一致だったことと、この他の MPS への所属は矛盾しない。

generic 点の完全な全150/350・S61 graph 証明書は
`CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_GENERIC_MIXED_20261107.json`。
その全134式と両全45 graph の確認を、新たな全成分非所属と読み替えない。

## 10. 独立算術と適用依存

`code/classification_translated_minimal_quartic_review.py` は producer・producer engine・CAS を import しない。
別作者の nil22 independent reviewer の Fraction 疎演算だけを再利用し、
全 D/G・47 cut minor・101 implicit minor・全 C150/350 の重みを別に生成した。
全 bounded implicit 多項式と全350残余を再展開し、(2)を Groebner に依存せず検証した。
実 u 曲線の全 Jacobi・modular 場、全 u≠0 の modular zero・eigenframe・crossing rank、
generic 点の D102/G48 も独立に確認した。
実測0.488秒、最大RSS 18,432 KiB。

独立証明書は `CLASSIFICATION_HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_MIXED_REVIEW_20261107.json`、
契約は `jobs/hyperplane_quartic_translated_minimal_mixed_review.json`。
数学的完全性は本稿の全25係数への橋と重み安定 IFT の論理、
MPS 帰属は既済半単純型20の全支持定理とその外部 MPS 定理を依存として区別する。
別稿の A/B 被覆を合わせると、全型10の支持は既存 MPS へ入る。
全分類・他の冪零型は未完成である。
