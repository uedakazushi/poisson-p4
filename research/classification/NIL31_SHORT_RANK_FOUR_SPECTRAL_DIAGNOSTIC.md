# nil31短鎖 r=2/3：一つの rank四SS12 spectral frame 全体の障害

2026-11-08、PRIMARY E、`/root/quartic_review`。
marked 重複H・rank二半単純層の閉包を排除した後、
別の実幾何として rank四 SS12 companion frame を一回だけ検査した。
**この明示 frame 全体は target に達しない。全 r2/3 E所属は未解決。**
generic SS12 全像・別の spectral frame の閉包非所属とは主張しない。
追加の profile、epsilon次数、normal jet 又はGBは実行しない。

## 1. genuine route change と全 parent 源

target は canonical(15)の t=d1。
modular 場に次の rank四変形を選んだ。

    Nε=y∂x+z∂y+(-4ε³x+4ε²y−εz)∂z+(1+εw)∂w。

linear eigenvalues は(-2,-1,1,2)ε。
zero-eigen hyperplane h を短い Jordan(2) block 側に保つので、
pure sl2 deformation が zero hyperplane を zへ送る欠点を避ける候補だった。
しかし modular field の存在を Poisson curve の存在と取り違えない。

実 spectral frame は

    X=z−εy−2ε²x, Y=z−4ε²x, Z=h+εw,
    W=z+3εy+2ε²x, h=h、det=12ε⁴。

inverse は

    x=(3X+W−4Y)/(12ε²), y=(-3X+2Y+W)/(6ε),
    z=(3X−Y+W)/3, w=(Z−h)/ε。

全 source は公開90係数線形解の SS12 九parameter源を使う。
ε·Fε^-1_*P_SS(n) を旧座標の候補とした。
この九次元 linear space が全て Poisson なのではなく、必要十分な七真正関係

    n4(n7+n8−1), n2n4, n1n4−n3(n7+n8),
    n0(n7−1), n6(n5−n8), n2(2n5−1),
    n0n6+n1(1−n5)−n2n3

を全係数のまま保持する設計である。
今回の不整合はより弱い線形 endpoint 条件だけで得られたため、
これらの nonlinear support 計算、候補 actualcurve の350/Pf検証は実行する段階に達していない。

## 2. finite window と一回の正確な負の結果

Fε=M·diag(1,ε²,ε,1,ε)。h-Casimir quadratic tensors の relative weights は[-3,3]。
ε·inverse pullback は pole高々二。
regular oldcurve を ε^-1 forward した parent係数の pole は高々四。
parent tail ε³以上は old negativepowers/endpointに寄与しない。
従って n_i の Laurent[-4,2]、全63係数は endpoint の必要情報を尽くす。
七 nonlinear関係の高次までこの truncation が尽くすとは仮定せず、
literal Laurent curve を作れた場合に全関係を恒等検証する十分 ansatz として契約した。

gate PASS 後の一回の exact QQ 実行は0.138218秒、最大RSSは raw report参照。
全 negativepower と endpoint の行列は298×63。
全行列/RHS/RREF/constant-frame全10basisを公開 reportに保存した。
linear rank60、augmented rank61で不整合。
しかし結論に必要なのは次の一行ずつだけである。

* row53：Pxy の monomial hy、epsilon exponent0。全63係数が零、RHS −1。
* row91：Pxz の monomial hz、epsilon exponent0。全63係数が零、RHS +1。

各 unit row functional が coefficient matrix を消し、非零 RHS を残す。
chosen minor/prime の rank だけから全支持を推測したものではない。

## 3. 任意 Laurent/rational parent coefficients への手証明

上の二つの monomial は truncation の都合で消えたのではない。
元 SS12 全源の long-three indices X,Y,W 内の括弧は

    PXY=(n5+n7−1)XY,
    PXW=n1h²+n2YZ+(2−n7−n8)WX,
    PYW=n6hZ+n7WY。

他の括弧は Z-index を含む。inverse の ∂Z=ε^-1∂w なので、
それらは旧 xy/xz coefficients を作れない。
long 内で h を持つのは n1h²,n2YZ,n6hZ だけ。
Y=z−4ε²x は旧 y を含まず、Z=h+εw。
従って旧 Pxy に monomial hy は存在しない。

また inverse tangent vectors は

    ∂X=(4ε²)^-1∂x−(2ε)^-1∂y+∂z,
    ∂W=(12ε²)^-1∂x+(6ε)^-1∂y+∂z/3。

その wedge の xz成分は1/(12ε²)−1/(12ε²)=0。
ゆえに n1h²/n2YZ は旧 Pxz の hz を作らず、n6hZ は h²/hwだけを持つ。
旧 Pxz に monomial hz も存在しない。
全 n_i を任意 meromorphic functions としてもこの欠落は同じで、
Euler correctionや endpoint gauge の問題でもない：frameはhを固定し、全源がh-Casimirだからである。

target は Pxy に −t hy、Pxz に +t hz を持ち、rank四条件 t!=0。
従って **この fixed spectral frame の任意 whole SS12 parent curve** が不可能。
どの Laurent tail や nonlinear parameter relation も欠けた二monomialを補えない。
これを別 frame/eigenratio にまで拡張しない。

## 4. 再現と未解決の範囲

研究用 wrapper は増やさず、公開 input JSON 内の exact program を一回実行した。
入力・全source・parameter products・frame/polebudget・resource/stop契約は
`jobs/nil31_short_rank_four_spectral_curve.json` に固定。
raw reportは `CLASSIFICATION_NIL31_SHORT_RANK_FOUR_SPECTRAL_CURVE_20261108.json`。
実行は成功した negative diagnostic であり、actualcurve成功ではない。
SymPyはexact線形/多項式表現用、上の最弱不整合証明は手でも確認できる。

全 shortnil r2/3 source の E所属、任意 MPS 閉和からの非所属、componenthood は未解決。
同じ companion の parameter次数・profile を変えて再試行しない。
次手にはこの monomial欠落を破る別の幾何的親の理由が必要である。
