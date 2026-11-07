# 対数閉包の残る二極型 x²yz・x²q：全源の独立再構成と非優越

2026-11-07、課題G。対象・入口は `UNIMODULAR_LOGARITHMIC_SUPPORT.md` 第5–8節。
本ノートは、既済の実L112、x³y、x²y²を再計算せず、最後の二極型を閉じる。
数学的還元、QQ(u)上の有限証明書、その独立算術再構成を区別する。
一般五点の尽くしも、全成分数の上界89も前提にしない。

## 1. 全低次延長を失わない計算の意味

体積形式で双対なbivectorを `P=C+L+Q` とする。
元の全134式に対応する正規化されたsymplectic延長は、必ず
`div P=0` と `Pf(P)=非零定数` を満たす。
以下では `div L=0`、全三次Pfaffian係数、全二次・一次Pfaffian係数を残す。
したがって元のJacobi条件を満たす延長は必ずこの源へ入り、
余分な解を許すとしても空性・非優越の結論は弱まる向きにしか働かない。

全24個の一次係数を最初から自由に取り、全行列と完全kernelを保存する。
translation `L -> L+dQ(t)` の全imageと非零minorを計算し、
そのminorが保証する真の線形sliceのみを使う。
元のPfaffianが非零定数なら、translation後の定数値も常に非退化なので、
この操作でsymplectic点を失わない。
その後も六個の定数係数を全て取り、二次Pfaffianから可逆に解ける分だけを解く。
残る全係数と定数Pfaffianを公表する。

各表のrankは標本でのrankではない。係数体QQ(u)での完全kernelと
同じ大きさの非零minorによる上・下両側のrank証明である。
全表示分母とminorの零点は捨てず、別の完全源として再構成する。
元134式への橋は `C0_S61_SLICE_SCALAR_AUDIT_20261106.md` に従う。
本ノートから局所nilpotent構造の一致や元134式の新たな既約分解は主張しない。

## 2. x²yz：独立配置の全パラメータでsymplectic延長なし

部分分数補題より、分母 `D=x²yz` の閉有理一形式は

`λ dx/x + μ dy/y + ν dz/z + d(f₁/x)`, `λ+μ+ν=0`

と書ける。x,y,z,f₁が従属なら既済の三変数pullbackである。
独立ならf₁=wとでき、最高部は

`Q=(X+aN)∧(Y+bN)`,

`X=x∂x-2y∂y+w∂w`, `Y=y∂y-z∂z`, `N=x∂w`

となる。既存ノートの `a=2μ-λ`, `b=ν-μ` はresidue平面から全(a,b)平面への
可逆線形写像である。X,Y,Nは可換かつtrace零。
wのscaleにより、b≠0ならb=1,a任意、b=0,a≠0ならa=1にできる。
a=b=0は既済の五辺型であり、この補集合を捨てていない。

独立checkerはこの二つのvector fieldのwedgeからQを作り直し、
producerの六成分との一致、`div Q=0`、`Pf Q=0` を確認した。

| 全源 | 一次rank・pivot | translation rank | 正常L | 定数rank・minor | 定数源C |
|---|---|---:|---|---|---|
| b=1,a任意 | 20, 8 | 4 | 0 | 6, −2 | 0 |
| b=0,a=1 | 18, −4 | 4 | `(n₀x+n₁w)∂y∧∂z` | 5, 2 | `a₀∂y∧∂z` |

全表示分母・minorはaに関して非零定数であり、例外値はない。
両源で全Pfaffianが恒等零となる。従ってこの独立配置の全型にはsymplectic延長がない。
線形形式が従属する補集合は既済pullbackの非優越に接続し、その全配置の空性は主張しない。

算術入力は `CLASSIFICATION_UNIMODULAR_LOGARITHMIC_DOUBLE_LINES_20261107.json`。
独立結果は同名の `DOUBLE_LINES_REVIEW` report。
checkerはproducerをimportせず、標準ライブラリFraction(u)と疎多項式だけを使う。
全行列・完全kernel・translation slice・六定数消去・全残係数を再構成した。

## 3. x²qの全quadricを九源へ還元する証明

qを既約二次とする。部分分数補題から

`ω=λ dx/x+μ dq/q+d(f₁/x)`, `λ+2μ=0`

を得る。f₁がxに比例する場合は既済の因子付き最高部であり、
独立ならf₁=yを座標に取る。
`θ=x²qω` を直接微分すると

`dθ=2q dx∧dy+(y+4μx) dx∧dq-x dy∧dq`。

μ=0は既済の実L112源である。
μ≠0ではyのscaleと全体の非零scaleで

`γ=2q dx∧dy+(y+x) dx∧dq-x dy∧dq`                 (1)

へ帰着する。非零scaleはsymplectic性も一次評価rankも保つ。

(1)はx,y上のJordan場 `V=x∂x+(y+x)∂y` を使って記述できる。
Vを保つ線形座標変更には

`x'=a x`, `y'=a y+b x`, `a≠0`

がある。以下のqの還元は、このcentralizer、z,wの自由な線形変更、
z,wへx,yを加えるshear、およびqとγの共通非零scaleだけを使う。
従って(1)の形を失わない。

### q|x=y=0 のrank二

z,wの変更と平方完成で `q=zw+F₂(x,y)` とできる。
`F₂=A x²+Bxy+C y²` と書く。

* C≠0なら、許されたshear `y -> y+(B/2C)x` によりxy項を消す。
  q全体をCで割りz,wを再scaleすると `q=zw+y²+u x²` となる。
  `u=(A-B²/(4C))/C` は任意の複素数として全て残す。
* C=0,B≠0なら、同じshearでx²項を消し、scaleして `q=zw+xy`。
* C=B=0,A≠0なら `q=zw+x²`。
* A=B=C=0なら `q=zw`。

最後は可約qだが、余分な源として含めても結論を強めない。
とくにu=0を他型へ吸収したとして捨てることはしない。

### q|x=y=0 のrank一

zの平方完成により `q=z²/2+w m(x,y)+F₂(x,y)`。
m=0なら三変数pullbackである。
m=αx+βy≠0について、β=0ならwのscaleでm=x。
β≠0ならJordan centralizerのyへのxのshearとwのscaleでm=y。
wへx,yの一次式を加えて、F₂のmで割れる部分を全て消す。
m=xなら残りはc y²、m=yならc x²。
c=0とc≠0を分け、後者をq,z,wの非零scaleでc=1へ送る。
これで次の四源を得る。

`z²/2+xw`, `z²/2+xw+y²`, `z²/2+yw`, `z²/2+yw+x²`。

### q|x=y=0 のrank零

`q=z m(x,y)+w n(x,y)+F₂(x,y)`。
m,nが従属なら、z,wを変えて三変数pullbackとなる。
独立ならz,wの線形変更でmixed部をxz+ywとし、
z,wへのx,yのshearで任意のF₂を消す。
従って残る唯一の源は `q=xz+yw`。

以上で、任意q、全restriction rank、全mixed係数、全binary quadraticを尽くす。
補集合のpullback、μ=0、従属f₁は既済の全源定理へ接続した。

## 4. x²qの全係数証明書と唯一の例外

既存runner `code/classification_unimodular_logarithmic.py --double-quadric` を拡張利用した。
旧reportは旧producerの全文bytesとhashを保持しており、書き換えていない。
新reportは `CLASSIFICATION_UNIMODULAR_LOGARITHMIC_DOUBLE_QUADRIC_20261107.json`。
契約は `jobs/unimodular_logarithmic_double_quadric.json`。
全九源と必要な一例外、計十源の計算が0.6823秒で完了した。

| q | 一次rank | 一次pivot | 定数minor | 結論 |
|---|---:|---|---|---|
| zw+y²+u x², u≠0 | 18 | −256u² | 4u | Pf恒等零 |
| zw+y², u=0の全源 | 18 | −64 | 2 | 同上 |
| zw+xy | 18 | −16 | 1 | 同上 |
| zw+x² | 18 | −256 | 4 | 同上 |
| zw | 18 | −16 | 2 | 同上 |
| z²/2+xw | 18 | 54 | 3 | 同上 |
| z²/2+xw+y² | 18 | −72 | 3 | 同上 |
| z²/2+yw | 18 | 12 | 1 | 同上 |
| z²/2+yw+x² | 18 | 192 | 4 | 同上 |
| xz+yw | 16 | −180 | −3 | 全primitive一次評価rank≤3 |

全translation rankは4。
最初の九つの表示では正常一次自由度は2、定数自由度は1であり、
producerの全Pfaffianを文字通り再展開すると恒等零になる。
全minor・全kernel分母・定数消去分母の唯一の零点は、最初のfamilyのu=0である。
上表の独立した全源がこれを保持する。他の例外値はない。

## 5. xz+yw源：全支持で一次評価rank≤3

順序(12,13,14,23,24,34)で最高部は

`Q=(0,xy,-x²,xy+y²,-x²-xy,wx+3wy+3xz)`。

正常一次自由変数はn₀,n₁,n₂,n₃、定数自由変数はa₀,a₁。
全残Pfaffianの中に実際に

`-2n₀²/5`, `-(2n₀+5n₁)(12n₀+25n₁)/50`

がある。被約支持では第一式からn₀=0、第二式からn₁=0。
非退化性や他の残式を捨てた仮定は使わず、実際の全残式の帰結である。
この代入後の全正常一次部は

`L=(0,0,-n₂x/2,n₃y,-n₂y/2-(n₂/2+n₃)x,n₂z+n₃w)`。

従ってz,wの二つの評価列は同じ一方向∂z∧∂wへ入る。
Qの最初の五成分もx,yだけに依存するため、
戻す全translation微分dQ(t)にも同じ性質がある。
よって**元の全一次評価** `L:V -> Λ²V` のrankは高々3。
定数自由変数や未消去の残式がどの値であっても、この上界は成立する。

体積双対とDarboux線形座標変更は一次評価rankを保つ。
symplectic延長ではPfaffianの一次係数零により `C∧L=0` なので、
対応する一次二形式βは既にCに関してprimitiveである。
したがってこのrank上界はprimitive βの一次評価にも移る。
既存の一般底W16にはinjectiveな評価の例がある
（`REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` 第14節）。
rank≤3は真の閉条件であり、この源の全像も像閉包もW16を優越しない。
これは一点の微分rankから像次元の上界を推測した議論ではない。

## 6. 独立確認の範囲と残る全分類義務

`code/classification_unimodular_logarithmic_pole_review.py --double-quadric` は、
producer/CASをimportせず、独立L112 checkerのFraction(u)・疎多項式演算だけを再利用する。
上記九quadricのリストをchecker内に別途固定し、(1)からQを再構成する。
全十源について、以下を文字通り確認した（1.5787秒）。

* 全一次・定数係数行列、完全kernel、rank minor。
* 全translation image、真slice、sliceとtranslationが全kernelを張ること。
* 六定数を保持した可逆消去、全二次・一次・定数Pfaffian係数。
* 全minor・表示分母の零点がu=0だけで、その独立源が存在すること。
* rank零源の実際の二平方式、n₀=n₁=0後の全translationを含むrank≤3。

結果は `CLASSIFICATION_UNIMODULAR_LOGARITHMIC_DOUBLE_QUADRIC_REVIEW_20261107.json`。
default実行はx²yzの二源を同じcheckerで検証する（0.2959秒）。
両実行は専用job契約を `scripts/classification.py gate` で通過後に行った。
幾何還元・supportからの帰結は本ノートの数学的証明であり、
算術checkerが外部の全foliation分類や部分分数補題自体を証明したとはしない。

既存のx⁴、x³y、x²y²と本ノートのx²yz、x²qを合わせると、
第5節の閉integrating-factor incidenceが与える**全対数閉包**は、
全低次延長を含めてGを優越しない。
これは対数型からの新たな成分が一般底を支配しないという結論である。
残るE3の実Borel軌道閉包、特殊底の被覆、最終全成分尽くしは別義務であり、
全分類や上界89の証明完了は主張しない。
