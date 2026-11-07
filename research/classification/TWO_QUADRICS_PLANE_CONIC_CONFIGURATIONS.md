# 平面二conicと一直線：全七free配置と全六磁気源

2026-11-07、課題E_reduced_nonnormal。
実因子 `D=H q₁q₂` の本質変数数が三の場合を扱う。
q₁,q₂は相異なる既約二次式なので、平面内では滑らかなconicである。
全支持を七配置（うち一つは一パラメータ族）と全六磁気係数へ帰着する。
一般の四・五本質変数を持つ配置の被覆は本ノートの対象外である。
七配置中、CPLUS族とそのcontact3境界は既存CPLUS_EXTへ入る。
さらに二固定配置の全閉包を23次元全成分と独立証明した。
contact4_secant_at_contactは全六磁気源をβ=0へ送って既知有限閉和へ覆った。
他の二固定配置の成分帰属は現時点では別義務である。

## 1. fiber依存を落とさずfree平面五次式へ帰着する

本質変数数は二以下ではあり得ない。二変数の二次式はC上で可約となるからである。
H=hとし、base変数をh,x,y、fiberをz,wとする。
h-Casimir代表Pはaffine h=1上で次数二以下のstrict Poissonである。

`J=∂xq₁ ∂yq₂−∂yq₁ ∂xq₂` と置く。
Jが非零ならq₁q₂と互いに素である。
実際、どちらかの既約qᵢがJを割れば、次数比較から `J=cqᵢ`。
他方qⱼの二変数Hamiltonian場はaffine線形かつ発散零であり、
h固定で三変数へ斉次化してもtrace零である。
これが非退化三変数二次形式qᵢを非零weight cでscaleすることは
`2tr X=3c` に反する。従ってc=0となり、J≠0に反する。

J=0なら、Hに対する接/非接の二正常形
`q₁=y−x²` 又は `q₁=xy−1` において直接係数比較すると、
二次以下のq₂は `q₂=aq₁+b` である。
相異なる既約性からa,b≠0。
normalなaffine線形場 `Xᵢ=P#dqᵢ/qᵢ` は
`q₂X₂=aq₁X₁` を満たす。
q₁,q₂の互いに素性と次数からX₁=X₂=0となり、非退化開上のdqᵢ≠0に反する。
したがってJ=0の分岐はrank四では空である。

q₁,q₂のnormal性より、両qᵢが `Pxy J` を割る。
Jが互いに素でPxyの次数は高々二なので、**Pxy=0** が従う。
従って全構造は

`P=U∧∂z+V∧∂w+B∂z∧∂w`、

`det(U,V)=非零定数×q₁q₂`                             (1)

となる。U,Vは初めはz,wに依存してよい。
z=w=0へ評価した二次anchorとEuler場Eの行列式は非零定数×hq₁q₂。
Saitoの判定から平面五次曲線hq₁q₂はfreeで、基底次数は1,2,2である。

このfree基底をC[h,x,y]からC[h,x,y,z,w]へ拡張する。
元の二次anchorをこの基底で書くと、次数二の二anchorの係数は定数であり、
Euler係数はh成分零から零となる。従って元のU,Vのfiber依存も消える。
一定のGL2 fiber変更を使えば、選んだfree基底のh-horizontal二anchorそのものへ戻る。
`[U,V]=aU+bV` のa,bはbase上の一次式であり、全Jacobiは

`B_z=−a`, `B_w=−b`

と同値である。よって**全六磁気係数を持つ全源**は

`P=U∧∂z+V∧∂w+(-az-bw+R₂(h,x,y))∂z∧∂w`、             (2)

`R₂=m₀h²+m₁hx+m₂hy+m₃x²+m₄xy+m₅y²`。
逆に(2)は全R₂でPoissonとなり、実因子はhq₁q₂のままである。
この還元は先にfiber係数を零とした小切片ではない。

free exponents(2,2)のHilbert–Burch分解から全Tjurina数は

`τ=4²−2·2=12`。                                       (3)

以下ではこの必要条件を使い、最後に完全Saito核で十分性も確認する。

## 2. contact budgetと全候補

sを三成分間の相異なる交点数、tを三成分の共通点数とする。
三成分は全て滑らかなので、局所Milnor積公式とBézout
`I(H,q₁)+I(H,q₂)+I(q₁,q₂)=2+2+4=8` から

`Σμ=16−s−t`。

(3)とτ≤μより `s+t≤4` が必要である。
Hは各conicに対して接線（一交点）又は割線（二交点）である。
次の四場合は全配置を尽くす。

* Hが両conicに同一点で接する。二conicのその点での交差数は2,3,4。
  2なら残りの交点分割は1+1又は2、3なら残りは1、4なら他の交点はない。
* Hが両conicに異なる点で接する。二conicの相異なる交点は高々二。
  接触分割4と2+2の二conicは、共通接線が共通接点でのものだけである。
  残る分割3+1には、異なる接点を持つ共通接線が一つある。
* Hが一方だけの接線。他方との二交点に接点が含まれないなら、
  二conicはH外の一点で交差数4を持つ。
  接点が含まれるなら、その点での二conicは横断し、H外で交差数3を持つ。
* Hが両方の割線。共通三重点が零ならbudgetに反する。
  三重点が一つなら、二conicの全交差数4がそこに集中する。
  三重点が二つなら、二conicの交点分割は2+2又は3+1である。

共通接線の上の主張は、`q=xz−y²` に対する次の正常形で直接分かる。
分割4は(q,q+x²)、3+1は(q,q+xy)、2+2は(q,xz−uy²)、u≠0,1。
共通接線のdual方程式を比べると、4ではH=xだけ、2+2ではH=x又はzだけ、
3+1ではH=xのほか `H=x+8y+16z` を得る。
これらの正常形は、qのP¹正規化上の交点divisorをPGL2で動かして得る。
残るqの倍数係数は、3+1と4では二点又は一点を固定するtorusで規格化でき、
2+2では一パラメータuとして全て残す。

同一点の共通接線H=xでは、q₁=qと置いた任意のq₂は

`q₂=a xz−y²+cxy+dx²`、a≠0

へscaleできる。a≠1ならqを保つshearでcを消せる。
残るdが非零なら1へscaleし、零なら分割2+2となる。
a=1,c≠0ならshearでdを消してcontact3、a=1,c=0,d≠0ならcontact4。
同一conicは除く。従ってこの場合の連続parameterも捨てていない。

これで全部で次の十候補となる。q=xz−y²、u≠0,1とする。

| ID | H | q₂ | free(2,2) |
|---|---|---|---|
| shared_tangent_generic | x | uxz−y²+x² | 可、全u |
| shared_tangent_contact22 | x | xz−uy² | 不可 |
| shared_tangent_contact3 | x | q+xy | 可 |
| shared_tangent_contact4 | x | q+x² | 不可 |
| contact31_common_bitangent | x+8y+16z | q+xy | 可 |
| contact4_other_tangent | z | q+x² | 可 |
| contact31_tangent_at_simple_intersection | z | q+xy | 可 |
| contact4_secant_at_contact | y | q+x² | 可 |
| contact22_chord | y | xz−uy² | 不可 |
| contact31_chord | y | q+xy | 可 |

## 3. 三つの候補の一様排除

順にshared_tangent_contact22、shared_tangent_contact4、contact22_chordでは
次の非零な**一次gradient syzygy**が全parameterで存在する。

`(4x,−y,−6z)`, `(0,x,2y)`, `(x,0,−z)`。

いずれも `δ(hq₁q₂)=0` を直接満たす。
free次数1,2,2のmoduleでは、次数一の場はEulerの定数倍だけである。
そのgradient syzygyは零なので矛盾する。
従ってこれらの全源は空であり、特別なuやその閉補集合も残らない。
余分なgeneric kernelのpivot因子を追わずに、この一様な弱い障害で全配置を排除する。

## 4. 全七配置の完全Saito核と全六磁気源

`CLASSIFICATION_TWO_QUADRICS_PLANE_CONICS_20261107.json` は、上表の九つについて
全18二次gradient係数を取り、全線形行列・完全kernel・非零rank minorと
全Saito行列式 `det(E,δᵢ,δⱼ)/(Hq₁q₂)` を保存する。
十番目は第3節のglobal一次syzygyで一様排除している。
七つの可の行は全てkernel次元二、Saito行列式は非零定数（parameterの有理関数）である。
全表示分母・minor・Saito係数の零点を保持した。
対象内の唯一の必要な例外はshared_tangent_genericのu=1/2で、
別の全kernel次元二と非零Saito行列式を保存する。
不可のcontact22_chordもu=−1で完全に再構成した。
u=0は可約conic、u=1は同一conic又は上表のcontact4へ既に接続されている。

独立Fraction(u) checkerは全行列、完全kernel、全rank minor、Saito行列式、
全例外と三つの一次syzygyを再構成した（0.419秒）。
報告は `CLASSIFICATION_TWO_QUADRICS_PLANE_CONICS_REVIEW_20261107.json`。

`--plane-sources` modeは各gradient場δを
`δ−(δH/H)E` によりh-horizontal化し、実際の座標変更でHをhへ送る。
得た全二anchorからa,bを解き、(2)の全六磁気係数を作った。
七配置と一つの有効例外、計八の完全源を
`CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json` に保存した。
全四strict Jacobiと実Pfaffianを全parameter上で確認した（1.158秒）。
別のFraction(u)実装は、座標逆行列、全kernel anchorの水平化、
全六磁気係数、全四Jacobiと実因数分解をproducerから独立に再構成した（0.222秒）。
報告は同名 `PLANE_SOURCES_REVIEW`。

## 5. CPLUS族とcontact3源の全支持は既存CPLUS_EXTへ入る

既存CPLUS_EXTの一般退化式は

`h{h²+2h(x−z)+2(x+z)²}{h²+2ch(x−z)+2(x+z)²}`。

一般cでは二conicはHに同一点で接し、そこで相互交差数2、H外で二つの横断交点を持つ。
従って本ノートのshared_tangent_generic配置族に属する。
既存 `legacy/v08/sections/cplus_component.tex` はこの実二parameter族の一般点で
被約supportが唯一の24次元成分CPLUS_EXTとなることを証明している。

第1節の全source還元により、その既知一般点は本ノートの全六磁気源のPGL飽和像に入る。
parameter uの既約開上のfree-kernel/frame bundleと全磁気A⁶は既約である。
その像閉包を含む全空間成分は、上の唯一成分点を通るのでCPLUS_EXTである。
従ってgeneric配置族の**全六磁気源**がCPLUS_EXTへ含まれる。
u=1/2の別chartもkernel rank一定の同じframe bundle内であり、垂直な未知源を付け加えていない。

contact3源については、H=x、q₁=q、

`q₂(t)=q+xy+t*xz`

を使う。t=0はcontact3、一般t≠0はshared_tangent_genericである。
t=0でkernel次元二かつSaito行列式非零、一般tでも同じなので、
零点を含む開上でkernelはrank二vector bundle、Saito frameはその開部分となる。
したがって任意のt=0のframeを、t=0で非零なpivot minorを使って正則に延ばせる。
第1節の全磁気式(2)はそのframeの関数として正則であり、
任意の六磁気係数を固定したまま全Poisson族を作る。
一般fiberが閉CPLUS_EXTに入り、その閉性からt=0の**全磁気源**も入る。
これは最高部族の稠密性だけで低次延長を外挿する議論ではなく、
完全kernel bundleと明示的全磁気式による実Poisson族の存在を用いる。

## 6. 残る具体的義務

後続 `TWO_QUADRICS_PLANE_COMPONENTS.md` は、全150係数の真正横断支持と
独立Fraction再構成により、contact31_common_bitangentとcontact4_other_tangentの
全source閉包を相異なる23次元全成分と証明する。
前者は二次、後者は固定された真正sliceの三次初項で支持原点を得る。

`TWO_QUADRICS_PLANE_EVEN_SUPPORT.md` はcontact4_secant_at_contactの全六磁気源を
parameterの反転なしの座標変更で正確なβ=0代表へ送り、既知閉和へ被覆する。

未帰属の平面源は次の二つで、それぞれ全六磁気係数を保持している。

`contact31_tangent_at_simple_intersection`, `contact31_chord`。

実sourceの全周囲Jacobi・軌道・parameter接方向、真正横断障害又は既知成分への実変形を
使って各閉包を全空間成分として認識又は既知成分へ被覆する義務が残る。
この二源は既に二次・一回の必要な三次支持を検査し非孤立だったため、
次は高次jetではなく構造的全source変形へ切り替える。
一つのsource点のPfaffian因子型や低い軌道rankだけでは成分性を主張しない。
非平面な1+2+2の全rank配置、既約非正規quintic、他の未完了分類義務は別に残る。
