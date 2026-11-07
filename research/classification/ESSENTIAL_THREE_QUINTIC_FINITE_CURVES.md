# 本質三変数の既約退化五次：全八曲線候補への有限還元

2026-11-07、PRIMARY E。全canonical Saito定理と全磁気除去定理を入力とする。
その帰結として対象の平面五次Fは有理cuspidal、全μ=τ=12、free exponents(2,2)である。
本稿は**必要な曲線候補を八つに尽くす**。全八候補のSaito核・Poisson構造の成分帰属は別義務であり、
八個の新成分を主張しない。

Nambaの原分類本文を未読のまま尽くしの証明にせず、以下のelementary linear-series、
Cremona、dual curveの議論で全parameterと補集合を処理する。
FZの原文と誤記の確認範囲は `ESSENTIAL_THREE_QUINTIC_CURVE_FRONTIER.md` に区別してある。

## 1. 重複度一〜四と、既に閉じた全linear-series枝

既約五次の尖点重複度は四以下で、全τ=12から少なくとも一つ特異点がある。
重複度五なら五次式は二変数のbinary形式になり可約となる。

最大重複度四の全族は `ESSENTIAL_THREE_QUINTIC_CURVE_FRONTIER.md` の直接係数比較で
`y⁴z+x⁵+C x³y²+D x²y³` へ還元される。
全二次gradient核はC≠0で次元一、C=D=0で一次対称性を持つ次元三、
C=0,D≠0でのみSaito(2,2)となる。
全D≠0は一つの型 `x⁵+x²y³+y⁴z` へ移る。

最大重複度三の全族は `ESSENTIAL_THREE_QUINTIC_MULTIPLICITY_THREE.md` により
E8+A4、E6+A6、E6+A4+A2の三型へ入る。
monomial E8+A4型は一次対称性からfreeでなく、他の三型はそれぞれ一つのprojective候補である。
全P1 sectionを保持した証明なので、残る連続parameterは無い。

最大重複度二の一尖点は `ESSENTIAL_THREE_QUINTIC_UNICUSPIDAL_A12.md` によりA12一型、
二尖点は `ESSENTIAL_THREE_QUINTIC_MULTIPLICITY_TWO_BICUSPIDAL.md` によりA8+A4一型へ入る。
これらの証明は結線の接触の補集合、basepoint枝、全有限Puiseux条件を保持している。

従って残るのは、最大重複度二で尖点が三つ以上の場合だけである。

## 2. 最大重複度二では尖点数は四以下

尖点pからのprojectionは正規化P1からP1へのdegree三のmapであり、
Riemann–Hurwitzの全ramification次数は四である。
p以外の各尖点は少なくとも一のramificationを与える。
従って尖点数は五以下となる。
五尖点なら全δ=6から一つはA4であり、それをpに取る。
その接線接触は四又は五なので、projectionはp自身でもramificationを持つ。
残る四尖点と合わせて少なくとも五となり矛盾する。
よって尖点数は四以下である。

全δ=6を二重度二の尖点へ分配すると、三尖点・四尖点の候補は

    A8+2A2,  A6+A4+A2,  3A4,
    A6+3A2,  2A4+2A2.

以下で前三候補の二つと、最後の候補を排除する。

## 3. A8+2A2は明示Cremona mapで排除する

A8尖点をp=[0:0:1]、接線をx=0とする。
local parameterでy=r²、x=r⁴+a6 r⁶+a8 r⁸+a9 r⁹+…、a9≠0とできる。
接触四のleading係数はglobal linear scaleで1に取った。
接線の残交点tはsmoothで、交点重複度は一である。

明示quadratic Cremona involution

    ψ1(x:y:z)=(x²:xy:y²−xz)

を取る。base pointはpだけで、pのnormalization上で三sectionに共通するorderは六である。
元degree5 sectionを代入したdegree10から六を引くので、imageはquarticである。
ψ1はx≠0でbirational isomorphism、他の二つのA2 cuspはその開にあるので保たれる。

pのimageをchart Y=xyで見ると

    X/Y=x/y=r²+…,
    Z/Y=(y²−x)/(xy)=−a6+(a6²−a8)r²−a9 r³+… .

従ってimageはA2 cuspとなる。
接線x=0上の残smooth点tのimageはt'=[0:0:1]で、image curveはt'でsmooth、
接線はX=0である。この線はpのimage cuspも通る。
image quarticは三つのA2 cuspを持ち、δ=3が算術種数を使い切る。

pのimage cuspからprojectionするとdegree二のP1 mapになる。
他の二cuspsで二つのramificationを持ち、さらにsmooth t'でその接線によりramificationを持つ。
Riemann–Hurwitzの総数二に反する。従ってA8+2A2は存在しない。

## 4. A6+A4+A2も明示Cremona mapで排除する

A6 cuspをp=[0:0:1]、接線をx=0、A4 cuspをq=[1:0:0]とする。
pqはy=0であり、p,qでともに接線でない。
pqの残交点sはsmoothで、p,qとは異なる。

次のquadratic birational mapを使う。

    ψ2(x:y:z)=(xy:y²:xz).

normalization上の共通factorはpでorder四、qでorder二、他では零なので、imageはdegree四となる。
pでy=r²、x=a4 r⁴+a6 r⁶+a7 r⁷+…、a4a7≠0とすれば、chart Y=y²で

    X/Y=x/y=a4 r²+…,
    Z/Y=x/y²=a4+a6 r²+a7 r³+…

となり、imageはA2 cuspである。
qでもgeneric transverse座標をy=r²と取れば、z/yはconstant＋even r²＋nonzero odd r³となり、
A4 cuspはA2 cuspへ移る。
第三のA2 cuspはexceptional二直線x=0,y=0に乗らないためそのまま保たれる。
従ってimageは三尖点quarticである。

smooth残交点sはs'=[0:0:1]へ移り、imageでの接線はY=0である。
qのimageは[1:0:c]なのでこの接線上にある。
qのimage cuspからのdegree二projectionは、他の二cuspsとs'の三点でramifyし矛盾する。
ゆえにA6+A4+A2も存在しない。

この二つの議論は、FZのCremona記述をexplicit mapsと全exceptional交点で確認したものである。
未確認のgeneric平面の一例だけを排除した議論ではない。

## 5. 2A4+2A2はdual quarticのgenusで排除する

degree五・四尖点・全multiplicity二では、Plücker公式からdualの次数は

    5*4−Σ(μ_p+m_p−1)=20−(12+4)=4.

各A4 cuspのlocal expansionはx=r²、y=c4 r⁴+c5 r⁵+…、c5≠0である。
dualのslope/interceptから、c4≠0ならdual cuspもA4、δ=2、
c4=0ならdualは(3,5)でδ=4となる。
二つのA4からdualにδ少なくとも四を与える。
同じdual点に写る場合は枝間交点のδが加わるので下界は減らない。
quarticの算術種数三に反する。
これはFZ原文の「dual cubic」という誤記を用いない補修である。

## 6. 三A4の全degree5 sectionsから一型を得る

三cuspsの逆像をs=0,t=0,s=t、target imagesを三座標点にする。
各coordinate sectionは他の二cuspsでそれぞれorder二を持つので、全sourceは

    X=t²(s−t)²(A s+t),
    Y=s²(s−t)²(s+B t),
    Z=s²t²(C s+D t),   C+D≠0                     (4)

に入る。X(p)とY(q)の非零係数を1にnormalizeしただけで、C,Dはprojective pairとして残す。
各cuspのlocal二座標のorder二・三vectorのdeterminantが零になるA4必要条件は、順に

    BC+(2B−1)D=0,
    AD+(2A−1)C=0,
    AB+2A+2B+3=0.                                  (5)

(C,D)≠(0,0)なので、最初の二式の行列式は

    −3AB+2A+2B−1=0.

最後の式と合わせてAB=−1、A+B=−1となる。
従ってA,BはT²+T−1の二根であり、二組はp,qの交換で同値である。
D=0なら最初の式とB≠0からC=0となり不可能。
Dを1へscaleするとC=−A−2で、C+D=B≠0となる。

このsourceはbasepointを持たず非直線imageのdegreeは五である。
各三点にmultiplicity二のcuspがあり、(5)からδは少なくとも二ずつである。
全δ budget六から正確に三A4で、他の特異点は無い。
従って全三A4 quinticは一つのprojective型で、(4),(5)はその具体representativeを与える。
有理係数の別parametrizationを使う場合も、その型との一致を次の有限算術で確認する。

## 7. A6+3A2は唯一A6 quarticのdualである

dual次数は同じPlücker公式で四である。
A6 cuspは接線接触四で、local expansion x=r²、y=c4 r⁴+c6 r⁶+c7 r⁷+…、c4c7≠0を持つ。
dual slope/interceptのeven項を消すと最初のodd項はnonzero r⁷なので、dualにもA6 cuspがある。
quarticのδ budget三を使い切り、dualは一尖点A6 quarticである。

一尖点A6 quarticの一意性も全sectionで確認できる。
接線sectionをX=s⁴、他をY=s²(t²+b st)、Z=t⁴+d st³+f s³tに取る。
target shearで省略した項はX,Yの倍だけであり、全sourceを保持する。
A6の最初のodd五が消える条件はd=2bである。
source tをt+b s/2へ移し、target shearを戻すとb=d=0にでき、

    (s⁴:s²t²:t⁴+f' s³t)

となる。f'=0はconicへのdegree二mapでquarticでない。
f'≠0はdiagonal scaleで−1へ移せるので一つのPGL3型である。
source/targetの交換により同じcurveは

    (s⁴−st³:s²t²:t⁴)

とも書ける。その定義式はQ=(y²−xz)²−yz³である。
chart t=1のWronskianは−2(8s³+1)で三相異なるsmooth inflectionを持つ。
したがってdualはA6 cuspと三つのA2 cuspを持つ実degree五curveである。
char0のbidualityから、全A6+3A2 quinticはこの一つのprojective型となる。
その具体dual parametrizationは

    (−2st⁴ : t²(4s³−t³) : −s²(2s³+t³)).             (6)

これは∇Qへquartic parametrizationを代入した三項から共通factor t⁷を除いて得られる。

## 8. 全八候補と次の算術契約

以上と既済linear-series証明を合わせ、対象の全curveは次の八型へ入る。

| 型 | 具体degree5 representative又はparametrization |
|---|---|
| 一尖点(4,5) | F=x⁵+x²y³+y⁴z |
| A12 | (s⁴t:s²t³−s⁵:t⁵−2s³t²) |
| E8+A4 | (s⁵:s³t²:t⁵+st⁴) |
| E6+A6 | (s⁵+s⁴t:s³t²:t⁵−st⁴) |
| E6+A4+A2 | (s⁵+s⁴t:s³t²:t⁵+3st⁴) |
| A8+A4 | (s⁵+s⁴t:s²t²(s−t):t⁵−3st⁴) |
| 3A4 | (4),(5)、A²+A−1=0、B=−A−1、C=−A−2、D=1 |
| A6+3A2 | (6) |

ここまでが全parameter・全補集合を含む必要候補の有限化である。
全候補をfreeと数値だけから認定しない。
次の一有限jobで各parametrizationのimplicit quintic、全degree四relationの不存在、
全十八二次gradient核とSaito determinantを確認する。
どれかがfree条件を満たさなければ、候補の除去理由を保持する。
全磁気除去定理はfreeと確認した各curveの全Poisson源を一つのactual23軌道へ送る。
各軌道閉包の周囲成分性・既知帰属は別々に必要であり、成分数を八へ固定しない。
