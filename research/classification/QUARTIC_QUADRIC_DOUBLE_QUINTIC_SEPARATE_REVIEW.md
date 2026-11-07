# rank四quadric上のdegree四曲面と二重quinticの独立読解

2026-11-07、PRIMARY G。producerと異なるrootが全幾何証明と既存whole-centre分類への適用を読解し、受理した。
新しいCAS・minor・jet計算は不要である。

rank四quadric三foldはsmooth quadric Y=P1×P1のprojective coneである。
vertexを除くUはTot(OY(1,1))で、projective bundle resolutionのnegative sectionを除くPicard計算は
Cl(Q)=Z²、H=(1,1)を与える。normal Qからcodimension三一点を除いてもClは不変である。
class(a,b)のeffective Weil divisorはU上の非零sectionを与え、
p_*OU=⊕n>=0 OY(−n,−n)によりa,b≥0が必要。
a=0又はb=0ならsectionは一方のP1のbinary formで、degree四のintegral surfaceにはならない。
degreeはHとの交点a+bである。

従ってintegral quartic Sのclassは(1,3),(2,2),(3,1)。
(2,2)はCartier2Hなので、normal schemeのheight一primeのdivisorial idealを使ってscheme divisor自体が2Hとなる。
0→OP4→OP4(2)→OQ(2)→0 とH1(OP4)=0でambient quadricへliftでき、SはCI22である。
ここはcycleの一致を任意の非CM thickeningへ置き換える議論ではない。
非CIでは残二classだけとなる。

Qがquintic Dを割らないならD|Qはclass(5,5)のeffective Cartier divisor。
S一般点のambient normal regular parameters(Q,t)におけるD∈(Q,t)²をQへ制限するとt²で割れ、
D|Q≥2S。残classは(3,−1)又は(−1,3)となりeffectivityに反する。
従ってQ|D。multiplicity・conormal・零因子のordinary/vertical仮定は一切要らない。

全degree四surfaceの既独立分類はrankを使わないので再利用できる。
S22/S13の非CI全四中心orbitのquadricは表示式から全rank四と直接確認できる。
S04のclosed-secant内はCI22、外はsmooth rational quartic curveのcone。
そのcurveのcontaining quadricがrank三なら、quadric surface coneのCl=ZR,H=2RからC=2HがCartierとなり、
ambient二次式とのCIなのでpa=1、smooth P1に矛盾する。rank二以下も非退化integralityに反する。
よって全外点・全collisionでrank四が保証される。
Veroneseには本quadric議論を流用せず、既全tangent像22の別証明を使う。

以上で全非CI quartic scroll/cone支持は既約quinticの二重支持として空。
既whole projection分類の境界を全て保持した適用であり、一般中心から特殊中心への外挿ではない。
全degree四rank二の統合には、normalCI22空性・non-normalCI22像24上界の別独立receiptを合わせる必要がある。
全Eや上界の完成は主張しない。新規性は未監査。
