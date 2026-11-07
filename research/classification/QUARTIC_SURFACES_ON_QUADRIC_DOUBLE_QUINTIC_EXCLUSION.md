# rank四quadric上の非CI quartic曲面：二重quinticはquadric因子を持つ

2026-11-07、主義務G、producer plane_collision_review。
**命題（初等幾何のproducer証明）。** rank四quadric三fold Q⊂P4に含まれる
integral非退化degree四曲面SがCI(2,2)でないとする。
quintic hypersurface DがSの一般点でmultiplicity二以上なら、QがDを割る。
従ってSは既約quinticの二重曲面supportにはなれない。
Poisson restrictionのrank、Sのnormality、零divisor Rの型を仮定しない。

## 1. projective quadric coneのclass groupとeffectivity

rank四Qはsmooth quadric Y=P1×P1のprojective coneで、vertexは一点v。
normalityはhypersurfaceのS2と、singular pointがcodimension三であるR1から従う。
U=Q\{v}はprojection pによりTot(OY(1,1))と同型である。
例えばresolution P_Y(O⊕O(1,1))のtautological bundle ξの全五sectionがQへ写し、
quotient Oに対応するnegative section E0だけをvへ潰す。
E0のclassはξ−p*OY(1,1)なので、これを除いたPicard groupはPicYとなる。
normal Qからcodim三一点を除いてもWeil class groupは変わらず、

`Cl(Q)=Pic(U)=Z⊕Z, H=(1,1)`。

各hyperplaneでvを避ければYのgraph sectionとなり、class(a,b)のsurfaceのdegreeはa+bである。

effective Weil divisorのclass(a,b)について、その非零regular sectionをUへ制限する。
pはline bundleのaffine projectionなので

`H0(U,p*OY(a,b))=⊕_(n>=0) H0(Y,OY(a−n,b−n))`。

従ってeffectiveならa,b≥0。逆にこの領域ではbase sectionのpullbackからeffective divisorを得る。
a0の場合はn0だけが残り、sectionは第二P1のbinary formとなる。
そのdivisorはruling planesの和であり、degree四のintegral surfaceにはなれない。
b0も同様。

## 2. degree四・非CIならclassは(1,3)又は(3,1)

a+b=4かつa,b正なので、候補は(1,3),(2,2),(3,1)。
もしSのWeil classが(2,2)=2Hなら、Cartier divisor2Hとlinearly equivalentなのでS自体もCartier。
normal schemeのheight一primeはdivisorial idealであり、同じcycleだけの別のnonCM schemeをここで選んでいない。
従ってSはOQ(2)のnonzero sectionのzero schemeである。

hypersurface restriction exact sequence

`0→OP4→OP4(2)→OQ(2)→0`

とH1(P4,OP4)=0により、そのsectionはambient quadric q'へliftする。
ゆえにS=Q∩(q'=0)はschemeとしてCI22となる。
仮定よりこの候補を除き、class(S)は(1,3)又は(3,1)だけである。

## 3. double条件とquinticのclassは両立しない

QがDを割らないと仮定する。
そのrestrictionはOQ(5)の非零sectionであり、effective Weil divisor D_Qのclassは(5,5)。
Sの一般点はQregに入り、ambientの正則normal parametersを(Q,t)に取れる。
D∈(Q,t)²からD|Q∈(t²)、従ってWeil divisorとしてD_Q≥2S。

class(S)=(1,3)なら残りD_Q−2Sのclassは(3,−1)、逆なら(−1,3)。
どちらも§1のeffective coneに入らない。矛盾。
従ってQがDを割る。一般normal multiplicityが二を超える場合も同じ必要条件に入り、
ordinary nodeやvertical成分に関する仮定は無い。

## 4. 既存whole projection分類への正確な接続

`DEGREE_FOUR_RANK_ZERO_SURFACE_EXHAUSTION.md` §1の外部分類入力はrankを使わない。
その全integral・非退化degree四surfaceはACM del Pezzo又はminimal-degree Y⊂P5の外点射影である。
ACM像は全CI22。Yの全四型と全centre補集合は同note §§2–5で既に別作者確認済み。
その後のrank零source証明だけを本命題へ置き換え、分類自体は再計算しない。

**S22/S13の非CI全四orbit。**
`PROJECTED_QUARTIC_SCROLL_ISOLATED_MODELS.md` §§1–4は全centreを分類し、
ACM/CI centres以外をnode/cusp各一orbitへする。
その全四modelのcontaining quadricは

```
S22 node: x1x4−x2x3;
S22 cusp: x0x2−x1²+x1x4;
S13 node/cusp: x0x4−x1x3.
```

全てrank四である（それぞれ四変数のinvertible quadratic matrixと一radical）。
非CIであることも同whole-centre分類に含まれる。
従って各modelの全rank二restriction枝、全R、全special parameterに本命題が適用出来る。

**S04の全連続族とcollision。**
centreがclosed Sec(C4)内なら既存whole criterionからACM/CI22へ入る。
外ならimageはsmooth rational quartic C⊂P3のconeである。
`PROJECTED_QUARTIC_CONE_RANK_ZERO_EXCLUSION.md` §1のrankを使わない証明により、
Cのunique containing quadric Q0はrank四smooth quadricである。
rank二以下はintegrality/非退化性に反し、rank三cone上ではCl=ZR,H=2RからC=2HがCartier、
ambient quadric sectionへliftしてpa一となりC≅P1に反する。
従ってSのcontaining cone Qはrank四三fold。
この主張はgeneric centreだけでなくclosed Sec外の全binary-quartic collision型を含む。
Sec内のtangent/secant limitsはCI22補集合へ明示的に戻す。
rank三containing quadricの未確認枝を勝手に本命題へ含めていない。

**Veronese。**
centreがclosed Sec内なら全CI22、外なら唯一smooth projected Veroneseであり、
その全tangent Poisson源は既独立確認済み `PROJECTED_VERONESE_TANGENT_NONDOMINANCE.md` によりPGL像≤22。
本命題をquadric無しのその枝へ適用しない。

projectionの境界がなおintegral・非退化degree四surfaceなら、全分類で再び上のいずれかに入る。
一般centreのclosureだけで全境界を覆ったとは推論していない。

## 5. degree四rank二Gへの帰結

全normal CI22 rank二空性は独立受理済み。
全non-normal CI22 rank二像≤24は別producer packetが完成し独立review中。
この二つとsmooth projected Veronese全tangent像≤22、本命題による残projection全排除を合わせれば、
全integral非退化degree四supportのrank二G像は二十四以下となる。
non-normal CI22の独立review前には、この最後の統合結論をそのreceipt付き確定claimとはしない。
本命題そのものはその計算に依存しない。

hyperplane内degree四surfaceは既済のprime CI symbolic-square次数論でH因子となる別枝である。
本稿は全rank二点のE帰属、成分数、一般五点又は分類全体の完了を主張しない。
新規性・先行文献での同じlemmaの有無は未監査。
