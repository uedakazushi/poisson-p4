# μ(F)=6 rational cuspidal quartic：r2 marked geometryの有限全リスト

2026-11-07、PRIMARY E。producer手証明、独立確認前。
`PLANE_HYPERPLANE_QUARTIC_MARKED_ENTRANCE.md` のmu六枝を全projective geometryで処理する。
FはA6,E6,A4+A2,3A2の四型、H∩Fは高々二点である。
本稿は全marked pairのリストであり、各sourceの存在・cofactor型・磁気商・成分帰属は
別の一括Saito/source検証に残す。mu五の非cuspidal profileは別義務のまま。

## 1. A6：全quarticと全r2 line

cusp逆像をs=0に取る。tangent h=0の接触四からh=s4。
transverse sectionをP1のt→t+lambda*sとbase shearでy=s²t²へ直す。
第三sectionのh,y成分を除けばw=a s³t+b st³+c t4、c!=0。
local h−c*y²の最初の奇数項はb*s5、次にa*s7となる。
A6であるためb=0,a!=0。非零scaleで全quarticは

```
[h:y:w]=[s⁴:s²t²:8s³t+t⁴],
F=(hw-y²)²-64h³y
```

へ移る。全global familyを正規化したもので、残るmodulusを捨てていない。
general lineはpullback A s4+B s²t²+8C s³t+C t4。
cuspを通るC=0では、r≤2はh（r1 tangent）又はyだけ。
A,Bとも非零ならs²(A s²+B t²)の三つのrootが生じる。

cuspを避けるC!=0でA=0なら∞が単根であり、残るcubicのtriple rootはlinear coefficient零と
nonzero constantに矛盾する。従ってr2ならA!=0。
affine r=s/tの二root partition(3,1)ではlinear coefficient零からv=−u/3、

```
(r-u)³(r+u/3)=r4-(8/3)u*r3+2u²*r2-u4/3.
```

coefficientを比較するとu³=1。
三smooth flex lineはr→zeta*r,zeta³=1という実base automorphismで同一marked orbit。
u=1を選べば H=h+2y−w/3。
(2,2)ではrootがu,−uとなりr³ coefficient零からC=0に反する。
従ってr2の全marked orbitは **y、h+2y−w/3** の二つ。

## 2. A4+A2：全四r2 line

全quarticの正規化と全root比較は `PLANE_HYPERPLANE_QUARTIC_A4A2_MARKED_LINES.md` による。

```
[h:y:w]=[s⁴:s²t²:12st³-6t⁴],
F=(hw+6y²)²-144hy³.
```

全r2 marked lineは

```
y,
w,
w-6y,
h-(32/3)y+(128/81)w.
```

r1はhだけで、Ga9power3の既済geometryへ戻る。

## 3. 3A2：全quarticと七line、三marked orbit

三cuspは非collinearである。joining lineが三cuspを通ればintersection六となるため。
正規化上の三逆像をs=0,t=0,s=tへ移す。
各joining lineは二cuspへtransverse、pullbackはその二点を二回ずつ持つ。
よって全quarticは非零base scalesを除き

```
[h:y:w]=[s²(s-t)²:t²(s-t)²:s²t²],
F=(hy+hw+yw)²-4hyw(h+y+w)
```

へ移る。座標を置換する実S3が三cuspを全て置換する。
general lineのaffine pullbackは

```
A*r4-2A*r3+(A+B+C)*r2-2B*r+B.
```

cusp h-axisを通るA=0では、二つの残rootがdoubleになるB*C=0（joining lines）か、
B+C=0による∞のtriple root（cusp tangent）だけがr2を与える。
S3で全cuspを処理すると、joining lineはh,y,w、cusp tangentはy−w,w−h,h−y。

cuspを避けるdouble-doubleではroot u+v=1,uv=B/A、constantから(B/A)²=B/A。
cuspsを避けるのでB/A=1、さらにC/A=1。唯一のsmooth bitangentはh+y+w。
triple-singleならv=2−3u、linear/constant比較からu²(u−1)²=0となり、cuspを避ける仮定に反する。
smooth flexは無い。従って全r2 geometryは七lineで尽くされ、S3を戻すと代表は
**h、y−w、h+y+w** の三marked orbit。
r1は無い。cusp tangent contact三、smooth hyperflexも上の比較で存在しない。

## 4. E6：二つの全quartic orbitとmarked line

cusp逆像をs=0、接線sectionをY=s4、transverse sectionをX=s³tに取る。
第三sectionからX,Yを引き、P1のt-translationでst³係数を消すと、全quarticは

```
[X:Y:W]=[s³t:s⁴:t⁴+a s²t²],
F=Y³W-X⁴-a X²Y².
```

a=0とa!=0は別に保持し、非零ならactual diagonalでa=−6に直す。
Hがcuspを通るtransverse lineの全familyは `PLANE_HYPERPLANE_QUARTIC_E6_TRANSVERSE_SOURCE.md` で処理済み。
tau12条件はH=X（他のline X+bYではsocle coefficient bが非零）を強制する。
a=0では全quadratic anchorが比例してframeが無く、a!=0のH=Xは全六magでrank一cofactor。
したがっていずれもtype6 rank二nil22残余には入らない。

cuspを避ける一般Hは B Y+A X+C W、C!=0。
pullbackはB r4+A r3+aC r2+Cである。
finite triple rootuならlinear coefficient零からv=−u/3、constant/r²比較により

```
a*u²=-6.
```

a=0には解が無い。a=−6ではu=±1、r→−rのbase automorphismで同一orbit。
u=1の代表は **H=8X−3Y+W**。
double-doubleならrootはu,−u、A=0、coefficient比較によりB/C=a²/4。
a=−6の代表は **H=9Y+W**。
∞の単根・二重根caseも保持すると、B=0,A!=0は残るcubicのtriple rootに矛盾し、
A=B=0はa!=0でr3となる。a=0ならsmooth hyperflex r1でありtau十三なので対象外。
E6のa=0には残るr2候補は無い。a!=0では上の二つだけがtype6未filter候補。
r1 cusp tangentYはa=0でtau十三、a!=0で既済Ga9linear-w geometryへ戻る。

## 5. 一括source検証へ渡す全十一代表

base変数を(h,y,w)へ統一し、次を**同じ有限batch**で扱う。

| F | marked H | 代表数 |
|---|---|---:|
| (hw-y²)²−64h³y | y、h+2y−w/3 | 2 |
| (hw+6y²)²−144hy³ | y、w、w−6y、h−32y/3+128w/81 | 4 |
| (hy+hw+yw)²−4hyw(h+y+w) | h、y−w、h+y+w | 3 |
| y³w−h4+6h²y² | 8h−3y+w、9y+w | 2 |

全r2代表はlocal unionがsimple ADEの和となりtau十二だが、それだけでfreeを宣言しない。
元の全18gradient unknown、Saito determinant、actual H→h GL3と逆、全六magnetic、
canonical affine cofactorのrank/Jordan/constant imageを直接確認する。
E6のa0 no-frameと非零a rank一controlも別行で保持する。
sourceの全orbit/magnetic正規化や全親membership、個々のambient componenthoodは別の次課題である。
