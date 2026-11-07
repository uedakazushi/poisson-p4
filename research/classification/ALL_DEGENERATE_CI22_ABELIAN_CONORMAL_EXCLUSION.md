# 全memberが退化するintegral CI22 pencil：rank零・既約quinticは空

2026-11-07、主義務G。
S⊂P4をintegral・非退化なCI(2,2)、defining pencilの全memberを退化quadricとする。
**命題。** S上で零となるrank四Poisson tensorの実Dは、既約被約quinticにならない。
以下は全symmetric singular-pencil型を含み、generic pencil又は一つのparameter点を選ぶ議論ではない。

## 1. 使うcanonical formの一次資料と正確な範囲

F. De Terán–A. Dmytryshyn–F. M. Dopico,
*Generic symmetric matrix pencils with bounded rank*, 著者公開PDF
<https://gauss.uc3m.es/fdopico/papers/jst2018.pdf> を2026-11-07に直接読んだ。
§2、Theorem2.1（PDF page5/index4、直前page4のblock定義）は、
**任意のcomplex symmetric pencil**をcongruenceでregular symmetric Jordan blockと
singular blockM_dの直和へ一意にするcanonical theoremを述べる。
`M_d=[[0,L_d^T],[L_d,0]]` のsizeは2d+1で、M0は一×一zero blockである。
regular/singular部分を同時に含む全型の定理であり、同論文のgeneric closure theoremを
全pencilへ誤って拡張して用いてはいない。
同箇所が参照する原典はR. C. Thompson,
*Pencils of complex and real symmetric and skew matrices*,
Linear Algebra Appl.147(1991),323–371,
<https://doi.org/10.1016/0024-3795(91)90238-R>。
原典の全証明を本作業で再証明したとはせず、直接読んだall-pencil theoremを外部入力として明記する。

以下、size五なのでM_dはd=0,1,2だけである。
integral・非退化Sを含む非零quadricはrank三以上であることを常に使う。
rank二以下ならprime I_Sがlinear factorを含み、非退化性に反する。

## 2. 全singular-blockの有限場合分け

### 2.1 M2だけの場合

標準形は
`q1=x0x3+x1x4, q2=x1x3+x2x4`。
CIはplane{x3=x4=0}を含む。
二次元integral CIが二次元planeを真に含むことは出来ず、Sの非退化性に反する。

### 2.2 M1とregular二次元部分の場合

M1の二式はx0x1,x0x2。
regular二次元部分が二つの相異なるeigenvalueを持つなら、
pencil PGL2と同時のx1,x2変更により全形は

`q1=x0x1+x3², q2=x0x2+x4²`。                         (1)

単一固有値でdiagonalizableなら、regular部分を消すpencil memberが
x0*linearとなりrank二で不可能。
単一Jordan2ならcanonical formは
`q1=x0x1+2x3x4, q2=x0x2+x4²`。
これもplane{x0=x4=0}を含むのでintegral非退化CIではない。
従ってこのblock型で残る全モデルは(1)だけである。

### 2.3 M0一つとregular四次元部分の場合

これはcommon radical一次元を持ち、q1,q2は四base座標だけに依存する。
regular四次元pencilの非退化memberでQ,Mを用いる。
任意の固有値λについてrank(q2−λq1)≥3なので、四次元Mの各eigenspaceは一次元。
従って各primary partは一つのJordan blockである。

全rank零sourceは `Π=P_base+∂c∧V` と分かれる。
conormal-abelian条件はbase constant skew matricesにだけ作用し、
vertex側は自由である。
`REGULAR_PENCIL_CI22_ABELIAN_CONORMAL_NONDOMINANCE.md` と同じcommutator計算を四次元に適用する。
r≥2ではskew centralizerが零なのでP_base=0、従ってΠはdecomposableでrank二以下。
r=1では単一Jordan4で、base normalizerはgrading一次元だけ。
従って `P_base=alpha(q1,q2)*A` であり、
`Π∧Π=alpha*(alpha*A∧A+2A∧∂c∧V)`。
実Dはquadratic alphaで割れるか零なので、既約被約quinticにならない。

### 2.4 M0一つ、M1、regular一次元の場合

common radical座標を除く全形は
`q1=x0x1+a*x3², q2=x0x2+b*x3²`、(a,b)≠(0,0)。
regular部分を消す `b*q1−a*q2=x0*(b*x1−a*x2)` は非零rank二quadricであり、不可能。

### 2.5 common radical二次元以上の場合

残る三次元以下のpencilは、あるnonzero memberでrank二以下となる。
三×三のhomogeneous determinantは、非零ならP1上にrootを持ち、恒等零なら全memberがrank二以下である。
二次元以下では最初からrank二以下である。
従ってこの全場合も不可能。

以上でsize五の全singular-block分割を尽くした。
M1+M0+M0は§2.5にも含まれ、他の複数singular blockはsize又はcommon-radicalの場合に入る。

## 3. 唯一のM1+regular二固有値モデルの全source

(1)は実際にintegral CIである。
x0≠0 chartはx1,x2を消去するirreducible smooth openであり、
x0=0 boundaryは{x0=x3=x4=0}のlineだけ。
complete intersectionのunmixed性により別のsurface componentは無く、
generic reduced性からidealはprimeとなる。

全ΠはKoszul/Eulerによりq1 A+q2 B。
secant/Frobeniusによりconormal idealはabelianで、prime CIのfree conormalから
`h_A=dq1 A dq2∈span(q1,q2)`、Bについても同じ条件を得る。

Aijをi<jの十係数とすると、h_Aの全式は

`−A01*x0*x2+A02*x0*x1−2A03*x2*x3+2A04*x1*x4`

`+A12*x0²+2A14*x0*x4−2A23*x0*x3+4A34*x3*x4`。

このquadraticがspan(q1,q2)に入るには、x3²,x4²の係数比較からA01=A02=0、
他のmonomial比較からA03,A04,A12,A14,A23,A34も零でなければならない。
A13,A24だけが全て自由に残る。
従って全必要sourceは正確に

`Π=(a*q1+c*q2)∂1∧∂3+(b*q1+d*q2)∂2∧∂4`。

全imageは固定四次元空間span(e1,e2,e3,e4)に入り、実Dは直接

`D=−x0*(a*q1+c*q2)*(b*q1+d*q2)`

となる。
これはJacobiを課す前から既約被約quinticではない。
係数の一部が零の場合もDが零となるだけであり、別のrank四irreducible枝を作らない。

## 4. 全CI22に対する到達点

全member退化pencilのrank零・既約quintic枝は上の全場合で空である。
非退化memberを持つpencilの全Jordan型は別noteでG像dimension≤24まで処理済み。
従って、独立review後は**全integral非退化CI22 rank零曲面**を未知G一般元から排除出来る。
この結論は全CI22上のPoisson tensorが零という主張でも、各可約D源の全E帰属でもない。
