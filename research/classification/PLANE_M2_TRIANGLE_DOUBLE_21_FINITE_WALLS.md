# triangle double h=L^2M：全sourceは十fixed壁へ

2026-11-07。主義務G/E。`PLANE_M2_TRIANGLE_DOUBLE_21_GEOMETRY.md` の
二fixed源の空性及びoneparameter generic geometryを、有限壁の全列挙へ sharpen する。
新しいmatrix/profile/GB/jet又はpoint sampleを実行しない。
このpureproofは独立確認前で、十壁のfullhigher G boundはまだ別義務。

## 1. 原wholeP0と全normalization補集合

f=xyz、normal rootu、Z=x partial_z、k20=0、

\[
 Y=\operatorname{diag}(-2r-1,r+1,r),\quad
 X=Y+\operatorname{diag}(-1,2,-1),\quad
 P^{u,v}=\kappa uv,\quad\kappa=3r+1.
\]

actualG2は-3(r+2)x^2yu^2なのでr=-2はm三以上の別義務。
全P1のbad u-rowはc v^2 partial_x、irreducible actualDにはc nonzeroが必要。
normal root degree、off-field sourceのA/B/C除去は既全actuallinear逆shears/scalesを持つ。
以下はr又はcを選んだsampleではなく、任意の全sourceに対する必要条件である。

canonical strict derivationの全weightsは(x,y,z,v;u)について

\[
 (-r-2,\;2r+4,\;2r-1,\;-3r-1;0).
\]

r not in {-2,-1,0,1/2}では全linear u/v off-blockをactualbase shearsで消せる。
bad poleはr nonzeroでdense u nonzero上のrational x'=x+lambda v^2/uで消せる。
lambdaはnonzero constant cのnonzero倍を5rで割ったもので、zero補集合を無断で落とさない。
rational shiftをglobal PGL同値とは呼ばず、polynomial invariantの必要形を読むためだけに使う。

## 2. r nonzeroのpolynomial invariantはS=ux+lambda v^2の非負powerへ

S=u x'、Torus weightsをmu_x、mu_y、mu_z、nu_vと書く。
元Fのnormalization後のLaurent表現はC[u,u^-1,S,y,z,v]に属する。
u最小negative exponentが残ったと仮定し、そのnumeratorをu零に制限する。
Sはlambda v^2へ移る。
固定y/z powersの二monomials S^a v^bとS^a' v^b'がそこでcancelするには2a+b=2a'+b'。
両方weight零なら更に

\[
 (a-a')(\mu_x-2\nu_v)=5r(a-a')=0.
\]

r nonzeroなのでa=a'、b=b'となり、異なるmonomialsのcancellationは無い。
従って最小negative powerは存在せず、**全invariant FはC[u,S,y,z,v]に入る**。
この議論はgeneric rだけでなく、linear-normalizer chart内の全r nonzeroに有効。
未知係数のgenericity、flatness又はunproved extensionは使わない。

## 3. normal order二がhigher monomialを小さくする

deg(u,S,y,z,v)=(1,2,1,1,1)。
monomial u^q S^a y^m z^n v^bは
q+2a+m+n+b=5、normal orderはq+a+b。
leadingS is uxなので、異なるmonomialsの最小normal項は互いに消えない。
指定G2=unit*x^2y*u^2により、唯一normal order二の項はunit*yS^2。
他の全項はnormal order三以上だから

\[
 a+m+n\le2.
\]

従って全higher monomial weightを列挙するには、下の十種類だけで足りる。
u^5は常に許され、a=m=n=0でb positiveならr=-1/3。

| (a,m,n) | bの全範囲 | weight零の全r（重複可） |
|---|---|---|
| (0,0,0), b positive | 1..5 | -1/3 |
| (1,0,0) | 0..3 | -2,-3/4,-4/7,-1/2 |
| (0,1,0) | 0..4 | -2,3,1/2,1/7,0 |
| (0,0,1) | 0..4 | 1/2,-2,-3/4,-4/7,-1/2 |
| (2,0,0) | 0..1 | -2,-1 |
| (0,2,0) | 0..3 | -2,-7,3,1 |
| (0,0,2) | 0..3 | 1/2,3,-2,-1 |
| (1,1,0) | 0..2 | -2,1/2,0 |
| (1,0,1) | 0..2 | 3,-2,-1 |
| (0,1,1) | 0..3 | -3/4,-2,1/2,0 |

各entryは
(-a+2m+2n-3b)r+(-2a+4m-n-b)=0の一式から直接得る。
これが全wallの必要完全性であり、任意の有限r samplesからの推測ではない。

## 4. 壁外全空と二つの別補集合

表に無いrではF=unit*yS^2+A u^5のみ。
A=0は可約なので、actualdiagonal scalesでF=y(ux+v^2)^2+u^5へ戻せる。
既前noteの全AR0/1/2証明から、vertexを戻したquadratic Casimir bivectorはrank二以下。
従ってwall外の全rank四源は空。

r=0ではrationalshiftを使わない。
semisimple weights(-2,4,-1,-1;0)にcommuteするnil部は
N=c v^2/u partial_x+d v partial_z、c nonzeroだけ。
normalleading F2=unit*x^2yu^2へのN像のorder三はnonzero*x y u v^2であり、
z方向はF2へ作用しない。他のhigher項から消せず、この壁も全空。
r=-2は前述のm三以上補集合として保持する。

結局残る全sourceは厳密に次の十rだけに含まれる。

\[
 \boxed{r\in\{-1,\tfrac12,-\tfrac13,-\tfrac34,-\tfrac47,
 -\tfrac12,3,\tfrac17,-7,1\}.}
\]

この十値では原P0へ戻し、全30P1、全9P2、badC零/nonzero及び全higher Jacobiを保持しなければならない。
とくにr=-1,1/2のresonant linear-normalizer chartは上のrational chartで除いており、
そこでoff-blockを消したsliceを代用しない。

## 5. 次のONE finite whole-source testへの含意

二fixed bad-y/z源は別pureproofで全空、oneparameter源はこの十wholeP0以外全空。
よって十P0それぞれの**全fullhigher source**を一batchで判定することは全源分類の有限手順である。
各P0に全30/9を付け、actual filtered normalizerのfiberを戻した像上界24以下、
全sourceactualDの必要因子、又は新Ore認識のいずれかで十分。
十値のselected higher点だけのJacobi検査で全wallを閉じてはならない。
この一batchが未閉なら追加profile/jet/GB列ではなく、残るwhole-source方程式を構造数学で扱う。
concurrent、nonreducedR、m3、他支持は別義務として維持する。
