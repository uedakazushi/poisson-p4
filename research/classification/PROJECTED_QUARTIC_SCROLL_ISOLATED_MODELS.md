# smooth quartic scrollの全isolated-singularity projection：四モデルとrank零quintic排除

2026-11-07、主義務G。
**命題。** smooth quartic minimal-degree surface S(1,3)又はS(2,2)⊂P5を一点からP4へ射影する。
imageがCI22となるcentreを除く全surface orbitは、下の四モデルで尽くされる。
各モデル上で零となるrank四Poisson tensorの実Dは既約被約quinticにならない。
nonCMなextra sectionを含むfull150源を計算し、Jacobi前のfixed quadric factorで排除する。
S(0,4)coneからのprojectionは別の明示補集合である。

## 1. S(2,2)の全centre orbit

normalizationP1×P1のO(1,2)埋込みを
`(s*u²,s*u*v,s*v²,t*u²,t*u*v,t*v²)` とする。
centreはC²⊗Sym²C²のprojective点である。
rank一なら、binary quadraticが平方の場合はcentreがsurface自身に入るので除かれる。
平方でないrank一centreはSegreP1×P2に入り、imageはACMとなる。
この最後の事実は、直接一次資料を読んだM. Brodmann–E. Park,
*On varieties of almost minimal degree I: Secant loci of rational normal scrolls*,
<https://www.math.uzh.ch/d.php/publication/551> Theorem6.2(b)(iii)で確認した。

rank二centreでは二つのbinary quadraticがpencilを作る。
Grassmann(2,3)はdual P2であり、そのPGL2 orbitはannihilator quadraticのrankで二つだけ。
nondegenerateならpencilはbase point freeで、代表をu²,v²に取れる。
rank一ならpencilはあるlinear factorを共有し、代表をu²,u*vに取れる。
fiber PGL2がpencilのbasisを全て変えるので、これらは実際の二whole centre orbitである。
Aut(P1×P1)のdimension六に対し、各centre orbitのdimensionは五、四、stabilizerは一、二である。

それぞれcentre `(1,0,0,0,0,1)`、`(1,0,0,0,1,0)` を選び、
quotient座標を次の表のS22 node/cuspとする。
異なるcentreがAutで移れば、そのquotient imageはPGL5で移る。
従って一つのcentreで全orbitのsourceを計算することは標本検証ではない。

## 2. S(1,3)の全centre orbit

六座標を `(U*s,U*t,V*s³,V*s²t,V*s*t²,V*t³)` とする。
centreをbinary linear block a=(a0,a1)とbinary cubic block C=(c0,c1,c2,c3)へ分ける。
Aut groupはbase PGL2三次元、relative fiber scaling一、
`U→U+V*(alpha*s²+beta*s*t+gamma*t²)` の三unipotent parameterでdimension七。
unipotent部分はcentreのfirst blockへ

`a0→a0+alpha*c0+beta*c1+gamma*c2`,

`a1→a1+alpha*c1+beta*c2+gamma*c3`

と作用する。この二×三Hankel matrixは、Cがpure cubic powerである場合だけrank一となる。
C=0ならcentreはdirectrix S(1)に入り、surface上なので除かれる。

Cがpureでなければrank二で、first blockを全て零に出来る。
残るnon-pure binary cubicのroot型は三相異根又は二重根+単根だけなので、二orbitである。
代表C=(1,0,0,1)と(0,1,0,0)を取り、node/cuspモデルを得る。
centre orbit dimensionsはそれぞれ五、四、stabilizerは二、三。

Cがpureなら、一unipotent方向でfirst blockの一成分を消せる。
残りが零ならcentreはsurface上、非零ならrelative scaling等で一orbitに直る。
この集合は正確にJoin(S(1),S(1,3))\S(1,3)であり、
上記Brodmann–Park Theorem6.2(b)(ii)によりimageはACM。
よって非ACMな全centreは上のnode/cusp二つで尽くされる。

## 3. ACMなdegree四codimension二imageはCI22である

これはsourceの例外を捨てるための仮定ではなく、既済CI22枝へ接続する矢印である。
ACM coordinate ringを三つのlinear parametersで割ると、
embedding dimension二、length四のstandard graded Artin ringとなる。
そのHilbert functionは(1,2,1)だけなので、imageのI2は二次元。
integral・非退化性から二quadricにlinear common factorは無く、degree四CIがimageをscheme-theoretically尽くす。
全CI22のrank零G排除へ接続出来る。

## 4. 四つの全normalization parametrizationとgeneric normal basis

| model | 五座標q | Q=F2 | F3 |
|---|---|---|---|
| S22 node | (s*u²−t*v²,s*u*v,s*v²,t*u²,t*u*v) | x1x4−x2x3 | x0x2x3−x1²x3+x2x4² |
| S22 cusp | (s*u²−t*u*v,s*u*v,s*v²,t*u²,t*v²) | x0x2−x1²+x1x4 | x1²x4−x2²x3 |
| S13 node | (r*s,r*t,s³−t³,s²t,st²) | x0x4−x1x3 | x2x3x4−x3³+x4³ |
| S13 cusp | (r*s,r*t,s³,s*t²,t³) | x0x4−x1x3 | x2x4²−x3³ |

S13ではrのweight二、s,tのweight一を使うdominant polynomial parametrizationである。
normalization全域をweighted projective planeと同一視したのではなく、dense V≠0 chartのU/V=rを使う。
pullbackが全global sectionの零性を判定するにはdominanceで十分である。

F2,F3は各dense chartでnormal idealを生成する。
S22 nodeはx2*x3、S22 cuspはx2*x4、S13二型はx3*x4を反転すれば、
base ratioからF3を得て同じintegral surfaceを切り出す。
F2,F3のglobal CIが含むresidual componentをSと同一視していない。

## 5. full150、extra十一番目のsource、正しいprojective bracket

全150 homogeneous quadratic係数を保持し、全十tripleの `q∧Π(q)=0` を課す。
S22二型では280行、S13二型では全weighted monomialによる300行である。
全四型でrank115、kernel35。
Euler kernel二十四を確認してcanonical quotientを取ると、全projective sourceは十一となる。
quadric×constant skewの十次元だけをansatzとして使っていない。

既約quinticならnonplane secant/Frobeniusによりconormal idealはabelian。
rawΠがEuler商でしかS上に零でないことを保持し、local generators F2/x0²,F3/x0³のbracketを使う。
その正確なhomogeneous numeratorは

`K=x0{F2,F3}−3F3{F2,x0}−2F2{x0,F3}`。

全五first partialをparametrizeして比較する。
S22は全330行、S13は全360行であり、各full sourceへのrankは次の通り。

| model | rank零projective source | conormal rank | 全必要源 | 全Jacobi row rank |
|---|---:|---:|---:|---:|
| S22 node | 11 | 9 | 2 | 1 |
| S22 cusp | 11 | 9 | 2 | 1 |
| S13 node | 11 | 9 | 2 | 1 |
| S13 cusp | 11 | 8 | 3 | 3 |

各matrix、complete kernel、nonzero rank minorを保存した。
必要条件をgeneric chartで課した後、polynomial partialの零性はintegral S全体へ延びる。
chart補集合を別のrank四候補として失ってはいない。

## 6. 全四sourceのactualDに固定Q因子がある

全残源について、全350 Jacobiと全126 actual anticanonical係数を直接構成した。
`CLASSIFICATION_PROJECTED_QUARTIC_SCROLL_FACTORS_20261107.json` は
四つのimmutable source reportをhashでbindし、各々について

`D_actual=Q*explicit_cubic`、exact remainder零

をQQ polynomial divisionで確認・保存する。
この恒等式はJacobi前の全二/三parameter sourceに成立するので、全irreducible quintic枝は空である。

S13 nodeではJacobi idealにz1²、S13 cuspではz1²,z2²が入り、
それぞれ元の全350行によるexact functionalも保存した。
従ってそれらの全被約Poisson源ではD=0となる。
この追加情報はG排除に必須ではなく、fixed Q因子だけで十分である。

## 7. 公開packetと未処理のcone

reusable producerは `code/classification_projected_quartic_scroll.py`、
四source reportは `CLASSIFICATION_PROJECTED_QUARTIC_{S22_NODE,S22_CUSP,S13_NODE,S13_CUSP}_20261107.json`。
各々が当時のcomplete producer bytesを埋込み、後のmodel/field名追加に合わせて旧証明書を書き換えていない。
各modeのjobとfactor jobは `jobs/projected_quartic_*.json` にある。
S22 cusp report末尾のlegacy proseにnodeという語が一つ残るが、model名・parametrization・全sourceはcuspを一意に定め、
factor reportと本noteはその正しいdataを使う。
全四のsource算術は各約6.1–6.5秒、factor/functionalsは約0.096秒で完了した。

残るminimal-degree normalization S(0,4)はrational normal quarticのconeである。
centreをvertex shiftでbaseP4へ移し、binary quarticとして扱える。
base centreがsecant variety内ならsecant locusは一dimensionalでACM/CI22枝へ入る。
外ならsmooth rational quartic curveのconeをimageとするisolated nonnormal vertexが残り、
binary-quartic orbitの一parameterとcollision型をまだ保持する。
この最後のcone family又は全degree四曲面を本packetだけで排除したとは主張しない。
