# 残plane二H＋cubic三源：全prior142比較と全成分性の独立純査読

2026-11-08、`/root/generic_exhaustion`、PRIMARY E/component audit。
**独立純読解・別手比較 PASS。** 新CAS、point、operator、minor、profile又は変形jetは無し。
producerのf37cb193は保存時T37を条件付きとしており、そのbytesを変更しない。
今回受理済みのwhole T37単純component定理と三targetの非代数葉を合成して、その最後の条件を満たす。
catalogue/Git/mainは編集しない。schemeのsmoothness／unique local branch／nil厚みは未主張。

## 1. 全三閉像とexact23を独立確認

source IDsは

```
node_flex_tangent_through_flex,
cusp_flex_secant_through_flex,
cusp_cusp_tangent_cusp_secant.
```

base(h,x,y)でcubicは順に

```
C1=h*y²-x*(x-h)²,
C2=h*y²-(x+h)³,
C3=h²*y-(h-x)³.
```

全R6 sourceと全GL2 Saito anchor frameを保ち、そのPGL飽和閉包をX1,X2,X3とする。
全sourceのactualDはhxCiで非零、sourceは既約だから各Xiは既約閉集合。
元fullsource theoremはfiber-dependent二次anchorもSaito moduleで戻しており、
選択した磁気sliceだけから全Xiを仮定したのではない。

固定projectiveDでは二次Saito二本の定数GL2 frame四係数とR6が全tensorを尽くす。
一意のh-Casimir liftのmixed coefficientsでframe、fiber-fiberのbase quadraticでRを復元出来る。
projective tensor scalarの一次元だけを引き、fixedD fiberは9次元。
二factor h/xを個別に保つconnected base stabilizerは
δh=ah,δx=bx,δy=cy+dh+exと書ける。
各Ciの全係数比較を別に行うとa=b=c,d=e=0、Euler scalarだけが残る。
二H交換は有限なのでbase PGL3 stabilizerは有限、orbit8。
essential3 subspace Gr(3,5)はdim6、Dorbit14と合わせて各Xiはexact23。
これは旧全sourceF48下界とも一致するが、新rank計算はしていない。

## 2. immutable142比較の漏れを読取検査

producer projection b662a0afは142個の相異ID、旧5369e9ddの126 IDを全て含み、差は16。
差は新三irreducible-quintic23と十三marked-multipleH25である。
read-only Git blob702baab17cc659cb6583dc5ea470311dccd875fdの全catalogue bytesのSHAが
32a37a4ec226f386fe81e0898dba3c5a0c4ab3c8b68216144bf69d97c7e2b3e6と一致し、
142各行のgeneric rank/projective dimension/actualD fieldがprojectionと一致することを確認した。
動くcurrent catalogueの件数を過去142へ戻したのではない。
15input manifest572bda73の全byte hashも一致する。
この機械読取は台帳の整合性だけであり、以下の数学をCASで検証したという意味ではない。

## 3. 全priorの必要な群と非包含の論理

|群|件数|確認した非包含の根拠|
|---|---:|---|
|rank二|6|rank≤2がclosed|
|rank四dim<23|18|Xi exact23|
|MPS23|29|Xi wholeD essential≤3のproper Grass像、同次元親はSNC essential5点を持つ|
|非MPS23でgenericDの型が異なる|46|包含なら等号、dense genericDのprime degree/multiplicity/既約性が異なる|
|既plane二H＋cubic23の九配置|9|同次元等号ならdense配置像が交わるが、内在的marked配置が異なる|
|nonreduced又は五Hの高次元閉因子親|20|proper repeated-H又はfive-linear multiplication像、XiのDはreducedでirreducible cubicを持つ|
|最大factor次数≤2の24親|7|proper multiplication像、irreducible cubicはその像に入れない|
|T39/T40|2|X1/X2の旧同点tangent24、X3のIunit fullfaithful graph恒等的ideal(1)|
|T21|1|proper ordered-two-factor modular/contact全source条件違反|
|Ga9 cuspidal24|1|proper marked fourthpower/quadraticcone/binarydiscriminant全境界条件違反|
|T36/T38|2|proper order二zero/R22条件及びwholezero-surface条件に各Xiのdense源が反する|
|T37|1|今回の単純component low-degree surface条件とtarget no-algebraic-leaf|

和は142。上の各群はprojectionのID partitionと一致する。

同じdim23の比較は**全閉像包含なら等号**を使う。
ある特殊磁気点が既知成分と交わらないという主張へ強めない。
MPS23のSNC点は既fullsourceで全edge零のactual tensor、単なるabstract graph名ではない。
proper Grass像でXiのwhole essential≤3を保つため、SNC点一つでも同じ閉集合には入れない。

追加十三25については各whole sourceが元multipleH仮定h²|Dを保ち、そのmarked-factor flagをproperに忘れる。
したがって全componentがrepeatedlinear因子を持ち、XiのreducedDとは両立しない。
M1/M11/B1/B3/V1のscheme smoothnessを仮定せず、既成分性・次元・閉因子条件だけを用いる。

## 4. 高次元親／全boundaryの範囲を別に監査

T39/T40の第一二源では同じ旧witnessのD101からprojective tangent24。
25次元componentはその点を含めないので、wholeXi包含も不可能。
第三sourceのD96/tangent29へこのargumentを流用しない。
全Iunit R6→whole QQ(lambda) normalizer/actual inversesとfaithful global134/45 graphの
独立receipt7136d980のconstant残差を使い、wholeXiがT39/T40へ入る可能性だけを除く。
I0のboundary個別外性は主張しない。

T21のproper factor flagはtargetの二simpleHとirredC3で正確に二Hへ戻る。
H交換も両必要条件を交換するだけである。
全source条件はC|x0がh³、C|h0がx³の倍を要求するが、三Ciは満たさない。
磁気係数に依存せず、zero磁気もこのflag違反を回避しない。

Ga9 cuspidal24は単なる一witness tangentを使わず、全marked factor/p/jet/discriminant flagを使う。
唯一可能Hについて第一二はquadraticconeがZ(X±Z)でH²でなく、
第三はbinarydiscriminant h²x²に三重根がない。
zero leading coefficient／discriminant0もproper incidenceに含む旧全証明を読解した。

T36の全rank四boundaryには「あるorder二zeroのhorizontal cubicがquadric-pencil閉coneR22内」
というproper条件がある。各targetの全候補q/dense磁気における非零one-formは
invariant irreducible C3を持ち、非constant quadric ratioはこのC3上で定数に出来ないためR22外。
T38のwholezero-surface条件と、各dense lambda源のzero dimension≤1の比較も正しい方向である。
ここでも有限lambda例外/I0/追加zero-planeの個別membershipを除いたとはしない。

## 5. T37最後の条件を今回満たす

source側のexact theoremは
`T37_LOW_DEGREE_LEAF_SURFACE_SIMPLE_COMPONENT_SCOPE.md` と独立receipt
`T37_SIMPLE_COMPONENT_LOW_DEGREE_LEAF_AND_GA9_SEPARATE_REVIEW.md`。
actualDのmultiplicity一成分上で、D自身がsmooth、πrank二の一般点はdegree≤5 invariant surfaceを持つ。
三targetのD=hxCiはreduced、Ciの一般点はh*xunit、smooth/rank二であり、このscopeに入る。
nonreduced componentのreduced supportだけへ誤適用しない。

全磁気正常化後のnode源は、tauによる正規化x=tau²,y=tau(tau²−1)上で

```
omega=dq+J*dlog(tau)+(I/2)*dlog((tau-1)/(tau+1)).
```

cusp-flex源はx=t²−1,y=t³上で

```
omega=du+(I+J)/2*dlog(t-1)+(J-I)/2*dlog(t+1).
```

第三cusp接線割線源はt=1−x,y=t³上で

```
omega=du−3(lambda+1)*dlog(1−t).
```

これら三式を印刷された全restriction bracketsから別に手再導出した。
第一二の(I,J)≠0、第三のIunit/lambda≠−1では少なくとも一つの留数が非零。
free fibre方向のHamiltonian invarianceにより、integral algebraic invariant surfaceがあれば
q又はuはbase C(t)の有限代数拡大でmeromorphic functionとなる。
projective normalization上でexact differentialのresidueは零だが、
pullback logarithmic residueは正のramification index倍され非零、矛盾。
整数比／有理比のresidueを別に除外せず、全有限coverで成立する。

したがって各Xiのdense sourceにT37の必要surfaceは無く、Xi⊄T37。
元I0/J0/lambda−1のclosed磁気boundaryはそのままXiに残す。
このdense noncontainmentだけでwholeXi親の排除には十分である。

## 6. 全ambient成分性と相互区別

任意のXiを含むP_redの既約成分Zを取る。
§§3–5の全142比較によりZはknown142ではない。
既独立 `KNOWN142_SEVEN_CANDIDATE_GLOBAL_COVER_SEPARATE_REVIEW.md` は、
全rank四／rank二点の有限closed被覆から、未採録全componentのprojective dimension≤23を証明している。
Gのbeta-image次元を全Eへ外挿した上界ではない。
従ってdimZ≤23。一方Xiは既約閉／exact23なのでXi=Z。
これで**三Xiはいずれも全Poisson支持の既約成分**である。

第一のunique cubicはnodal、他二つはcuspidal。
第二の二factorlinesはcuspを通らず、第三の二factorlinesはcuspを通るので、
各denseD配置のorbitは互いに異なる。
等しいXiのdense constructible D画像は同じ非空openで交わるはずだから矛盾し、三成分は相異なる。
Ga9linearw/constantDOUBLE候補のgenericDはそれぞれ1+4／nonreduced2H+C3で異なり、
同じ23次元閉包とは等しくなれない。constantTRIPLEは22次元なので等しくない。
後でrootがそれらを成分と認識しても三Xiとの重複加算は生じない。

全R6 source／GL2 frameのclosed像なので、I0/J0と全lambda例外を含む全点も新三成分へ入る。
これをもってspecial boundaryが既知親と交わらないとはしない。
新成分性は被約支持の極大性だけであり、smoothness・localunique・一般nil厚みは依然不明。
旧D101/101/96や局所二次/三次nonclosing結果を変更しない。
rootだけが新三ID・count・主稿を統合する。
