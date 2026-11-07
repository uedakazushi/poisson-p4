# R=2L+M・m3 squarefree：全十九源の別作者検証と二残余

2026-11-07、別作者 `/root/quartic_review`。**全十九源の独立算術・全finite geometryはPASS**。
G全完了の意味ではない。十六whole SF枝は空、一c-dependent枝は全像≤24、
`DIAGONAL_-2_-2` と `NIL_SPLIT` の二whole源は本receiptでは未解決のまま。

入力はactualD被約既約quintic、restriction rank2 plane、zero cubic R=2L+M、L≠M、
normal位数3、generic squarefree binary cubic G3。
R=3Lの既済二モデル排除は再計算せず、m3 repeated-root、m2、reducedR及び他Sは別義務。
producerの新Ore手証明もこの十九源packetのscope外である。

## 保存された正しい実行と独立実装

Cloud resetで失われた138f21ee…を同じreport identityとして偽復元しない。
**NEW** `CLASSIFICATION_PLANE_DOUBLE_LINE_M3_SQUAREFREE_REPLAY_20261107.json`
（bf8be1a3…）、実行source95981c68…のexact bytesを入力とした。
producerのfresh runtimeは22.097秒。歴史runtimeをmemory内だけ戻した旧hash照合receiptはprovenanceであり、
本独立検証はfresh値と現存する全mathematical dataを使う。

checker `classification_plane_double_line_m3_review.py` はproducer/CASをimportせず、
既別作者のFraction・疎多項式演算だけを用いる。
有限job gate PASS後4.412秒、42,376KiBで全十九を完了した。
新report `CLASSIFICATION_PLANE_DOUBLE_LINE_M3_SEPARATE_REVIEW_20261107.json`
には実行reviewer全文、全input hash、全case rank/minor/判定を保存する。

全sourceはP1 raw36/Euler6=30、P2全9。生の全basisとEuler quotientのrankを再構成した。
producerのmixed列順d,b,e,cは独立basisの明示置換で合わせ、変数名による推測をしない。
全actual equationsは次を含む。

* global K三未知の**全30係数affine系**、rank3/augmentedrank3、unique K。
* 全150 P0係数、全350 Jacobi、実Euler-Pfaffian零。
* tracefree mixed fieldsからの3fβを再構成し、root geometryからの全connection numeratorとも照合。
* 全180×24 parallel PDE、完全kernel及び全membership annihilator。
* 全350×30/350×9 operator、全24×30 actualD3 map。
* weight1のJacobi+parallel条件の全matrix、全K1/K2、小行列式、完全actual shear。
* 全D3image matrix、完全image basis、その**全係数不定元**上のbinary-cubic discriminant。

各rankは非零minorの下界と独立完全kernelの上界を合わせて証明した。
補助annihilatorも全kernelを殺すこととrankを検証し、chosen subsetへのmembershipで代用しない。
D3はΠ0/Π1のcrossだけに依存し、Π1²/Π0Π2は次normal orderなので対象から漏らしていない。

## K手予測の明示訂正

旧 `PLANE_DOUBLE_LINE_M3_CANDIDATE_LIFTS.md` のaffine quotient bracketの後二項符号は誤りだった。
正しくは
`{z/b,w/b}=Pzw/b²+w Pbz/b³−z Pbw/b³`。
独立実装はこの手予測を使わず、globalEuler-Pfの全K系から解いた。

| candidate | 正しいunique global K |
|---|---|
| C=a²,g2 | 2z² |
| C=a²,g3 | 0 |
| C=a²,g4 | −2z² |
| C=b²,T=a²b² | 0 |

全四つともglobal Kが存在し、全P0Jacobi/Pf零。
第四のhistorical kind文字列 `cdependent_inconsistent` は旧予想名であって数学判定ではない。
第四をK不整合で除外せず、後記の**全actualD3image**でSF空と判定する。
旧手note bytesを改変せず、erratum2778dbcc…を明示的に依存とする。

## 全root-cover/irregular/resonance幾何の独立読解

U=P2\(L∪M)=Gm×A1、π1=Zなのでdegree3root-cover monodromyは1、C2、C3の全三つ。
t=s^d、d1/2/3のKummer coverでsplitする。
C[s,s^-1,r]はUFD/Pic0、primitive root vectorsの全pair determinants/overall factorはunit。
r-degreeの和≤2、ただ一つpositiveなら他二constant independent rootsとのunit determinantsに反する。
従って全rootはr-independent、又はdegrees(0,1,1)。
後者の唯一degree0rootはdescendし、global primitive factorのbase degreeが正なら残quotientはbase degree≤1で
c²項を持てないため矛盾。よって**constant normal root**である。
connected C3は全degree同じなのでr-independentだけ。split/1+2を捨てない。

r-independentではtwisted stable cubicのLie stabilizer零から全connection一意、dr成分零。
全P0はc-independent suspension、ordinary β=(C/t²+D/t)dt。
parallel discriminantのdet² characterでtrC0。
nonzero semisimple Cは三finite Riccati rootsのpair log-difference poleとt^-3係数が矛盾。
infinity rootがconstantならE項零、残二rootのdifferenceにdouble poleが残り同じ矛盾。
nonzero nilpotent Cでも三finiteの場合E項をpole比較で消し、linear Riccatiへ戻す。
以後のfinite rootsがC*上poleを持たないことは**regular linear ODE**から従うのであって、
一般projective Riccati rootがinfinityを横断しないと仮定していない。

一般linear解はc0/t+d t^e、c0≠0。e−1はnonzero log term、e0のlog係数もalgebraic性で零。
三finite rootsならw³係数t^p、p≤2に対するz³係数c0³t^(p−3)をe≠−1では消せず矛盾。
ゆえにinfinityがconstant root。残二根の全polynomial条件はp2、splitならe0、C2ならe−1/2のみ。
constant normal unipotentによるresidueの余分off-diagonal消去もactual linear変換である。
この全二nil caseをpacketは保持する。

C0ではfinite projective monodromyがnilpotent residueを排除し、constant GL2でM1 diagonal。
squarefree regular normal cubic supportは必ず(0,2),(0,3),(1,2),(1,3)の一pairを含む。
全p,q∈{0,1,2}との二weight equationsを独立に解き、**全36 intersections**とnormal swap後の全13pairsを再生した。
scalar cases4,1,−2、全追加monomialsもparallel PDEの完全kernelに残す。
これは有限点probeではなく全eigenvalue-support条件の必要十分な列挙である。

c-dependentではG3=z(Az²+2Bzw+Cw²)、collision divisor C/TはU外だけに支持。
Cはdegree2のa,b monomial、Tはdegree4monomial。
B_c≠0 linear、C|B_c²からC=abは不可能、C=a²又はb²。
B0のdivisibility後、actual c-shear/scaleでB=ac又はbc。
C=a²ではg2/3/4、C=b²ではg0/1/2をまず保持し、fβのglobal polynomial条件が後者g2だけを残す。
全T非零amplitudeはconnection dlogから消えても、parallel kernel/actualD3imageに全て保持した。
flatβからglobal K存在を無断に推論せず、前節の全affine K系を別確認した。
以上で13+2+4=十九のwholefirst-symbol normal formsを尽くす。

## 最弱の全support結論と二残余

十六caseは完全K1上の全D3image discriminantが恒等零。
imageの一つのbasis vectorでなく、全独立image parameterの同時多項式を展開した。
任意のfull nonlinear Poisson solutionはこのlinear必要源に入るから、指定SF枝は全空。
C=a²/g3,g4及びC=b²candidateもこの理由で排除される。

残c-dependent C=a²/g2はK1dim10、K2dim0。
P0の全七非零bracketsはactual c→c+h b変換で恒等的に保存される。
独立実装もgenerator b∂cの全raw tensor Lie derivative零を確認した。
これはlocally nilpotent linear flowなので全finite hに対するactual invarianceである。
前のB10 parameterizationのc-choiceに追加の一free fiberを与え、
既に引いたnormal-frame/c-scale fiberとは、cとbが独立でnormal nil amplitudeが非零のため独立。
従ってwholefirst-symbol B≤9、全higherは各P1上empty又は固定kerP2のtranslateなので≤10。
全像は `9+10+6−1=24`以下。旧reportのB10/upper25/false flagを改変しない。

| 未解決whole source | K1 | K2 | 旧coarse上界 |
|---|---:|---:|---:|
| DIAGONAL_-2_-2 |20|6|36|
| NIL_SPLIT |15|6|33|

両方を大きいrank2部分のdimensionだけで成功にしない。
本receiptは十七枝の最弱supportを受理し、二つのOre構造課題を明示保持する。
全nonredR m3SFの完了、全G/E、component追加はまだ結論しない。
