# 全本質三変数既約五次の八曲線有限化：別作者pure-proof review

**独立PASS。** 八つの**具体**sourceと尖点型の独立PASSは
`ESSENTIAL_THREE_QUINTIC_EIGHT_SOURCE_SEPARATE_REVIEW.md`で既に分離した。
本稿はさらに、対象となる**任意の**curveが八型に入る必要分類を監査する。
Namba本文、FZの未確認の分類文、Moeのsecondaryな一覧を尽くしの根拠にはしない。
新component数、全P4分類の完了は本稿の主張ではない。

## 1. 仮定と全分岐

既受理canonical essential3/Saito theoremから、Fは複素絶対既約plane quintic、
有理・全特異点unibranch・全mu=tau=十二・free exponents(2,2)である。
算術種数六だから全delta=六。以下で使う外部幾何はBézout、正規化のgenus公式、
Riemann–Hurwitz、characteristic零のPlücker公式とbidualityであり、曲線一覧の引用ではない。

重複度五の点があればFはbinary quinticとなりC上可約なので、最大重複度は二、三、四。
重複度二のunibranch germはWeierstrassの平方完成でA_(2k)、delta=kとなる。
重複度三の接線接触はBézoutから四又は五で、coprimeな一Puiseux pair
(3,4)又は(3,5)、すなわちE6/E8、delta三/四である。
重複度四は(4,5)、delta六となる。

## 2. 最大重複度四：全十八係数の直接比較

p=[0:0:1]、接線y=0とし、全F=y4*z+f5(x,y)である。
x5係数は既約性から非零。x->x+lambda*yでx4*y、zのlinear shearでxy4,y5を消すと
F=y4*z+x5+C*x3*y2+D*x2*y3となる。零非零C,Dをともに保持する。

全二次syzygy (A,B,H)をz次数で展開すると、最高次から順に
b2=a2=0、B1=-h2*y/4、A1=0、C*h2=0を得る。
続いてA0=k*y2、B0=3h2D*x2/16+p*x*y+q*y2、
H1=-4p*x-4q*y、H0=-(3Ck+3Dq)*x2-2Dk*x*yとなり、残る二条件は

    5k+2Cp+9h2D2/16=0,    3Dp+2Cq=0.

全十八unknownを使った三角比較であり、残るrowは上式から恒等零である。
C非零ならh2零、kだけが自由で二次核一次元となり、free(2,2)でない。
C零D零なら(0,y,-4z)のx,y,z倍の三次元であり、これもfree(2,2)でない。
C零D非零だけが二次核二本を持ち、producerのU_D,V_DとSaito determinantに一致する。
Dの三乗根によるactual diagonal変換で全非零Dは一つの(4,5)代表となる。

## 3. 最大重複度三：二尖点を用いる全linear series

重複度三尖点pのほかに尖点qがある。重複度三の点二つは結線との交点六で不可能なのでqは二重度二。
結線はどちらの接線にもなれず、全交点は3p+2qである。
正規化でp=s0、q=t0、targetでp=[0:0:1],q=[1:0:0]とすると、
全sectionはtarget shearとscaleだけで

    (s5+A*s4*t : s3*t2 : t5+B*s*t4+C*s2*t3)

へ入る。三sectionはbasepoint-free、線形独立であり、degree五が素数なので
非直線imageはquintic、mapはbirational。
p以外、q以外のcuspはr=s/t!=0で2r+A=0、3+2Br+Cr2=0を満たすため高々一つである。

A零ならE8、p以外のcuspはqだけで、delta budgetからA4。
その必要条件C零を課すとB零又は非零。
B零のmonomial curve Y5-X3Z2はEuler以外の一次log場を持ち、free(2,2)に反する。
B非零はactual diagonal変換でB一となり、E8+A4一型。

A非零ならE6。qを残る最大deltaのcuspとして選び直せばqはA6又はA4でありC零となる。
これはC非零のsourceを捨てる操作ではない。
qの局所奇数五係数は1+AB。
AB=-1ならqのdelta三以上とpの三でgenus budgetを使い切り、E6+A6。
それ以外はqがA4なので残るdelta一のcuspが必要であり、その存在条件がAB=三。
従ってE6+A4+A2となる。A非零のsource diagonalでA一へ戻せ、ABは保たれる。
これで全最大重複度三が三型へ尽きる。

## 4. 最大重複度二、一尖点と二尖点

一尖点ならA12。接線接触は四で、残りのsmooth交点をt0へ移すと
全三sectionは

    X=s4*t,
    Y=s2*t3+C*s3*t2+A*s5,
    Z=t5+B*s*t4+D*s3*t2+E*s5.

他の項はX,Yのtarget shearで消した項だけである。
奇数五、七、九、十一の消滅条件を順に解くと
B=2C、D=2A-C3、E=-AC2、A2*C=0。
この全二parameter残式の有限恒等式は独立checkerでも検算する。
A零では全sectionがtを因子に持つためbasepointとなる。
従ってA非零、C零。source sの三乗scaleで全A非零を-1へ直せ、A12一型となる。

二尖点のdelta分配は(5,1),(4,2),(3,3)の三つだけ。
高いdeltaのpは接線接触四。結線が接線でない場合、残るsmooth交点をs=tへ置き、
全sectionは

    X=s5+A*s4*t+B*s2*t3,
    Y=s2*t2*(s-t),
    Z=t5+C*s*t4+D*s3*t2.

pがA6以上ならB零、A非零。qがA4以上ならD零である。
A8+A4ではT=AC+2A+1=0、K=C-AとしてK(C+3)=0。
K零はA=C=-1となり三section共通因子s-tのbasepoint枝。
残る全枝はA一,C=-3であり、A8+A4一型となる。
A6+A6ではT=0とAC+2C+1=0から同じbasepoint枝だけで、空。

A10+A2ではDを自由に保持する。奇数五、七を消すと
C=-2-1/A、D=(A+1)2(A-1)/A3。
次の奇数九係数は2(A+1)2/A4となり、A=-1のbasepoint枝しか残らない。
結線がqのA2接線となる残chartも保持する。
その全Y=s2*t3、X=s5+A*s4*t、Z=t5+C*s*t4+D*s3*t2でA,D非零。
奇数五、七はC=-1/A、D=-1/A3を強制するが、奇数九係数2CD=2/A4非零で矛盾する。
結線がpの接線というchartはBézoutで不可能。
したがって二尖点全sourceはA8+A4一型だけであり、tangent/残交点の補集合を落としていない。

## 5. 三尖点以上の全数値候補とCremona排除

二重度二cuspからのprojectionはP1のdegree三mapで、全ramification四。
他の各cuspが一以上を与えるのでcusp数は五以下。
五つならdelta分配(2,1,1,1,1)であり、A4を中心に取るとその接線接触四又は五により
中心自体にもramificationがある。四つの他cuspsと合わせて五以上となり矛盾する。
従って数は三又は四だけで、全delta六の分配は
A8+2A2、A6+A4+A2、3A4、A6+3A2、2A4+2A2である。

A8+2A2ではA8のpを[0:0:1]、接線x0とし、
psi1=(x2:xy:y2-xz)を使う。唯一basepointはp、normalization共通order六なので
degree十から六を引いてimageはquartic。psi1の合成はx3*(x:y:z)であり、x非零の開で実birational。
pでy=r2,x=r4+a6r6+a8r8+a9r9+...、a9非零からimageはA2となる。
他の二A2はBézoutによりx非零の開にあり保たれる。
接線x0の残smooth点はimageでsmooth、接線X0はpのimage cuspも通る。
そのcuspからのdegree二projectionは他の二cuspsとsmooth点でramifyし、総数二に反する。

A6+A4+A2ではp=A6=[0:0:1]、q=A4=[1:0:0]、p接線x0、結線y0とし
psi2=(xy:y2:xz)を取る。base共通orderはp四、q二、他零でimage quartic。
明示inverse(X2:XY:YZ)との合成はxy2*(x:y:z)で、xy非零の開で実birational。
pはx/yとx/y2のorder二/三からA2、qもyとz/yからA2となる。
第三A2はexceptional二直線に乗れないので保たれる。
結線の残smooth点はimageでsmooth、接線Y0はqのimage cuspを通る。
再びそのcuspからdegree二projectionのramification三で矛盾する。
いずれのargumentも四次imageの三A2がgenus三を使い切るので、imageの衝突や余分な枝を隠さない。

## 6. dual curveと残る三A4・四尖点

2A4+2A2のdual次数はPlückerで20-(12+4)=四である。
A4のlocal x=r2,y=c4r4+c5r5+...、c5非零から、dual slope/interceptは
c4非零なら再びA4(delta二)、c4零なら(3,5)(delta四)となる。
二つのA4がdualにdelta少なくとも四を与え、同じdual点なら枝交点分も加わる。
quartic genus三に反する。この排除はFZの「dual cubic」という誤記を用いない。

3A4では三image点はBézoutから非共線。
正規化の三点を0,∞,1、targetを三座標点とすると、全sectionは
X=t2(s-t)2(As+t)、Y=s2(s-t)2(s+Bt)、Z=s2t2(Cs+Dt)、C+D非零。
各点の二次/三次jet vectorのdet零というA4必要条件は
BC+(2B-1)D=0、AD+(2A-1)C=0、AB+2A+2B+3=0。
初めの二式のdet零と最後の式からAB=-1、A+B=-1。
D零はC零を強制して不可能なのでD一へscaleするとC=-A-2。
Aの二根はp,q交換により同値であり、Z-scaleは
(-A-2)(-B-2)=一から正しく戻る。
このcanonical seriesはbasepoint-freeであり、三点のdelta各二以上が全六を使い切るので
正確に3A4である。独立source receiptの有理代表は直接3A4と認識済みなので同じ一型へ入る。

A6+3A2のdualもdegree四。
元A6はprojective tangent接触四なのでlocal y=c4*x2+c6*x3+c7*x^(7/2)+...、c4c7非零。
Legendre slopeがxの単位倍から始まるためinverseを取り、interceptの微分=-xを使うと
dualの最初のfractional exponentも7/2、すなわちA6となる。
quartic genus三を使い切るのでdualは一尖点A6 quartic。
その全sectionはX=s4,Y=s2(t2+b st),Z=t4+d st3+f s3t。
A6条件d=2bの後、source t->t+b*s/2とtarget shearでb=d零へ直り、
(s4:s2t2:t4+f' s3t)を得る。
f'零はconicへのdegree二mapでquarticでなく、f'非零はscaleで-1にできる。
したがって一つのPGL3型。Q=(y2-xz)2-yz3の具体parametrizationと全gradient代入を
独立に検算するとdualは公開A6+3A2 representativeそのものである。
bidualityにより元の全quinticも一型となる。

## 7. 結論とreceipt

全必要curveは(4,5)、A12、E8+A4、E6+A6、E6+A4+A2、A8+A4、3A4、A6+3A2の八型。
全nonzero parameterのscale逆と全basepoint/tangent補集合を上で保持した。
各具体curveの存在・複素絶対既約・free/Saito/full source・正確なcusp型は
先行の独立explicit-source receiptで確認済みであり、この必要分類との循環はない。
accepted canonical essential3 theoremと全磁気inverse、curve変換のPGL5への延長により、
対象全Poisson sourceは八つのactual23軌道の和へ入る。
その閉包の全成分性・prior相違は別ambient receiptを使い、ここで新件数を変更しない。

checker `code/classification_essential_quintic_curve_proof_review.py`、SHA256
`3325dfc8f671dd0eeb3d50cea1b4bb756788114fc02d6cdb364abaa1bf1c0b15`。
job `research/classification/jobs/essential_three_quintic_curve_proof_review.json`、SHA256
`3dee10ac13a398c014bb62841505521c9f6447a17c79f727346e1773f3654d18`。
gate PASS後の有限恒等式再現は0.083318553秒、最大RSS11904KiBでPASS。
report `reports/CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_CURVE_PROOF_REVIEW_20261107.json`、SHA256
`7d84273bf719c4d41e99ae17be45f603f2d54dfad573cbdf55c85e75ed3ca68d`。
この小さい計算だけを全curve分類と呼ばず、上の全source幾何と併せて受理する。
explicit-source checker/reportの旧bytesは全て凍結のまま、別fileで実行した。
依存manifestは `inputs/essential_three_quintic_curve_proof_review_manifest_20261107.json`。
