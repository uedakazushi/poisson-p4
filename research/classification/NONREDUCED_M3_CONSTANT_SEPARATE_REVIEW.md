# 非被約R・m三constant double/simple roots：全十四源と二十二moduleの別作者確認

2026-11-07、`/root/quartic_review`。**全source幾何、独立Fraction算術及び全G非優越を受理する。**
対象はrank二plane、actualD被約既約quintic、R=2L+M又は3L、normal位数三、
generic normal cubicがdouble+simpleで、両primitive normal rootのbase degreeが零である全枝。
他primitive degrees、generic triple root、m二又は全G/Eの完了は結論しない。
C-exact三源は弱い全G上界だけで、全点既成分帰属は依然未証明である。

## 1. 不変入力と独立実装の範囲

入力はstable18file `plane_nonreduced_m3_constant_modules_manifest_20261107.json`（7124268e…）。
producerの三段の実行sourceはそれぞれROOTS36bb402d…、ORE77edcaaf…、MODULES8ca389d5…に全文固定されている。
mutable current runner又はcurrent producer jobsをinput identityに使わない。
新checker `classification_nonreduced_m3_constant_review.py` はproducer/CASをimportせず、
既独立sparse Fraction/QQ(tau)演算だけで全sourceを再構成した。
finite job `jobs/nonreduced_m3_constant_separate_review.json` gate PASS後、同一の不変reviewerを二modeで実行した。

| phase | 内容 | 時間 | 最大RSS |
|---|---|---:|---:|
| sources | 全十四＋全alpha壁、全六Euler/Q/W及び全source壁 |84.926秒|75,164KiB|
| modules | 全二十二＋全alpha壁、全affine W十八 |20.582秒|72,348KiB|

それぞれSOURCE_SEPARATE_REVIEW、MODULE_SEPARATE_REVIEWの新reportに、実行reviewer全文と全input hashを保存した。
最初のall-modeはmodule側のSchouten規約の係数二・column constant・augmented(M|-C)への対応箇所で停止した。
これはproducerの数学的欠陥ではなく独立実装の規約相違であり、失敗source全文をCONVENTION_RECEIPTへ保存した。
そのrunを全成功とは呼ばず、補正後の上記二つの正式receiptだけを採用する。
ROOTS/OREはsingle cyclic、MODULESはbilinear Schoutenすなわち二倍を明示的に合わせた。

旧六源のnegative quadratic-rowspan判定は非閉鎖の発見用履歴で、今回の証明には依存しない。
そのfalseを消さず、同じ失敗診断を再実行しない。
今回必要なのは全source包含、strictQ、全actual normalizerと全affine moduleである。

## 2. 十四first-symbolの全幾何

constant normal GL2でG3=h z²wと出来る。
double/simple rootを保つ全connectionはdiag(p,q)、twisted weightからp=dlog(f/h)。
qは全flat O(-1) connectionであり、log性を仮定しない。

f=a²bではh=a^e b^(2-e)、e零/一/二。
全qはcanonical(dlog f)/3に
`alpha dlog(a/b)+gamma d(b/a)+delta d(c/a)`を加えた源。
二重lineの高位極は全linear-potential/aで尽き、d(c/a)を落とさない。
actual c→sc+tb+uaとtensor normalizationはalphaを保ち、
delta非零、delta零/gamma非零、両零を三whole alpha lineへ移す。
全inverseはnonzero diagonalとminus shearなのでe三値との積で九whole line。
固定planeのfull symbol source upperは10/9/8、finite alpha壁では一減。
normal rootsのframeをscalar connection特殊値でも残す安全なparameter上界であり、一般orbit次元の外挿でない。

f=a³ではq=dloga+d(P2/a²)、P2 modulo C a²の全五係数。
affine base GL2とtranslationsによるdegree二potentialの全型は
零、b、b²、b²+c、bc。
Hessian rank二、rank一/kernel linear零・非零、rank零/linear零・非零を全て含む。
nonzero square-root選択は複素数上のactual invertible変換であり、zero strataを除かない。
各potential stratum dimensionは0/2/3/4/5、fullsymbol upperは5/7/8/9/10。
これで十四全源を尽くす。

## 3. 全十四の算術と八閉鎖

全三K係数のglobal D2零を元全三十monomialsで確認した。
二重line K=delta(3e−4)zw、三重line K零が唯一であり、全350 P0 JacobiとactualD零を検証。
全350×30/J1、350×9/J2、24×30 actualD3 map、全actual六shearsを再構成した。
全23のline-annihilatorはrank23でtarget h z²wを殺し、元J1へ追加した全matrixを使う。
complete kernels、全行消滅、独立性、非零minors及び全alpha根を確認した。
G3はP0/P1のcrossだけに依存し、P1²及びP0/P2はnormal order四以上なので省略した項ではない。

actualD3像零なら全m三源は空。
それ以外の全像は `B+k1+k2+6−1` 以下であり、各固定P1上のP2は空又はkerJ2のtorsorの部分集合。
全nonlinear compatibilityを解かなくても必要上界として有効である。

これでC-exact e零/一/二は全alpha・全壁で22/24/23以下、
LOG e零/一及びTRIPLE ZERO/SQUARE/RANK_TWOは全actualD3像零で空。
旧六粗上界falseは不変に残し、実成分の存在とは解釈しない。

## 4. 六全残源のEuler/Ore包含

各complete必要K1と**fullP2九**を戻し、literal二段Euler gaugeで
`Pi=Q(a,b,z,w)+partial_c wedge W(a,b,z,w)`へ移す。
全係数でc非依存を確認し、Qの全六pair×十monomialsのcomplete image/rank/basisを検証した。
全original sourceのWはW0+十八dim源へ入る。
Wzのz²、Wwのw²を除く最後の二gaugeは正当なEuler gaugeであり、unknown sourceの制限ではない。

QとWを独立に取ることは**covering superset**でよい。
元sourceの相関を捨てて存在を認める推論はせず、実解は必ずこの大きい源の解に含まれるという必要性だけを使う。
全projective Jacobiをc次数で比較するとstrict[Q,Q]零とE_y wedge[Q,W]零が正確に必要十分。
global divergenceから余分なstabilizerを仮定していない。
全Q-familyの全四strict Jacobiも独立再構成した。
全alpha壁のcomplete K1とQ imageを再度戻したので、basisのpoleをsource消失と誤認しない。

E0/E2 EXACT_Bは全alpha及び全壁でQにconstant factor partial_a/partial_bを持ち、その方向でQ不変。
一般補題Q=v wedge Z、L_vQ零から、PfQ零でactualDはc非依存。
moduleへvをwedgeするとv(D)零も従い、essential高々三となる。
既約quinticはessential二以下になれず、既八essential-three成分へ全点入る。
plane又はvの選択次元を、既全ambient成分の23へ再加算しない。

TRIPLE SQUARE_PLUS_LINEARではstrict supportがq3(q1−q2)=q0q3零。
q3零は同じconstant factor。
q3非零ではq0零/q1=q2となりQのz-column零。
元G3 lineをここで戻すことが必要であり、Wz=rho zw+sigma w²のsigmaを捨てない。
全actualD3のa²w³係数は正確に3sigmaなので、target a²z²wへのline条件がsigma零を強制する。
他の係数はこのmonomialへ寄与しない。
従って全Piのz-columnがz倍で、actualDにz因子がある。
この枝だけはW-supersetへ移した後も元leading line制約を明示的に戻している。

## 5. 残Qの全reduced supportと特殊値

TRIPLE LINEARのstrict Jacobiには−2q1²というspatial係数があるため、complex reduced supportでq1零。
残条件はq2q3零、q3(q0−2q3)零。
q3零はconstant factor、非零はq2零/q0=2q3。
actual b→b+ha、diagonal、central shearでcohomを直し、全残central二係数の零/非零から四fixed Qへ尽くす。

E1 EXACT_Bの全Qは任意alphaで
`Qab=−uaz+F, Qaw=−uzw, Qbz=uz², Qbw=u(alpha+2/3)zw`、他零。
全G3=3u ab z²wなのでu非零を正当にnormalize出来る。
alpha=1/3で旧basisにpoleが出ても、元Q imageは同じ四dim空間である。
central zwだけをshearで除き、残z²/w²の全零・非零から四whole alpha line。

LOG e二はd=3alpha+2とし、d非零・extra alpha=−8/3以外でcohom r,B,C,D,Eとcentral Fを使う。
G=d B(r−E)+(d+3)CDと置くと、独立に再構成した二Jacobi polynomialは

`Jabz=−6rD*z³/d+(4−3alpha)G*z²w/d²+6(3alpha−1)BC*zw²/d²`、
`Jabw=(3alpha−4)rD*z²w/d+(2−6alpha)G*zw²/d²+27alpha BC*w³/d²`。

従って全supportはrD零、BC零、G零。
これはgenericな係数だけを反転した結論でない。
BCの係数対6(3alpha−1),27alphaはconstant gcd一、Gの係数対4−3alpha,2−6alphaもconstant gcd一。
alpha=1/3,4/3もこの同じ三式が全てを保つ。
fullG3=a²z²w(dr−3E)の非零条件を保つと、
B非零はC=D零/E=r/r非零/alpha!=1/3のB型、
B零はconstant factor又はD型、及びalpha=−5/3/CDE全非零のcritical CD型へ尽くされる。
全B,C,Dの零補集合を含む。

scalar alpha=−2/3では元complete Q-spaceから
Ar零、BC零、AC+Br−BE零、leading E非零。
A非零はD型、B非零はscalar B型、A=B零はconstant factor。
A,B両非零はE非零に反するので欠けた交差枝は無い。
scalar Bでは新w² central方向を消さず、zw/w²全零・非零の四fixedへ戻す。

extra alpha=−8/3では新Aを保持し、strict supportは
`rD−3AC/2=0, 2B(E−r)−CD=0, BC=0, A(r+E)=0`。
A零はgeneric eight-space、A非零はE=−rとleadingからr非零。
D=3AC/(2r)、B=−3AC²/(8r²)、BC零よりC零、従ってB=D零。
これがcritical A型の全支持であり、new cohomをgeneric式へ無断に押し込めない。

## 6. actual normalizersと全central零境界

LOGを保つactual diagonalsは(a,b,c,z,w)の五nonzero scalesとtensor sa/sc。
E1 EXACT_B及びTRIPLE LINEARではsb=saという正しい追加条件を使う。
全nonzero cohomのnormalization equationsは、これらscale charactersから同時に解ける。
例えばCD型の三unit C,D,Eには、sz=σsaD、sb=saD/E、sw=σsaDC/Eを取れる。
critical AのA,rにはsz=σsb r、sw=sz²/(σsaA)を取れる。
全inverseはnonzero diagonal inverse、必要な有限square roots及びminus shear。

全a,b linear-normal shearのcentral変化は
`Xa(v)−Xb(u)−(partial_a Qab)u−(partial_b Qab)v`。
独立に代入して次の全imagesを確認した。

* E1／TRIPLE LINEAR／LOG D：span(zw)。
* LOG B：span(z²,(1−Lambda)w²)。Lambda一はalpha=−2/3でscalar Bへ。
* scalar B：span(z²)だけ。critical CD：span(zw,w²)。critical A：span(z²,zw)。

image外の全central係数を残し、residual diagonal charactersで各々零又は一へ直す。
例えばD型のz²/w²は1/σとσsa²/sw²の独立characters、
scalar Bのzw/w²は1/σとsb/(σsa)である。
zero strataも残るので、central periodを全てunitへ固定したsourceではない。
変換後Wはfull十八sourceへEuler regaugeで戻り、Q/W correlationsの元sourceもその中へ含まれる。
これで十whole lineと十二fixedの全二十二moduleを尽くす。

## 7. 全二十二affine moduleと全alpha壁

各canonical Qはstrict Poisson。WをW0+全十八で動かすと元全350はWにaffine linear。
全350×18 matrixと定数、augmented(M|-C)、complete kernels、全minors、particular又はrank矛盾を別Fractionで確認した。
**homogeneous kernelだけに置換せず、particularと全kernel方向のactualDを検証した。**
全consistent源でactualDはdegree一・二・三のproper spatial factorを持つか零。
他はaugmented rankがmatrix rankより大きく不整合。

| 全canonical源 | 結果 |
|---|---|
| E1 EXACT_B四whole line | eta零はkernel七でz²w因子、eta一は不整合 |
| LOG D四whole line | eta零はkernel七でz²w因子、eta一は不整合 |
| scalar B四fixed | 全kernel四、z²w因子 |
| TRIPLE LINEAR四fixed | eta零はkernel七でz²w因子、eta一は不整合 |
| LOG B二whole line | 通常kernel四/z²w、alpha零はkernel六又は五でw又はzw因子 |
| critical CD二fixed | kernel四、z²w因子 |
| critical A二fixed | central零はkernel七/z²因子、central非零は不整合 |

全generic minor、augmented minor、kernel/particular分母及び保存factor-quotient分母の根を保持した。
全複素rootが報告された有理alpha値で尽きることをQQ[tau] exact divisionで確認し、各壁を元全matrixで検証。
LOG Bのalpha=−2/3は全四scalar Bへ、alpha=1/3はleading零へ明示的に戻す。
他壁はfresh affine matrixのfactor又は不整合で閉じる。
選んだ一点のfactorizationをwhole lineへ外挿していない。

## 8. 正確な全G結論と残余

八initial源のG上界／空性、三Ore源のessential三／factor、残二十二moduleのfactor／不整合を合わせ、
**非被約R・m三generic21・両primitive roots constantの全枝は空又は像≤24**。
有限個のconstructible像の閉包にもこの上界が及ぶので、未知25–27優越componentの一般元はこの枝にない。

constant-factor源だけは全pointwiseに既essential-three八componentへ入る。
C-exact三源は上界22/24/23を使うだけで、全E帰属は未閉。
他primitive row degrees、generic triple normal cubic、m二/reduced R/他surface及び全G/Eは別義務。
全component数又はglobal台帳をこのreviewerが変更していない。
