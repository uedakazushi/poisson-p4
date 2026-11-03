# CB原点conductor multiplicity4の独立承認

2026-11-03。対象は `AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63`。
**新しいCB実self-pair Q_b** の独立監査を入力に用いる。
CAのP(a)やCAのmultiplicity等号をCBへ転用しない。
全新規算術を再実行せず、accepted Q_bの全hashと今回のCB十生成元の対応を検証する。

Q_bの受理済み証拠は両source rank29、同じ点へ戻した接空間union rank37、
meeting edim21、G×b image rank20、generic全proper fiberが二つの被約点であること。
従って全CB completed imageは `T/(I1 intersect I2)`、normalizationは二smooth29の直和である。
二つの全射 `Bi->D=T/(I1+I2)` のfiber productとして、conductorは
`ker(B1->D) directsum ker(B2->D)`。これはnormalizerの二idempotentを掛けることで両包含を
直接証明でき、`A/conductor=D` はintrinsicな同定である。
旧R29との未検証の座標alignmentは不要。

先行のCB上界定理でconductorはdim20。新Q_bのirreducible G-saturation Zもdim20である。
Z上のgeneric長が1なら、他成分・nil supportを除いてchar0 regular openを取れるため、
稠密constructibleなG×b像のある点でedim20となる。しかしgeneric像のedimは全て21である。
従ってZのgeneric長は2以上。二旗だけから非被約性を仮定していない。

critical W/Vはt=0なのでnu=0。一方generic Q_bの二source arcは同じR値で異なるt値を持つ。
tが全quadratic algebraを生成するので、`nu=(t1-t2)^2/4 != 0`。
P0の一台全fiberとpropernessによりarcは原点へ延び、ZはW/Vと異なる第三のsupport branchとなる。
CB固有のgeneric V受理済み型と共通W型は各々smooth20でgeneric長1。
代数的local conductorのassociativityとcompletionのmultiplicity不変性から
`e(C)>=1+1+2=4`。上界4と一致する。

上界のCM20/e4 initial sourceからgr(C)へのkernelが非零なら、そのassociated primeは
unmixedなsourceのassociated primeに含まれるのでdim20・正のmultiplicityとなる。
e4同士の全射と矛盾するため、全initial equalityが成立する。
これによりactual CはCM20・edim22・e4となり、追加のembedded primeはない。
associativityの等号は三つのsupport枝、generic長(1,1,2)、第三枝の原点support multiplicity1を
固定する。excellent reduced equidimensional supportのmultiplicity1から、そのcompletionも
smooth一枝であり、formal minimal primeも正確に三つである。

全initialの二次二元と非零U³を持つ三次一元のactual liftはcomplete divisionで全残余idealを
生成し、三次式は二次idealのmaximal-ideal倍に入らない。
actual高さ2 CM残余idealのminimal生成数3からHilbert–Burch最後のrankは2、従ってtype2。
七actual linear functionを戻してmu(I)=10。
独立コードはCBの全150係数上の三liftと七linear元を用い、raw十行
18,19,28,29,60,67,72,46,95,63への遷移行列式1、外側140列の零、
CB固有cubic係数−4/243を確認する。

R29上のdepth lemmaからdepth I=21、`A=R directsum I*t` からdepth A=21。
`I subset m_R`、`nu subset m_R` よりodd productsはm_R³に入り、edim A=29+10=39。
また `m_A²=m_R*m_A` なのでm_Rはm_Aのreductionであり、finite R-rank2からe(A)=2。
coneはdim29/edim39/depth21/e2、projectiveはsmooth scalar因子を除いてdim28/edim38/depth20/e2。
元Bへの注入によるdomain性は維持され、conductorのnilpotentをCBのnilpotentと取り違えない。
generic Z上ではsquare-zero rank1 nilradicalとsmooth20 supportをopenで取り、formal smoothnessで
`C[[g20,eps]]/(eps²)` を得る。

未知binary cubicの混合係数、原点のhigher analytic normal form、actual CA/CB交叉D、
全特殊bでの全flag分類は結論しない。既受理の普遍的なregular-sequence relationsの再利用は、
CA self-pairをCBの二旗証拠へ置き換えることとは異なる。

独立証拠は `reports/AMBIENT_CB_ORIGIN_CONDUCTOR_MULTIPLICITY4_20261103_v63_REVIEW/`。
