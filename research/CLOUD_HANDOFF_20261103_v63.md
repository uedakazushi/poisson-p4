# v63 Cloud 継続記録 — 2026年11月3日

研究継続中。これは停止指示ではない。全未解決課題を完了したとは主張しない。

開始SHAは `4c3f623fb698f9adc6077015b25e2325b11072af`。
v62保存は `e187b40de68bce534391f406f616940cddd045a7`、
`9e82f3c5120380028990cdb58b07318866c21f4b`、
`66a80ab53dfcf0ff842c67140320d8da0226039d`。全て通常push済み。
このv63保存のSHAはGit履歴で識別する。

## v63 Y: 実接触長4・被約交叉の全gluing・半正規化

直前保存 main `57457c3acc785101d6a66e23aa6f8b95635458a5` は通常push済み。以下は過去時点を含む。

generic B18での実N∩Rを正確に決定した。全N sourceの27 IFTと二次Weierstrassを二変数三次まで全150係数で構成し、独立なLie-log計算でも全113 graph cubicの零を確認した。零一次接触h4と合わせ、実接触kも正確に4である。群十九因子を含む共通商はC[[g19,T]]/(T4)となる。

被約交叉Dredの全形式環は、既知のactual N全完成環と正則R完成環を、この明示的length4商でglueするfiber productとして決定した。actual N37生成元の中心化座標の非零値はu16=13T,u22=−25T²,u65=−12T,u112=−13Tのみで、有限gluing mapも保存した。被約Dのcone数値はdim28/edim37/depth20/e2、projectiveは27/36/depth19/e2。全非被約Dのdepthはcone19〜20のまま未決である。

被約Dの正確な半正規化もC[[g19,x9,T]]/(T*x1,...,T*x9)と決定した。二つの正則normalizerをreduced G19上でglueする有限subintegralな環であり、square/cube判定でseminormal。cone28/29/depth20/e1、projective27/28/depth19/e1。ここでのconductorは半正規化からnormalizerへのものに限り、元Dredのconductorを計算済みとはしない。

T_Dred=T_N37が確定したため、全Dとのcotangent差は正確に二次元となる。物理係数の基底classは[u46+u95],[u63]。各classには形式的nilpotent liftがあるが、生の線形函数自体をnilpotentとは主張しない。さらに全Dの真の接錐の被約台は一つの線形28平面、接錐環のnilradicalの六乗は零と証明した。これは元Dのnil indexや全initial idealの決定とは別である。

W16のSpin4切片を二形式U(2,1),V(1,2)として同定し、その既約混合判別式D(8,8)を構成した。quintic判別式からRes(V)²を除く正確な恒等式とnormalizer不変性を証明し、D(p0)上へ斉次次数16の正則有理不変式として降下する。極次数、24次多項式lift、NEW25との対応は未解決。NEW25方程式次数下界22と実方程式未決の状態は維持する。

非平方はadFBCJ全開でのQQ product41が既証。full cubicには実9×43判定が残り、全四次354→319の同span縮約と286×238の十分条件を独立検証した。これらは全三次57・四次全射の証明ではない。次はこの実残余と元Dのnil・exact depth、一般betaのδ_E∈J0判定を続ける。総下界58/rank4下界52、総上界未定。五つの凍結原稿・PDFは変更しない。

元四quadのZW商について、全aC開で20×35 product階数は正確に19と証明した。二つのquadに共通linear factorがあり、既知一本のsyzygyが常に残る。generic CIとcolon exact sequenceを合わせると、light四次quadratic productはgeneric146・全点上界146、weighted版はgeneric152・全点上界152となる。全点での下界146は未決であり、147を目指す前提は誤りと訂正した。QQ恒等式の保存117項には零一項があり、非零116項である。

さらにvertexの実binary resultantが非零となる明示的R9開では、raw四quadのCI4とpointwise146を証明した。R9は未知degree9多項式として保持し、全域で反転しない。残るR9=0と、四次での追加三次・四次関係の全射は別に続ける。

Spin4混合判別式は、一つの既受理NEW25点で厳密に非零だった。dense DAG、全β、QQ二次体frame、全170 circuit node、Sylvester判別式を独立再構築したため、この候補とp0冪で極を払ったどのnumeratorもNEW25方程式から除外される。Spin4構成自体は有効で、別の方程式候補を継続する。

## v63 X: adFBCJ全域のQQ product41とB18接触

直前保存 main `e92227d0afe6e27a778cd67d053321209b673761` は通常push済み。以下は過去時点を含む。

非平方の全adFBCJ開（P=ad表示ではaPFB0CJ開）で、元の42×84 product行列の階数が正確に41であることを証明した。最後のpair条件を全QQの有限商・subresultant・inverse恒等式だけで解決し、既知のgcd例外処理とbinary-D縮約から全raw24全射へ進んだ。元の六二次式への17 kernel方向と全known-unit座標変換も独立に再構成した。元productのimageはrank41 subbundle、cokernelはlocally free rank43である。有限体のunitをnonproper開へ移した議論や、未知因子の反転ではない。

この結論はd0/J0とより広いaFBC全開、full cubic57・四次全射性・ambient除去を含まない。product41とheavy7だけから出る三次下界は48であり、残る実9方向の特殊化を別に調べる必要がある。追加Δproduct・Δheavy chartでの9×43 bordered-minor circuitと、四次354→319行の同span縮約を保存した。部分的なrank結果を全三次・四次の証明へ昇格しない。

B18原点では、実R曲線の接方向がN sourceの微分像へ入ることを確認し、実N曲線の曲率を差し引いて全像イデアルの接触次数k≥3を得た。正規化sourceへの二次持上げは全kernel速度を許しても不可能であり、この失敗と像イデアルの接触は区別した。全81関係とtrace/discriminantの実際の二次曲率から、零一次N関係の接触は正確にh=4と証明した。kは3または4であり、対応する被約横断edimは19または18、線形nil次元は1または2となる。

primitive W16をSpin5のgamma-traceless vector-spinorへ全係数で同定し、二つの二次adjoint共変量を構成した。4×4 pencilの行列式五係数は全8次不変式の基底となる。これは新25成分の方程式同定ではなく、その次数下界22・実方程式未決という状態は維持する。

一般beta五点尽くし、全成分上界、B18原点の全相対環・exact depthは引き続き未解決。総下界58/rank4下界52、総上界未定。既存118/16/91/14/17頁の原稿とPDFは変更しない。

## v63 W: 第二embedded成分・全三次gate・有限例外のQQ除去

直前保存 main `7ea5cf639a99bd0d57a6d2b5d5aaad175c1a2958` は通常push済み。以下は過去時点を含む。

Q_b側のactual CA/CB交叉Dの全完成環を、全150係数のQQ(b)三次syzygy計算と独立Lie展開で決定した。非空unit openと群飽和上でC[[s20,z7,f,rho]]/(fz7,f³,f rho²)となる。nil moduleは滑らかなS20上free rank4、二乗rank2・三乗零で、embedded primary quotientの長さは6。既知P_a側のsquare-zero rank1 nil moduleとは異なる。両stratumの射影19次元閉包は相異なるglobal embedded associated componentである。全associated primeやb0/B18原点を尽くしたとは主張しない。

新25次元成分のfirstjet hypersurfaceの方程式次数下界を22へ改善した。全20次不変空間のcharacter次元36と、実dense source36点でのQQ非零36minorにより全20次制限の単射性を証明した。独立な整数word評価とBareiss determinantも一致した。実方程式と次数22の存在は未解決。

非平方では高rank比例を全known-unit開で除外し、残る有限gcd例外coverの全raw QQ全射性を証明した。有限Z_(32003) modelで全localizerのunit性を先に確定し、全672×1568作用の全射からNakayamaでQQへ移す。nonproper開のmodp計算をそのまま移す議論は使わない。別の補助5×21では第四unit pivotも厳密確認した。pair gateと全域のrank条件、元product41/full57/四次被覆は引き続き未決。

個別CA/CBの全形式型とseminormalizationは既証。B18原点のactual D全型、一般beta五点尽くし、全成分上界は別に継続する。総下界58/rank4下界52、総上界未定。凍結済み118/16/91/14ページ原稿は変更しない。

generic B18ではOLD共通29底を使って二成分の実150×39接埋込みを再構成し、T_CA=T_CB=T_Dの同じ39平面であると独立検証した。actual Dはcone次元28・edim39・多重度2でdepth≤20、projectiveでは27・38・2・depth≤19。さらに全schemeの強い群不変性とformal IFTで十九群因子を分離し、depthはcone19〜20、projective18〜19と絞った。横断環のdim9/edim20/e2とdepth0〜1は確定し、embedded vertex/socleが最後のdepthの差となる。個別完成環の抽象同型から相対配置を決めたものではなく、D_redの接平面やnilの可視性も未決である。

非平方の別のstructural証明により、六linear entryのrankは全aPJ開で正確に4と分かった。二つの残余係数は四linear変数だけへ依存し、binary Dのx0²係数はknown unit。従って補助五二次rank5はpair gateから従う。DB0のrank5自体も独立に証明済み。有限gcd例外の処理と合わせ、全raw全射の十分条件はA=x0=W0上の二quadric独立性という単一pair gateに縮約された。この最後のgateはまだ未証明である。

次に行う実課題は、(1) pair gateの構造的QQ証明、(2) B18原点の相対全idealとsocle、(3) 実42変数δ_E∈J0の判定である。pair gateのcofactor degree3/4探索、QQ scalar GB60秒、補助rank GB300秒はいずれも証明を与えず、失敗を非所属へ昇格しない。大きいQQ GBの同じ入力での反復は避け、既存の小さい構造式を利用する。全分類の完了条件・上界は引き続き未達である。

さらにDVR曲線とのfiber productを用いる新しい環論的議論で、B18原点のnilradicalがm²に含まれないことを証明した。被約横断環はedim≤19だがactualは20なので、nilには少なくとも一つ一次cotangentに可視なclassがある。被約接空間を単に二成分の接空間の和と仮定していない。正確な被約edim・特定nil座標・全環は未決。個別形式型と二つの交叉型を新規17頁日本語原稿 papers/ca_cb_local_normalization_ja.pdf に独立保存した。


## v63 V: 全形式標準形・seminormalization・全99 kernel

直前保存 main `2ce106c7425a114d98efef89704d9c07162e0852` は通常push済み。以下は過去時点を含む。

generic B18代表点で、個別被約CA/CBの全形式標準形を決定した。両者は新しいformal座標でA=R⊕I*t、R=C[[g19,z6,ν,w,r,s]]、I=(z6,ν−s(s−r),wr,ws,sr²)、t²=νとなる。三軸conductorのnil moduleと函数の剛性、G19-equivariant finite coverのproduct、最小Cohen表示、intrinsic pushoutを分けて証明した。全kernelは44 syzygy＋55積関係の99式。作者側のKoszul/Hilbert–Burch証明と独立監査に加え、Cloud M2 1.24.11のdefault eliminationでも全source kernelとの一致を確認した。

全formal germの特異・非正規・非CM・非seminormal locusも、三枝のreduced conductorに正確に一致すると証明した。seminormalization後の前三者はself枝だけに縮む。normalizer fiberの全局所stratumも記録した。

正確なseminormalizationもA_sn=A[wt,rt]=R⊕J*t、J=(z6,ν−s(s−r),w,r)と決定した。二新元は元Aのsquare/cube条件を満たし、domain R/Jとν=s²非零を使ってnormalization B全体のsquare/cube判定を証明した。A⊊A_sn⊊Bで、self-pairの二点gluingは残る。全81関係とnormalizer conductorのreduced node coverも保存した。

これらは個別完成環のformal同型であり、収束holomorphic同型や元のCA/CB相対配置を与えない。従ってgeneric B18原点のactual D全型はなお未解決。全Poisson分類、総数上界、一般beta五点尽くしも継続する。

Q_b側CBのdouble meetingと全81 fiber-product型も独立検証した。実交叉Dは一般に非被約だが、P_a側と相対接空間が異なるため、以前の全D型をラベル交換しない。alternate branch∩CAには閉長4・relative free4とisotropic quotientを全150係数の独立series演算でも証明した。98秒のsymbolic workerは再実行せず、完成済み入力を独立な短い計算へ移した。actual Dの全idealには最後のambient cubic条件が残る。common R29へのétale性だけで余分なambient関係を消せるという候補は撤回し、正確なcolon criterionと未検証点を保存した。

非平方の四quadricでは、単元a²Jだけによるexact row reductionと5×21の小さい条件へ移す証明を保存した。generic regular sequenceはQQで確認したが、全unit開でのgcd/constant-rank条件は未決。大きいQQ GBを繰り返さず、低次係数の小さい厳密問題を続ける。総下界58/rank4下界52、総上界未定。

## v63 U: 両原点conductorとactual Dのembedded成分

直前保存 main `b1f02523f7485255f5670c28cee4905367b8aa80` は通常push済み。以下は過去時点を含む。

generic B18原点の個別被約CA/CBで、conductor quotientの全初期idealと十生成元を証明した。七actual constraintのIFTとintrinsic G19作用からdim20/edim22/e≤4を得て、critical W/Vの長1ずつと実自己旗成分の長≥2をmultiplicity associativityで合わせた。e=4となり、unmixed CM初期模型からの非零kernelが排除される。両conductorはCM20/type2、actual Iは十行で生成される。個別CA/CB imageはcone29/edim39/depth21/e2、projective28/38/depth20。高次analytic normal formやB18原点actual Dの全環はまだ決めていない。

初期binary cubicも両方で正確に認識した。CAは(12/125)U(U−V/3)²、CBは(−4/243)U(U−9V)²であり、未知だった三次係数は解消した。これは全graded initialの認識で、高次formal normal formは別課題である。

新しい別族Q_bでCBの実自己旗二枝も全150係数と全微分で独立構成した。全proper旗被覆を両族で完了し、全非零parameterでP_aはCA二旗・CB/N一旗、Q_bはCB二旗・CA/N一旗となる。一旗側のCA/CBはsmooth29、Nはsmooth28。両族を取り違えず、成分数は増やさない。

P_a側CAの二枝のmeetingを真二次差からC[[s20,rho]]/(rho²)と決定し、全81関係で全CA完成環を提示した。generic unit近傍でcone29/edim38/depth21/e2、projective28/37/depth20/e2。二枝接空間の和37だけをwhole embedding dimensionとする以前の誤りは使わない。

さらにactual CA/CB交叉Dを同じ開で全次数で決定した。QQ(a)の全150真二次差でT_CB⊂T_CAを証明し、全81 kernelのCB制限とS全体の接空間条件から全idealがf*(f,u8)となる。D_hat=C[[s20,u8,f]]/(f²,fu8)、nil module=C[[s20]]f、associated primesはminimal Nとembedded S20の二つ。projective Dはdim27/edim28/depth19/e1。Sの閉包はglobal actual Dの新しいembedded associated成分(projective19)で、既知minimal R19とは異なる。B18原点や全primary分解へは広げない。

一般beta残余の実scalar所属delta_E∈J0、QQ product41/全57/四次被覆は継続中。raw e0のGF certificateは原wrapperの誤りと独立全24式確認を分離して保存した。三既存fieldのcofactor復元はQQ候補に至らず停止し、有限体結果をQQへ外挿しない。総下界58/rank4下界52、総上界未定。Cloud環境の短い再起動後も作業fileは保存され、セットアップ全体はやり直していない。

## v63 T: 実CA自己旗枝・一般点の全枝、次数20とactual scalar target

直前保存 main `3233e823d0f23551ca52bb319d5354bf8e719602` は通常push済み。以下は過去時点を含む。

CAで二つの異なる旗が同じ実Poisson構造を与える全次数の多項式直線P(a)を構成し、全150係数と全strict Jacobiを独立検証した。有限三次・四次候補から実際の枝へ進んだが、rankは全aで2である。全proper fiberと両source immersionにより、一般aの被約CAは正確に二つのsmooth29形式枝を持つ。meetingのedim21からprojective depth≤21で非CM。接空間の和37をwhole CAのedimと取り違えず、clean meetingや原点full idealは未決とする。

全直線は第一CA旗を通じて実Nに入り、明示的CB表示も得た。第二CA旗はN旗ではない。全N sourceの閉subbundleと全fiberを監査し、一般aではNがprojective27で正則と証明した。指定CB liftも全a≠0でimmersiveだが、CB全fiber型はまだ決めていない。N所属からR19所属は従わず、CA smooth R minus Nの既証定理と矛盾しない。

この族にはa≠0のみを反転する一様なPoisson Cremona(2,2)があり、有理中心は正確にC(H/Z³,ZW)。一般葉の閉包は通常次数7、4a²+1=0では5の有理coneとなる。全150係数の実sourceから既存LPB39・Log1112_28にも全直線が入る。成分数は増やさない。

新25次元成分のfirstjet超曲面について、完全な18次不変式空間への制限が単射とexact QQ九次minorで証明した。独立整数word計算も一致し、定義方程式の次数下界は20へ改善した。方程式自体は未同定。

一般betaの未知残余にはactual42上の新EulerRed determinant delta_Eを構成し、次数は正確に28と証明した。B空性はdelta_E∈J0という2838個の次数≤4生成元へのscalar所属と同値であり、polynomial cofactorには次数≥24が必要。所属そのものは未取得。raw四quadricの小さい有限体GBと47個のexact QQ seed所属、bounded QQ GB失敗はfrontierとして保存し、QQ product41/全57/四次被覆へ昇格しない。総下界58/rank4下界52、総上界未定。

## v63 S: V/W境界の全完成型と固定Darboux変換

直前保存 main `cd5002934c19f09a0016b1db55ff729796368792` は通常push済み。以下は過去時点を含む。

CAとCBのgeneric V境界で、個別被約成分の全完成環を独立監査して決定した。いずれもC[[g19,s8,t²,s_i t,t³]]⊂C[[g19,s8,t]]というprojective Ac8型で、dim28、edim37、depth20、多重度2。全81 kernelは正則列のKoszul関係と積関係で尽くされ、M2未実行を証拠に数えない。CAはV1稠密軌道、CBはmodulus Vを保つQQ(V)計算であり、CB全Vへは外挿しない。normalizationは非零actual Iが分数体を回復することを証明してから同定した。

共通のactual W曲線でも七一次と二WU/WV classの九次minorからgeneric全完成型を得た。さらにr^-1 diag(1,1,r,1,1)_*P_W=P_(W/r)を全150係数で証明し、非零Wの単一projective軌道によって全W≠0へ移送した。CA/CBは上のAc8型、NはC[[g19,s7,t²,s_i t,t³]]型でdim27/edim35/depth20/e2。normalizationはいずれもelementary subintegralである。W0のCA/CB全idealとactual CA/CB交叉Dは未決のまま。

同じW曲線にはWを反転しない固定polynomial Darboux自己同型と逆式を構成した。正確なCremona bidegreeは(5,5)。W0では二つのpolynomial Casimirにより全affine葉がA2となり、全有限値の葉閉包は次数5の有理coneである。B18の全特殊点へは広げない。

CA原点のordinary自己旗三次条件は全150式を保ってGmと決定し、有限三次witnessも確認した。実formal/algebraic枝の存在ではない。高次条件を続ける。非平方QQ product41と一般beta残余unitは未証明。homogeneous modular GBのlate16prime、early3primeと同じ93bit比較は係数復元のfrontierとして保存し、QQ証明へ昇格しない。総下界58/rank4下界52、総上界未定。

## v63 R: 忠実な単一scalar criterion、CA/CB odd欠損、厳密回路34

直前保存 main `5c5ec7a44f93b6291ee0bc0887bdcced033a90a4` は通常push済み。以下は過去時点の記録を含む。

CA/CBのgeneric B18実odd idealは、判別式nuを含まないことまで独立証明した。全actual generatorsの零一次二次spanはWU/WVのみでUVが欠け、nu liftの非零UV係数と矛盾するためである。非線形critical消去を保った三次計算も独立監査し、IのJ/mJへの像は十次元、mu(I)≥10と決定した。生成数上界と全I型・exact depth・実Dは未決。R19の二旗を個別CA/CB conductorと取り違える議論は使用しない。

一般affine残余は、一つのgeneric定数aによるmoment Schur局所化へ移した。この局所化はBに対して忠実平坦なので未知成分・nilを失わず、単一scalar idealのunit性は全未知残余Bの空性と厳密に同値である。重根も含むquartic normによりdeltaのa次数≤532を証明した。元全670式degree≤10、またはNormと134代入式からなる135生成元degree≤34という具体的な有限targetを保存した。unitそのものは未取得であり、一般五点尽くしは未解決。

同じ残余空性は、さらに具体的なvector所属 delta e0∈image Phi、またはscalar所属 delta²∈Jgと厳密に同値である。独立staircase検証によりformal target次数は7/14、135生成元上界は30へ改善した。実所属は未取得であり、cofactor次数上界を主張しない。以下の34は直前の有効な弱い次数上界である。

この構成に先立ち、全3350 scalar entryの一般二次rank179/total224と全45一次scalar spanを証明し、完全な十五Schur coverも構成した。scalar Fittingの等号をAのannihilationと混同しない。無限遠では全134 midpoint恒等式からaxis冪等元を作り、D(w)を保って正確に分離した。さらに全61変数の普遍degree3恒等式wG6^h=sum Ki^h Fi^hでI134⊂I135⊂I134:wを得た。全w-torsionの除去を主張しない。

非平方の新QQ scalar unit a P²Jを975項cofactorで厳密証明し、次の23×58行列を全6244-node回路で構成した。数学的に構成済みpivotは34、展開済みnative cacheは33である。67000-bit定数や巨大outer展開を避けるための厳密な保存形式であり、未知factorを反転しない。全点product rank下界37は別の既証定理で、QQ product41は未証明。homogeneous module GBの小さい有限体pilotは成功したがQQ復元・全検証は別に継続中。下界58/rank4下界52、上界未定。

## v63 Q: generic B18のN全完成環、大域五次表示、厳密有限体証明書

直前保存 main `0e87048ae6eb55900631e4581b0d26b9895070ab` は通常push済み。以下は過去時点の記録を含む。

generic B18軌道で、個別被約target Nの全完成環を決定した。全actual odd generatorsの七一次＋二二次classを元150係数とNeumann二次逆で独立監査し、J=I+(nu)の全九minimal classを得た。NakayamaによりI=J=(nu,s6,wv,v²)が全次数で成立する。projective完成環は正則環C[[g18,w,v,s6,t]]の閉部分環 C[[g18,w,v,s6,t²,s_i t,t³,wvt,v²t]]。全81関係のkernel等号を別のM2でも確認した。dim27・edim36・depth19・多重度2、domainでembedded primeはない。

Nのnormalizerは正則で射はsubintegral、seminormalizationもこのnormalizerとなる。conductor quotientはC[[g19,w,v]]/(wv,v²)で、そこにあるembedded nilとN domain自体を区別する。このgeneric完成germのSing/Nonnormal/NonCM supportはconductorの滑らかな被約台である。Nの全特殊B18点やglobal bad locus、CA/CB個別の全型、実CA/CB交叉Dは未解決である。

一般affine全134式を、normalR42底上の大域monogenic五次と134余り、5×670 moduleへ移した。S1(n),S2(n),H(n)は固定有限circuitであって自由parameterではない。自由formal45変数へ一時的に持ち上げた全670行は二次頭rank670で、定数cofactorの純一次vector spanは零と証明した。ただし非線形basechange後の実42変数障害に移さない。全beta有限性や一般五点尽くしはなお未取得。

Qa8/9の実T38 sourceを明示的四次体と全150係数の多項式curveで構成した。固有frame・全unit norm・全350Jacobiを独立監査した。固定frameでlevel−1 edgeとlevel5 seedの二つの単純経路は除外したが、P+やP2のT38全所属は未決である。

非平方について有限体で厳密証明したのは、p還元済み31-prefix残余10×53の全射性だけである。元product41への全31prefix素数整合性は未監査、QQ41も未証明。完全な183182項証明書と無効なpartial exportを区別して保存した。八素数CRTは係数の安定再構成に至らず、V27 archiveへ実行状態・全signature・再現手順を保存した。新しいQQ scalar unitとnative縮約は別に進める。下界58/rank4下界52、上界未定を維持して研究を続ける。

## v63 P: 一次像次数18、product下界37、実normal graphとN部分射影

直前保存 main `2f0d9f8` は通常push済み。以下は過去時点の記録を含む。

新成分の15次元primitive一次像を定める方程式の次数下界は18へ改善した。degree16全不変空間20次元への実像制限は、全二十QQ点と厳密20minorで単射、degree17は零。全physical正規化・全1200値・独立QQ行列式を監査した。degree18の存在・不在と実方程式は未解決。

非平方adFBCJ開の全点product rank下界は37へ改善した。既約quadを使うGrassmann incidenceとFano全rank分類、三変数shadow8＋12＋17方向による証明である。第33pivot後のunimodular行の存在は非構成的に従うが、実行済みpivot数33・24×58は変わらない。有限体で検証済みなのは縮約10行列の全rankであり、元product rank41への全prefix単元条件は未監査。有理数right inverse探索と明確に区別する。

一般残余のcompact divided normal incidenceは、全42方向chartで四つのfiber座標を多項式に解け、P41へのscheme閉埋込みとなる。実残余はs≠0の開部分である。元134式に属する全45minorから、raw空間には一般betaで被約孤立無限遠成分もある。三次元無限遠族との共存によりraw pure3/CM-pure3方針は偽である。選択した単純無限遠だけを除く実sheaf colon分離を得たが、全w-torsionと実affine空性は未決。modp nilpotenceだけをQQへ移さない。

B18の共通29実target内でCA/CB odd一次係数rank7と独立二判別式を確認した。raw Y28をrow61固定で動かすrank7と、adjusted etaを固定するrank6を区別した。Nの共通29座標への実部分射影像はordinary pinch C[[g26,nu,b,z]]/(z²−nu*b²)となる。これは全N targetや実CA/CB交叉Dの型ではない。全odd idealの高次lift・相対位置、P2実枝、全generic尽くしを続ける。総下界58/rank4下界52、上界未定。

## v63 O: 正軸の被約smooth性、normal残余、全B18での個別非CM性

直前保存 main `9c8f974cdf70d6ba1313993a70f4e8823d90da80` は通常push済み。以下は過去時点の記録を含む。

正の重み六軸の全非零点で、新25次元被約成分は滑らかである。全Poisson横断環の多重度上界5と新成分一般長さ5を結合公式で比較し、全多重度5・被約新成分多重度1・唯一の25次元極小成分を証明した。解析的分裂を先に仮定せず全branchを数える。低次元成分・embedded構造は未決で、全Poisson環は非CMのままである。

一般affine全134式はnormal座標42個の多項式環上、実際に五元で有限生成される。全z*n cross rank126とtiny plane行列を元式から独立監査した。実残余のnormal射影は閉埋込みとなり、2796低次数recipeのorigin saturationで定義できる。さらにnormal-only三次commutatorの42Jacobianがp-unitで、42式を加えたJ0は被約originと残余の直和、normal-only冪等det D84により残余=J0+(D84)。これはbeta底上有限性や残余unitではなく、一般五点尽くしはまだ未解決である。

B18におけるCA/CB/Nの正規化source全ramification schemeを元150係数から計算・独立監査した。CA/CBのcone完成環はC[[g19,w,u,v]]/(wu,wv,uv,u³)、NはC[[g19,w,v]]/(wv,v²)。source結果を実交叉Dへ転用しない。実odd idealの判別式上への像は六線形と上の小coreで、最小生成数は10/10/8。Nの正則底は自身の28cone変数で別に構成した。

有限R加群A=R⊕I*tのCM性ならIがprincipalになる矛盾から、個々の被約CA/CB/NはB18全体で非CM・非正規・特異であると証明した。bad locusへの包含であって等号ではない。共通R内のodd idealのnu方向liftと相対位置、実Dのnilradical・Ass・depthは未解決。

J0の二つのDixon単元も全row・witnessを独立Fractionで監査し、39安全pivot・18×52を受理した。最新primitive入力は876541698 bytes、SHA039ee56e033a9d576aae8e271483f34e10dd111e87f391eb740251dc90c9489a。全coker unitやJ0欠損有限性は未取得。

非平方は元二行split right inverseと次transportにより33安全product pivot・24×58に進んだ。全product rank41、全cubic rank57、四次被覆は未証明。再現用V24 archiveは小さい入力・witness・失敗限界を保存し、巨大生成行列はhashとrecipeで再生成する。下界58/rank4下界52・上界未定を維持し、残るP2枝、一般残余、B18実交叉、非平方の全被覆を続ける。

## v63 N: 新成分の一般長さ5と全完成環、独立原稿

直前保存 main `0bc26c274760d02af900dd4f2fee9aaef2e3e966` は通常push済み。以下は過去時点の記録を含む。

新25次元成分の全Poissonスキームの一般横断長さを正確に5と決定した。全350方向Jacobiの一般rank98、二法線augmented rank99、n111の閉横断Artin長さ5を独立監査した上で、実smooth branchへのsectionと有限加群のrank一致から平坦性を証明した。n111での全projective完成環は C[[s1,...,s25,u,v]]/(uv,u³,v³)。CM type2、nilradicalの三乗は零、被約成分はこの点で滑らかである。群方向だけの形式積を仮定していない。全境界や大域正規性は未決である。

一般退化因子3H+滑らかな接quadric、modular Jordan(3,1,1)、一般PGL軌道23・連結安定化群Ga、正しいshearを伴う二変数軌道被覆を全150係数で監査した。全非零正軸では全Poisson環は非CMである。新成分の被約環の非CM性は意味しない。正軸は全HquarticとT21も避け、除外済み旧rank4成分は17件となる。他MPS帰属と全正軸枝は未決。

第52成分の全61出力dense S61写像を独立再構成し、凍結旧51カタログを変えず補足indexを追加した。新一次像は15次元、像上fiberは一次元で、一般betaを優越しない。その既約超曲面方程式の次数は全Sp4不変式空間の独立検証により16以上と証明したが、実方程式は未取得である。

14ページの独立日本語原稿 papers/positive_weight6_component_ja.pdf を追加し、LuaLaTeX・M2・全82入力hash・数学的scopeを監査した。旧118/16/91ページ原稿と歴史的SOURCE_INDEXは不変。既存57との区別は証明済みだが文献上の優先権は未確認。

J0は一次cofactorで37安全pivot・20×54まで独立監査し、author側の追加二単元で39pivot・18×52まで保存している。全coker消滅や欠損有限性はまだ証明していない。既存の全cubic欠損次元≤1はJ0上の有限性を含意しない。非平方adFBCJ開では新shadow証明で全点product rank≥35を得た。一般affine残余E42のunit、B18実交叉のodd ideals/primary構造、全成分分類を続ける。総下界58、rank4下界52、総上界未定。

## v63 M: 新25次元全成分、総下界58

直前保存 `50bce81f14f00db48ebe4b49b268d9d50278725a` は通常push済み。以下より前の57/51は当時の下界・既存カタログを記録する。

新しい全25次元既約Poisson成分 R4_POSITIVE_WEIGHT6_25 を証明した。補助incidenceのN²条件を外すと三変数の全次数Poisson familyが残り、全140 affine・350 canonical Jacobiを独立QQ演算で確認した。PGL・scalar込みの像微分rank26と、全350から得た真四次初期式のactual局所上界25が一致し、飽和像閉包は実際の全成分となる。既存57成分との区別も独立監査した。総下界は58、rank4下界52へ上がる。これは全分類や58の等号、文献上の優先権を意味しない。

新成分のprimitive一次像は全七変数の函数体rankを別engineで監査して正確に15次元と決定した。S61成分は16次元なので像上の一般fiberは一次元であり、全W16を優越しない。従って今回の新成分は一般beta五点尽くしの反例ではなく、その未解決課題は継続する。一次像のhypersurface式とfiberの全型は未計算。

新成分の一般退化因子は3H+Qで、Qは滑らかなrank5 quadric、Hは接超平面である。正の重み六軸を含み、その極限Qa8/9も含むため、Qa8/9の全局所・完成極小素数も41以上となる。全Qaへの包含は未証明。generic reducedness・全局所枝・正規性は未決で、固定QQ二点でtangent27が出たことだけからgeneric nilpotentを推論しない。

正軸の全8normal形式模型を構成し、真の全初期式を五次まで厳密計算した。HF=(1,8,31,72,126,171)、五次下位idealの根基は唯一の三平面L3と二平面P2。新実familyが三次元を与えるので、正軸のactual横断次元3・projective局所次元25は等号で確定した。全初期idealや他の形式枝は未決である。

一般アフィン残余はさらに、一行をbetaだけで正規化した線形42行列の冪等det E42へ簡約した。実unitはなお未取得。generic B18の実CA/CB交叉は共通R28上の二つのfree rank2 normalizerの個別image kernelへ還元した。CHにより四行有限module M/Uとtarget作用で実全環を復元できるが、係数idealと全primary型はまだ未計算。J0は五整数単元pivotにより36pivot・21×55へ進んだ。研究と通常pushを継続する。

## v63 L: 大域idempotentと次数42残余・Darboux同定

直前保存 `01bc0d884b54e3027a4bd3cb83c13c32cecc5395` は通常push済み。

一般アフィン代数を全schemeとして A=E_known5×A_res に分離した。全134二次式のpolarizationとadjugateから次数45のtrace idempotentを得た後、全129残余を線形normal行列A129×42(q)nで表し、次数42の一式R42によりI_res=I+(R42)=I+(e)と厳密決定した。betaだけの逆元を使い、cyclic/Hや未知qで切らない。未知残余の空性はなお未証明である。別にbetaのみで全35変数124式へ縮約し、旧quintic開では元critical商の正確な5×680表示を得た。

全Qaの完成極小素数下界は40、a0,2/9,−8/9では41以上。全非零aの既知MPS40所属はT01–T38に加え、T39 iff2/9、T40 iff−8/9と完全判定した。a0の両楕円成分の非所属や未知の全枝は主張しない。

正の重み六軸は一つの22次元projective軌道で、全T01/T02/T03/T39/T40、CA/CB、二quadric成分、六つの普遍特異jet成分の外にある。全部で十五の既知rank4成分から排除したが、残るMPSとHquarticへの所属は未決で、新成分と数えない。全軸にはJacobian1の二次Darbouxと(2,3)Poisson Cremonaを構成した。さらに全Qaに同じa独立の三次三角Darbouxを得て、全aのPoisson正規形をa∂X∧∂W−2∂Y∧∂Z/9とした。a0の全affine葉は二CasimirのA2 fiber、非零aではA4全体が一葉となる。

非平方では二定数の短い単元pivotでfull26×60、J0では別の四定数単元pivotで22×56へ厳密縮約した。全各点rank・欠損有限性・quartic/ambient除去は未解決。保存後は一般残余のunit、正軸の全小transverse模型、実B18境界交叉と非平方縮約を続ける。下界57、上界未確定を維持する。

## v63 K: 全Qaの四十成分と交叉の全非零a型

直前研究保存 `64c5c455ac642cbeacbf59fc63fd638e80102712` は通常push済み。障害時の記録のみ `1986a0ec09113a291dee46cfd98a1a1f6c0ddead` をGitHub APIで保存し、復旧後に取り込んだ。

全有限Qaが実MPS成分T01–T38の全三十八に属することを証明した。可換冪の実canonical toric arc、十七QQ etale source、十六の単元Jacobian付き正のHensel sourceを独立監査した。全域CA/CB所属と合わせ少なくとも四十既知成分が全Qaを通る。新成分の追加や全枝の尽くしではなく、全体の下界57は変わらない。

全非零Qaでは個別CA/CB/Nの全pinched完成環をaだけの反転で決定し、共鳴a2/9も含めた。実交叉Dの全218関係・colonを保つ有限target平行移動も得て、Q1の主平方零nilradical、smooth20次元embedded prime、depth20を全非零Qaへ延長した。a0の全型と全Poisson枝は残る。

R19のgeneric B16軌道の全完成像はC[[u16,t,r,c,p²,pr,pt]]、quadric normalizer r²−25t²−2cp=0を持ち、dim19,edim22,depth18,mult3である。B16の全特殊点や実R∩Nの厚みへは延長しない。

全重み六閉切断は二special軸a=±8/9を持つ厳密primary idealで、totalはCM一次元、special fiberはdepth0。負軸は実T40。非平方では短いQQ単元pivotと有界次数lift不能を独立確認し、v8/J0の続行証拠も保存する。全各点rank・四次被覆、一般アフィン五点の全支持の尽くし、全大域primary分解は未解決である。

障害後に同じworkspaceへ復旧し、開始時main64c5c455から緊急引継ぎ1986a0ecを確認してfast-forwardした。M2 1.24.11/default QQ GB、FLINT3.1.0、SymPy1.14.0を軽量再確認し、セットアップは再実行していない。二度目の短い接続断も復旧した。新成果の通常commit/push後も研究を続ける。

## v63 全Rのnormal/CM locusとRat113の実prime

直前保存 `b5545eeadc928bb5a0cb11cb2d9906f7113eedc2` は通常push済み。

全Qaが実MPS成分T02とT03にも入ることを、二自由固有値の実toric arc・明示曲線・a像の密度により証明した。全Qaは少なくともT01,T02,T03,CA,CBの五つの既知成分を通る。全局所枝の尽くしや新成分ではない。

R19の全normal proper quadric模型を構成し、全被約R∩Nを既約B18（18次元）と決定した。B16⊂B17⊂B18は16/17/18次元で、R∩Nf=B17。Rは全体でR1を満たす。さらにR自身の全特異・非正規・非CM locusは正確にB17∪C15であり、C15は一つの15次元軌道である。

B17一般点のR完成環は二つの滑らかな19枝が17次元で交わる C[[u17,a,b,c,d]]/(ac,ad,bc,bd)、depth18。C15では全完成像環をC[[u15,B,E,F,p²,p³,Bp,Ep,Fp]]と同定した。dim19,edim23,depth16,mult2、normalizer一点fiber長2である。実交叉D=CA∩CB自身とは区別する。Dのminimal一般長はN1,R2、R上のDとDredのCM locusは正確にR\B18である。

全37模型内には実Rat113の正則A20商を全係数で構成した。全350式のcofactor C113によりa*C113=0を明示し、C113はQa0完成でも非零。全R37のa0非平坦性と厳密なprime peel I=N113∩(I+C113)を得た。残余は未分類である。

一般アフィンの全正次元支持を含む元45座標のcritical行列式Eと、凍結点を通る両立する五変数断面を得た。実Unit(I,E)は未取得。非平方では全27×61の厳密縮約と有用な失敗receiptをv6 archiveへ保存したが、各点rankや全四次被覆は未解決である。

共鳴点の通常初期式は全四次まで、377方向のcancellationを含めて厳密に検証した。その全てが二つの15平面で消えるので、四次以下からの12切断Artinian目標も不可能である。実blowupの別経路は解決した。厳密な19一次pivotと二法線の全二次を覆う3minorから、completed unsat商をQQ[[δp14]]上有限にし、r飽和strict transformのexceptional fiber次元≤13を得た。従って二つの15平面はいずれも真の全通常接錐には含まれず、交叉はproper closed coneで次元≤14となる。指定点や平面の全方向の除去、全模型の上界14/12とは区別する。全実Poisson枝、全一般アフィン五点の尽くし、成分上界は引き続き未解決で、下界57を維持する。

## v63 全a有限模型・実五点開・通常接錐の残余

直前保存 `6d82c007d497c5e6edbfdd006274d3133bc405fd` は通常push済み。

全a上で定義される37変数の多項式有限模型、aだけを反転する35変数模型、a(9a−2)だけを反転する34変数模型が完成した。元の全350Jacobiと150係数を保持し、frame微分を含む全群作用Jacobianも定数単元として検証した。古いframeの余分な例外は不要になった。

全Qaは実被約T01の特異点であり、その接次元は全aでaffine30以上・projective29以上と証明した。二つの可換行列source表示の実際の微分像の和を用いる。二つの局所枝や正規化fiberを数えた結果ではない。

一般アフィン問題では、全134式を保つ厳密monic五次・5×675表示を指定非空開で得た。さらにnorm開で全16方向を零にして、既知五点の実有限自由5代数を分離した。係数は有限逆行列レシピで保存し、その補集合は未解決である。

共鳴点の真の全通常接錐のη=z33非零開は、T39の線形接平面の被約滑らかなGm×A8に厳密に一致する。他の接錐成分はη0に限られる。一方、下位ideal J=(57Q,17G)自体には二つの15次元錐があることを全係数で証明した。従って旧12一次切断をJだけからArtinianにする目標は不可能であり、新しい真の四次以上の初期式を計算する。これは全Poisson局所次元15以上という意味ではない。

真のa0二方向だけを残す閉切断も、全350式の両向き証明書により(as,at)と同定した。a0平面は実際のRat113族に入り、全schemeは被約だが、この閉切断の非平坦性を全R37へは移さない。

成分数の下界57、全一般アフィンの尽くしと上界未確定を維持し、保存後も研究を続ける。

## v63 全Q_aのT01帰属と全パラメータの一次構造

直前保存 `0c8877d94f23ad136acbe41adcd5213bc0649e39` は通常push済み。

全Q_a族が実際のT01へ入ることを、全150係数を持つregular nilpotent表示とcanonical多項式退化で証明した。a=0も含む。CA/CBの全proper模型への直接持上げも合わせ、全aで実既約成分T01(25),CA(28),CB(28)が通るため、全Poisson完成環とその被約空間の完成環は非CMである。全枝の尽くしではない。

全126方向のJacobiと全GL5軌道をQQ[a]上のSmith変換で決定した。Jacobi rankはa=0で73、a=2/9で74、他は75。affine軌道rank17は全aで一定である。古い選択frameの余分な因子は内在的な軌道の例外ではない。高次の全a模型への改善は別の実作業として続ける。

全Rは既存LPB39だけでなくLog(1,1,1,2)28にも入る。元のlog builder全150係数の定数倍一致を確認した。全R\Nには二次Poisson Cremona、有理Casimirs、全函数体の逆式を構成し、一般葉の有理性を証明した。F非零の一般葉閉包は次数6の有理cone。R∩Nへの葉の分類は拡張しない。成分数の下界57、上界未確定は変わらない。

## 確定結果

v63では、共鳴点Q_(2/9)の全Poissonスキームの横断完成環を、35変数・350残余式の正確な有限同変模型にした。74変数の陰関数は正重みによる有限多項式であり、高次jetの打切りではない。17群切断と全係数の独立再構築、QQ密入力、全共鳴曲線を検証した。全枝の分解は未完了だが、既知12次元枝と9次元枝の差から全横断環は非CMと分かる。

重み6だけを残す五変数の閉切断は、u=a−2/9として I=(uη,(w0,w1,w2)(w0,w1,w2,η)) に厳密に等しい。全350式と十生成元の両向き恒等式を検証し、二つの被約枝、一般厚み4と1、CM・次元1・重複度5を決定した。この閉切断の枝を全35変数の枝へ読み替えない。通常二次障害57本は厳密な13平面を含むため、二次障害だけでは全模型の次元12以下を証明できない。

CA/CBの共通9次元一次像のrank3開では、全被約ファイバーは二つの埋込A10で、その交叉はA9である。閉properな旗入射を全点で構成した。β=0では全被約CA0=M0(13)、CB0=M1(13)、N0=R0(12)と同定した。指定Q_(1)で個々の被約CA/CB/Nの正規化ファイバーは長さ2の一点、各完成環はdomain・重複度2であり、全体でR1だが指定点で非正規・非CMである。Q_(1)とQ_(2/9)を区別する。

全134式の一般P31切断も被約五点と証明した。全17底上では、cyclic行列式を反転しない実annihilator det(T_t−I) と分母を払ったresultantを構成した。これらの零集合、固定二座標の一般rank、全アフィン五点の尽くしは未解決である。非平方のJ=5ad−6e−4=0では二次積rank41の欠損を有限へ落としたが、J≠0と残る四次欠損は継続する。

正規な実退化因子の認識・分類の証明を、再帰的な旧TeX入力のない日本語91頁の独立原稿へ統合した。Pymのpurely elliptic仮定との違い、幾何ファイバーの正規性という相対仮定を維持する。全分類・成分上界は未完成で、下界57を維持する。

## 再現の入口

- `reports/Q_WEIGHTED_FINITE_MODEL_20261103_v63/report.json`：全有限模型・正確な係数検証。
- `reports/Q_WEIGHT6_CLOSED_SUBSLICE_20261103_v63_REVIEW.json`：両向きideal恒等式。
- `reports/Q_RESONANCE_NON_CM_20261103_v63_REVIEW.json`：全模型の非CM性。
- `reports/Q_ORDINARY_QUADRATIC_CONE_20261103_v63_REVIEW.json`：二次coneの13平面という厳密な限界。
- `reports/AMBIENT_FIBER_BUNDLE_20261103_v63_REVIEW.json` と各零fiber・正規化・R1レビュー。
- `reports/S61_AFFINE_BASE_PLANE_SUPPORT_20261103_v62_REVIEW.json`：全P2底と一般P31。
- `reports/S61_GLOBAL_ANNIHILATORS_20261103_v62_REVIEW.json`：反転前の全実代数への所属。
- `reports/NORMAL_ACTUAL_DIVISOR_PAPER_20261103.json`：91頁・三回LuaLaTeX・210labels・参照/組版検査。

## v63 既知局所枝の実分解と全域R1

R\Nの全軌道はκ²一変数族と六単独軌道へ完全分類し、全点の交叉をordinary doubleと決定した。R∩Nは別の残余である。非平方では全U=aFBC非零上の三次row moduleを自由57、全関係核を自由6として明示した。product submodule自由41・抽象商自由16も証明したが、各点の像階数下界や四次被覆はまだ従わない。

直前保存 `153c22630f05f2b04dfd187ae8b3cb3aef82e204` は通常push済み。

全射影CA/CB/NのR1を証明した。proper模型の全幾何ファイバーを内在的に一意へ戻し、源の全微分・長さ1・finite/Nakayamaからrank2境界因子のtarget滑らかさを得た。Q1の非正規・非CM性は維持する。R19は有理で、一般正規形の完全軌道不変量はκ²、PGL軌道18次元。実交叉はその正規形開の全点でordinary double C[[u19,ε]]/(ε²)である。一方R全体は既存LPB39にも入り、Rの全閉点の全Poisson完成環は28/39次元minimal成分の差により非CMである。

全35模型の既知T39/T40枝について、実際の商をそれぞれQQの9変数多項式環と同定し、実350式の差商cofactorから明示的なscheme開のideal等号を得た。さらに正確な全ideal分解 I=N39∩N40∩K、K=I+(C39,C40) を証明した。全Rはa=2/9と−8/9で非平坦。Qstarの残る対象は明示的なI+(C39)であり、その全枝・厚みの分類はまだ残る。

全109変数slice上で群17minorがaのみの単元であることも全係数・重みにより証明し、既知A9開のcanonical affine26/projective25のregularityを正当化した。全一般アフィン尽くし、成分総数の上界は未解決で、下界57は変わらない。

## v63 全射影交叉の確定

全35模型には、固定平面上の全五偏導関数の140×66行列の明示的な66minor Δを構成した。Δ非零では実退化因子が幾何正規で、被約支持は既知T39に入る。Rstar曲線上ではΔ=非零整数×ε59。残る他の枝はΔ0内に限られるが、その尽くしは未決。

全35横断模型の局所次元については、実三次初期下位ideal内の13個のcofactor付き多項式とproper Projにより、無条件のQQ/C上界22を得た。これは有限体の標本からの推測ではない。既知最大枝の次元12に対する上界12、全枝の尽くしは引き続き未解決。

直前保存 `0713e279a55bcf81590242e3b908509a3cefa499` はGitHub mainへ通常push済み。

全射影CA/CBの被約交叉を完全に決定した。結果はN27だけではなく、階数二の19次元既約族R19を加えた正確に二つの成分 N27∪R19 である。33相対旗×二端点の全66交叉ベクトル空間をQQで尽くし、具体点が全proper N像の外にあることをcanonicalな内在条件で証明した。RはS61に入らないので、従来のS61一意性やQを通る交叉成分の一意性と両立する。Rは全Poisson空間の新成分ではなく、下界57は変わらない。

全CBにも非自明associated P1束を使うproper28次元模型を構成し、全rank≤2境界は26/27次元の二成分と証明した。Nの全rank≤2境界もCBの二端点の平方型として26/25次元の二成分へ認識した。これらの全域R1・正規化・scheme厚みは別である。

共鳴全35模型では、全350式の通常三次をQQで再構築し、真三次初期類17個とHF0..3=(1,35,573,5859)を決定した。二次で残った偽13平面は三次で二つの12平面の合併に切れるが、全局所次元12以下はまだ従わない。

一般アフィン問題には、全134式を保持する六超越パラメータのmonic unit判定を得た。unitなら全一般fiberは零次元、逆も成立する。これは正しい有限判定への帰着であり、実unitや全長5はまだ得ていない。7パラメータ87共役Wronski判定も別に検証した。実solverの時間・メモリ制約は別receiptへ保存し、失敗を空性に読み替えない。

## 後続保存の確定結果

指定Q_(1)の実交叉についても全49変数218関係の完成局所表示を決定した。nilradicalは主平方零(F)、全加群はC[[u11]]F、associated primesはNと滑らかな11次元embedded層の正確に二つで、S61の次元18・depth11。独立なNの94生成kernelと両向きideal/colon等号を検証した。環の分裂は主張しない。

直前保存 `c18e7420b145ae0b486dba679d018e3dcd0bcc86` は通常push済み。

S61の全被約交叉をNと同定し、Qを通る交叉成分の一意性を証明した。ただし全射影閉包のgeneric Poisson rank≤2だけの追加交叉成分は、S61では見えず未検査である。実際の交叉はQ_(1)で非被約であり、埋込associated primeを持つ。正規化の微分だけでなく、個々の被約targetの全cotangentを比較した結果である。

個々の被約CA/CB/Nの全完成環も決まった。S61でAc=C[[u10,s1..sc,p²,p³,s1p..scp]]と書けばCA=CB=A8、N=A7。全射影模型では次元28/28/27、embedding dimension37/37/35、depth20、多重度2である。共通一次像Zの局所環も七滑らか因子を持つpinched surfaceとして同定した。これらのR1はS61およびrank4開の主張で、全射影rank≤2境界への自動延長はしない。

共鳴35変数模型に加え、指定Q_(1)を含むD(a)(9a−2)開で全34変数模型を構成した。実残余296の非線形269項を全て残して単元係数で一変数を消すため、全次数のideal等号である。共鳴点では35変数を維持する。

全射影CAには9次元旗多様体上の滑らかなP19束というproper全射模型を構成した。これを使い、全generic rank≤2境界の被約成分は相異なる27次元・25次元の二つと決定した。全正規化、境界のR1、CBとの交叉は別課題である。

一般曲面の判定は、正確な9パラメータ次数問題と8パラメータunit問題へ減らしたが、その実rank/unit結果は得ていない。既存166/191模型と有限体空性をNakayamaで結ぶには、全20法線を忘れたp整な実アフィン商の有限生成性が不足することを明示し、逃走点の反例で誤った推論を排除した。

## 未完了と次手

一般アフィン未知曲面・曲線の除去、一般五点の尽くし、特殊一次像上の未知成分、
既知全成分の定義prime、全大域根基一致、成分交叉と非被約構造、全軌道分類を維持する。
細分化した次手は `research/NEXT_TASKS.md`。下界57、上界null。

上記の基底Z・個々の成分のpinched完成環・全S61被約交叉は独立監査を完了した。
Q1の実交叉の全localideal・nilmodule・全associated primesまで解決した。元の全Poisson模型の高次初期関係と、全射影rank≤2交叉を継続する。
通常三次以上の初期関係と、残る非平方QQ rowmodule、8パラメータWronski unit検査を継続する。

## 実環境と計算上の限界

`source /workspace/cloud-setup/activate.sh` を使用。M2 1.24.11、FLINT3.1、既存venv、LuaLaTeXを再利用。
M2 QQ LinearAlgebraは配列順序エラーを確認したためdefault GBを使う。
一台16GiB cgroupで重いsolverを調整する。一般係数のGBは6GiB cap/実時間上限で未完のものがある。
失敗は証明として使わず、同条件を反復せず構造削減・小さい証明書へ戻す。
非斉次modular標準基底の候補だけでは逆包含を未確認なので受理しない。
有限体のRational変換には分子分母を明示する。
巨大な展開・実行バイナリはignored workに残し、小さな再現ソースとhash付き証拠だけ保存する。
ユーザーから存在した `reports/20261102T142937.117200Z-all/` は未追跡のまま保持する。
ChatGPT Libraryへの保存成功receiptはなく、GitHub通常pushを実保存とする。
