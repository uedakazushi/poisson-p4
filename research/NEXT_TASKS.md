# 次の課題 — v63-s61 / 2026-11-03

## v63 AD / 2026-11-04: 十一の四次全域関係・全係数写像・joint Spin4体

直前保存 main `a0812a6f624600b54c214817126a5047c1622bad` は通常push済み。以下は過去時点を含む。

全Schouten/Bianchi作用素を実際に構成・検証した結果、直接の重み四部分は旧予測30行ではなく15行で、最高空間次数の84行は全て零だった。その33 genuine quartic定数blockから11個の既知単元pivotを得て、全15関係を354→319へ正確に輸送した。新11方向はaFBC全域、旧22との合計33方向はFULLUで独立。80方向まで47不足し、全四次全射は未解決である。さらに全aCで六quadの定数階数6を証明し、aFBCでの旧Z2欠損を除いた。

五source回路上でF4、五pencil係数、N12の全係数を小さい四shape変数の有理多項式として正確に得た。新点は使わない。独立な再構成で全七写像を確認し、六つの低不変式の普遍的代数独立性も同じ保存点のminorで証明した。この六生成元による24次41 monomialには、点なし係数行列のrank41から関係式がない。全24次空間は67次元なので、次数下界は24のままで、NEW方程式の不存在とはしない。残る方向を扱うには全degree24空間の具体的basisが必要である。

連結Spin4のjoint biform V(2,1)⊕V(1,2)については、完全normalizerと固定体の逆式から有理六変数不変函数体を決定した。W16へ必要な追加normal-lineのU/V交換は別問題であり、fullSp4/W16の有理性や全不変式環を主張しない。一次文献の限定検索も保存し、近接するsingle-biformの定理を直和へ誤用しない。

一般betaのδ_E∈J0、実full cubic57・全四次239、generic B18の全非被約Dとdepth19〜20、NEW定義方程式・target函数体index、全成分尽くしは継続する。総下界58/rank4下界52/総上界未定、五つの凍結原稿は変更しない。


normal formのcross-block edge積ghijにも正確な不変式G32を与えた。NEW上での全32次展開は不要で、七つの保存済み写像のB分母と次数から一意な最高B5係数の非零を証明した。G32候補を除外し、一般NEWがLaurent ghuv chartに入ることを確定した。T線形回収のL単元性や、NEWのmultisection次数の証明ではない。

## v63 AC / 2026-11-04: 固定四次欠損・有理source・単位frame接続

直前保存 main `f39105f754edf85d568f884930d7a4c0e88fe068` は通常push済み。以下は過去時点を含む。

固定319×239四次行列の欠損がaFBC開で余次元2以上となることを証明した。既存の相対自由65模型を、nilpotent kernel・CM性・正則cutの変更で今回の固定cutへ移す。正確な低重みHFの和65から全height1でweight4 quotientが零となる。全域独立80本のsource kernelがあれば全射を導けるが未構成。直接の重み四以下Schouten/Bianchiだけでは最大30新行の見込みで、旧22本と合わせても80には足りない。実際のoperator抽出と単位pivotをbounded jobで進め、全点239を先取りしない。

四次lazy APIと既存quad/cubic57のframe接続を正確に完了した。全33 genuine quarticの単位row rescalingも含み、旧H²係数は保持される。今後の実122×64係数計算はこのadapterを通す。未補正のT34/Tnを混ぜない。

新成分first jetに六sourceの明示的有理共変量回路を構成し、GL16・perspective4を加えた実30×26微分のrank25とsource kernelを独立に決定した。微分の情報、この一方向の有理積分、および実target商への有限indexを区別する。

元の被約Dredのconductor商について、被約環はC[[g19,w]]、nilradicalはG19正則係数環上free5、横断H0長さ5と決定した。これは元の非被約Dのnilではない。後者の全型とdepth19〜20、一般betaのδ_E∈J0、NEW方程式、全成分尽くしは未解決。総下界58/rank4下界52/総上界未定、凍結五原稿は変更しない。


このsource kernelの有理積分も完了した。五つのfirst integralとconicの全有理逆式からsource定数体は有理五変数体と決まり、t=2のsectionで実β・共変量を1032-nodeの五source回路へ縮約した。既知有限involutionの固定体も有理五変数体だが、実NEW target不変体との有限indexは未解決である。次は新点の高次数probeではなく、この小さな回路の低次数不変式を正確に展開できるかpreflightする。

普遍残余曲線側では全一般交点代数を八次monogenic表示と線形T回収へ還元した。これは全8次代数であって、選択した実成分が八次・既約・被約であることや、NEW multisectionの次数を意味しない。種数も引き続き未決とする。

## v63 AB: 元の被約conductor・embedded21境界・有限ΣのQQ四次階数

直前保存 main `972753a4c8131c35c50bb3f30e003c5c8edf3a09` は通常push済み。以下は過去時点を含む。

generic B18で元の被約交叉Dredの全normalization conductorを決定した。N側は実九odd係数を使う17生成上界、R側はT4である。Nの元のconductorとの差は、G19を除くと長さ4のcyclic module、annihilatorはHF(1,2,1)の完全交叉となる。全二変数source四次・18写像・even-base gaugeを独立検証した。非被約Dの全nilと正確なdepth19〜20は未決である。

既存embedded21支持がB18を通ることも、全多項式E(W,λ)と群作用で証明した。21小行列式から一つの正則形式枝を得たが、全germの滑らかさは主張しない。B18のnil支持はG19横断で少なくとも二次元であり、三つの曲線だけという模型は使わない。

非平方の六quadに全22 syzygy加群を構成し、像加群の局所自由146とpure四次商のreflexive64を証明した。周囲での特殊化rankは別で、孤立した非自由点はなお可能。既存proper有限Σ全体については一回のGF block rank2184とsplit有限NakayamaによりQQ quad78・元146/152を証明した。他のG2/vertex有限例外と実122×64全射、full cubic9×43は引き続き未解決である。

W16の固定特性基底における残余曲線を二つの明示的(2,2)式で表し、有限ordered拡大と相対符号quotient後のm次数上界8を得た。有理trace不変量mは五Inv8から独立である。NEW multisectionの次数や方程式の上界ではなく、種数も未決。NEW方程式次数下界24、一般betaのδ_E∈J0、全成分尽くしは継続する。総下界58/rank4下界52/総上界未定。


元のDred/conductorもconeでdim20/depth19/e1、projectiveで19/18/e1と決定した。有限長4の差を含む完全列のlocal cohomologyによる証明であり、元の非被約Dのdepthへは移さない。

実四次の次の係数計算には、今回構築・独立検証したlazy APIのframeを明示して接続する。APIはT34、既存quad/cubic57はT34*diag(Tden)で、Tden=[a,a,a,1,a,B/50,a,a]。これはaFBC上の単位変数変換で、全33 genuine quarticの単位row rescalingも必要である。未補正の二つのframeを混ぜて122×64を評価しない。API接続のexact receipt・独立監査と、Schouten/Bianchiによる全四次syzygyの新しい構造的検討を次に続ける。

## v63 AA: 全W非零の実交叉完成環・特性基底への優越・実122行判定

直前保存 main `75a9f688f5932c996c25e250943bc830656698bb` は通常push済み。以下は過去時点を含む。

全W非零における実交叉Dの全完成環を解決した。W1を明示的GL frameで既存critical-beta全sourceへ移し、全150係数・61座標・134関係、旧49変数のv4=v7=1点、全unit blockを独立に照合した。旧218関係と全colon・reductionを正当に局所化できる。nilradicalは主平方零で正則21次元上free1、associated primesはminimal Nとembedded21の二つ。coneのdim/edim/depth/eは28/37/21/2、projectiveは27/36/20/2。全W非零へ群移送し、W0=generic B18の未解決な全Dへは延長しない。

NEW25 firstjet像は全8次spectral base A5とF4を含む全特性基底を優越すると証明した。同じ既存QQ点だけで実source微分のrank5を得て、独立reverse ADと全共変量微分で検証した。従ってpencil係数だけで作る非自明な判別式など、基底のみの真のdivisorという候補を全て排除できる。NEWの未知条件は、残る一次元の普遍商fiberにおけるgenerically finite multisectionである。その次数・種数や実方程式は未決で、次数下界24を維持する。

全四次は追加DeltaH7・DeltaQ146 chart上で実122×64 circuitへ厳密に縮約した。全33 genuine quarticをunit2×2 H²操作で保持し、元四次cokernelと同型、全rankは175+rankC122である。低重み90×64は十分条件であり、full cubic57の点別rankを仮定しない。実rank64は未証明。次はこの実残余、generic B18の全nilとdepth19〜20、一般betaのδ_E∈J0判定を続ける。総下界58/rank4下界52/総上界未定。五つの凍結原稿は変更しない。

## v63 Z: 長さ13の実Noether fiber・二adjoint有限切片・次数下界24

直前保存 main `f9e849842bc1c19243acadebcfe8d5cd133af466` は通常push済み。以下は過去時点を含む。

generic B18の全非被約交叉Dに、物理28座標による有限な単射R28→Dを構成した。generic加群rankは2。一変数source四次と全87独立graphを独立Lie-log法で検証し、閉Noether fiberを正確にC[z,Y10]/(z³,zY_i,Y_iY_j)、長さ13、HF(1,11,1)、socle11と決定した。R28加群の最小生成元数は13である。閉fiberの冪零次数を元Dへ移さず、元Dのdepthはcone19〜20のまま、全nil構造も未決とする。

W16の二adjoint共変量対について四つの普遍関係を与え、一般pairからsourceを符号±だけを除いて復元することを証明した。実像は四式零点の一つの一般滑らかなcodim4成分で、有理Sp4不変函数体も一致する。regular semisimple Dと二つのtorus gaugeの明示的開では、四式idealは三つのgraphと七変数77項の一式に正確に書き換わる。これは6次元の有限軌道切片であり、NEW25にはさらに未知のdivisorが必要である。

NEW25 firstjet超曲面の方程式次数下界は24へ改善した。全Inv22の次元21に対し、既存36個のQQ点だけで制限rank21を正確なQQ小行列式と独立監査で確認した。実方程式や24次での存在は未証明で、24次probeは行わない。

非平方ではternary pencilの全leading caseをQQ恒等式で尽くし、共通因子例外が有限であることを証明した。B²を除いたpairの全域rank2は正しいが、初稿のB*S/B*U比較による五式syzygy推論は撤回し、別のERRATUMを保存した。追加gate G3またはintegral scrollを明示した修正版quad78補題は受理した。全三次57、全四次全射、一般betaのδ_E∈J0と全成分尽くしは継続する。総下界58、rank4下界52、総上界未定。五つの凍結原稿は変更しない。


全W非零の交叉接空間も全150物理係数で独立に決定した。CA/CB/Nは38/38/36、実Dは37で、被約D=Nとの差は正確に一つの線形nil classである。ここから新embedded成分や全完成型を推測せず、既存の全218式模型との実座標同定を次に進める。

scrollのE5非零開での幾何学的integralityと、E5零・vertex rank不足の有限交わりも証明した。これにより修正版の全gateを合わせ、全aPJ開でraw quad78欠損が有限と決まる。元のdegree3 colonを六Koszulで厳密に持ち上げると、FULLUで有限例外を除きlight146・weighted152となる。有限例外の実rank、残る122×64全quartic判定、9×43 full cubic判定は未解決である。

二種類のR28を区別する。今回のfinite13はconeの物理R29からfold行61を除きscale行106を残したR28である。旧rank4 union-module模型のR28は、fold行61を含むcone R29からscale行106を除いた射影底である。後者の二normalizerがfree rank2という性質を前者へ移さない。

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

## M以降: 下界58を基準とする

新25成分R4_POSITIVE_WEIGHT6_25の全成分性・旧57との相違は解決し、総下界58/rank4下界52へ更新した。新成分のgeneric scheme厚み・二normal Artin模型、全S61写像、既証15次元一次像の定義式と一次元generic fiberの全型、実軌道moduliと原稿を続ける。正軸L3を除去して局所上界24を証明する目標は実三変数familyに反するため再試行しない。一般beta五点の全尽くしは別課題であり、新成分の一次像優越性を次元16だけから推論しない。全generic affine残余E42のunit、B18実交叉のodd ideals/primary構造、非平方のrank・四次被覆は未解決のまま続ける。

## L以降の実行可能な残余

一般アフィンの既知五点は大域idempotentで分離済みで、未知schemeは厳密にI+(R42)。このunitを証明するか未知支持を構成することが次の有限目標。全beta-only35模型と129×42線形normal行列を再使用する。完了degree3 Macaulayが同じ中心で線形帰結なしなので、同値なD1leftinverse探索は繰り返さない。正軸は十五の既知rank4成分から除外したが残るMPS・Hquarticが未決。全physical350から得る8normal formal模型と高次obstructionを扱い、有限jetのconeを全枝と同一視しない。全QaのDarbouxと40/41下界は解決済み。実B18境界交叉と非平方の残余単元・四次被覆を継続する。

## K以降の優先入口

全Qaの実T01–T38所属と計四十成分の下界、全非零QaのCA/CB/N完成環・実交叉embedded構造、generic B16の全R完成環は解決した。これらを再試行せず、一般アフィン全支持の残余、非平方の実各点rankと四次被覆、B16特殊境界と実R∩Nの厚みを続ける。有限arc sourceのT39/T40必要条件を全成分非所属へ拡張しない。二つの重み六軸のうち負軸はT40、正軸の帰属は未決である。保存後も研究を継続する。

## v63で解決した段階と現在の課題

指定平方層Σの全一般β上の除去、既知51族の全S61写像と一次像次元、L1の交差次元27、L2の共鳴追加方向の全次数積分は完了した。v63でL2の指定完了条件である全35変数の有限残余idealも得た。正規な実退化因子の長い証明鎖は91頁の自己完結原稿へ統合・監査した。これらを未解決へ戻さない。

- **一般アフィン尽くし。** 全17底・1031記号の実作用と、反転前の大域annihilatorを使い、全t可逆支持を尽くす。P31切断の被約五点と、全一般ファイバーの五点は別である。固定二座標の実一般rank≤29なら曲面次数下界30と矛盾するが、上界はまだない。六パラメータmonic detector c^45−Σq_i c^i のexact unit検査は全fiber零次元と同値。既知五点を含む指定norm開では全16方向を零にした実自由5代数を得た。補集合H0とcyclic開外の全t可逆支持を扱い、全体の未取得unitを次の有限目標とする。8パラメータWronski射と87共役版は曲面判定の別経路で、射影center回避や次数30を仮定しない。
- **非平方の残余。** J=5ad−6e−4=0の二次積欠損はQQで有限と証明した。J≠0の全二次積被覆と残る別の三次・四次欠損を扱う。全Uの三次row module自由57・全六関係が完成した。30pivotを保つad-basisから27×61へ縮め、各点rank欠損支持と四次被覆を扱う。短い単元Bézout証人による実pivotを使う。巨大小行列の落ちを全周囲の障害の落ちと同一視しない。同じ条件の重いGB再試行を理由なく繰り返さない。
- **特殊一次像と交叉。** rankKβ=3でのCA/CBの全被約A10ファイバー、共通A9、β=0でのM0/M1/R0は完成した。全S61の被約交叉はNと解決し、Q1で実交叉の非被約性と埋込素idealの存在を得た。Q1での全実intersection完成ideal、主平方零nilradicalと11次元prime annihilatorまで決定した。全域のscheme構造と残るprimary表現を進める。全射影CAのrank≤2境界は27/25次元二成分と決まった。CBのproper模型と全rank≤2境界、全射影被約交叉N27∪R19は解決した。CA/CB/N全体のR1、R19の有理性・一般軌道不変量κ²とordinary double開まで解決した。R\Nの全軌道・ordinary doubleは全点で解決した。全被約R∩N=B18、R自身の全Sing=Nonnormal=NonCM locus B17∪C15、generic B17の二枝環とC15の全pinched完成環は解決した。実DのR上のCM locusとminimal一般長N1/R2も解決した。残る実R∩Nの厚み・embedded層、全大域primary分解と全finite正規化を進める。指定Q_(1)の長さ2正規化・R1・非CMと、共鳴Q_(2/9)を区別する。
- **全35/34変数模型の残余構造。** 共鳴Q_(2/9)の35変数と非共鳴Q_(1)の34変数の正確な有限DAGは完成した。被約枝・primary構造・全初期idealは未完成である。二次障害は13平面を含み、これだけで次元12以下とは言えない。真三次17類と57二次の74初期生成元は完成した。全35模型の局所次元上界22を得た。明示的Δ開の被約支持はT39に入るため、既知T39/T40を実schemeとして剥いだK=I+(C39,C40)、QstarではI+(C39)を扱う。残余全枝・厚みが未解決。全Qa⊂実T01が証明されたため、残余KにはT01も含まれる。T01を忘れて残余がCA/CBだけと予想しない。可換行列sourceの全微分・局所枝を調べ、四次以上の有限証明、上界12への改善、既知primeとの比較を使う。全aの37変数、aだけを反転する35変数、a(9a−2)だけを反転する34変数の全高次模型は完成した。真の全通常接錐のη開はT39接平面と解決した。η0でJ=(57Q,17G)には15次元錐があるため、旧12一次切断のArtinian目標はJだけでは不可能である。真の全四次初期式も377方向のcancellationを含めて検証し、両15平面が全てで残ることが分かった。四次以下の同じ切断は再試行しない。actual blowupの有限性・飽和により両15平面の一般方向は除去した。それぞれ真の接錐との交叉はproper closed coneで次元≤14となる。残るproper閉部分と他の接錐支持を扱い、全次元上界14/12へは自動延長しない。Qa0では実Rat113 A20 primeと明示a*C113=0が完成し、残余I+C113を扱う。重み6閉切断の二枝と全模型の全枝を混同しない。
- **原稿。** 新91頁原稿の証明・仮定・有限排除の範囲は監査済み。投稿上の新規性の網羅的調査と、他の候補論文への独立再構成は残る。凍結v62主原稿のhashを無言で変更しない。

## Cloud 再開の入口

[v63引継ぎ](CLOUD_HANDOFF_20261103_v63.md)、[初回再開](CLOUD_HANDOFF_20261103.md)、[動作・接続記録](../reports/CLOUD_RESUME_20261103.md)を使い、成功済みのセットアップをやり直さない。Cloud M2 1.24.11ではQQのAlgorithm=>LinearAlgebraを避け、確認済みdefault GBを用いる。全Cloud研究の開始SHAと各保存SHAは引継ぎにある。

## 作業の選び方

全 S61 の被約既約成分の尽くしを優先する。楕円二成分の一般次数1・全域の二次 graph、
T01 の次数3・S3モノドロミー、既知五点の三平面断面は証明済み。これらを未解決へ戻さない。
各結果を保存した後も研究を継続する。証明済みの有限課題だけで完全分類を完了扱いにしない。

## v60 からの具体的な前進

- 固定二座標で零次元となったコロンの長さと t の可逆部分を計算する。191生成表示の作用を用いて、無限遠の有限点をアフィン部分から区別する。
- 指定因子への交叉数は曲面があれば20以上である。実局所商長の和19以下、あるいは全指定点での周囲冪零性を証明する。非平方 a=0 と d=0 は完了した。平方 c=2a の有限値は中央商が正次元になる二次数体の層を先に扱う。
- 一般五点の尽くし、特殊一次ジェット上の成分、全体の個数と一般元の幾何を扱い、軌道・交叉・非被約構造の分類へ進める。

## 具体的な開始点

1. **全一般ファイバーの三平面への包含。** `TORIC_UNIVERSAL_THREE_PLANE_20261101` から一般 β 上の42本のアフィン一次条件を得る。これらが全134式の根基に属することを、局所化・有効なイデアル証明書・幾何学的認識定理で証明する。包含が示されれば三平面との交叉代数 A_T×K×K により一般ファイバーを尽くせる。既知五点の滑らかさや平面断面の長さだけから包含を推論しない。
2. **明示された T01 成分を全体の消去に使う。** 一般二捩れ同定、C9 の全域多項式性、多項式 X により、Δ≠0 上の実際の T01 は正規有限エタール三次被覆と一致する。全 T01 の飽和イデアルを既知二楕円 graph と合わせて使う。既知三成分の合併を与えるイデアルと全134式の根基が一般 β 上で一致するかを証明する。
3. **全解を保持するアフィン消去。** 135式は三次 Macaulay 閉包を反復しても新しい二次・一次帰結を生まなかった。既存の120〜180秒の汎用GB試行を同条件で再実行しない。三平面に適合する42法線＋3平面変数でのF4試行も600秒で未完。新たな幾何学的削減・飽和・高次関係を先に作る。未飽和斉次化の無限遠の旗をアフィン枝と混同しない。曲線の c4,c6,j の輸送保存則は、Δ≠0でもT02の反例があるので追加方程式にできない。
4. **特殊層と境界。** 楕円二 graph の交叉は旗円錐と原点の長さ10冪零部分まで解決したが、周囲の全 S61 局所環は未決。rank1 固定ファイバー、CA/CB の共通9次元像上の全枝、T01 の非半単純境界を進める。20代数型の模型は全T01境界を覆わないと証明した。β=0で正確に10次元のファイバーがあり、実際の飽和イデアルの例外部分を解析する。
5. **既知成分の座標カタログ。** T39/T40 の素イデアルは q−Q51、q−Q52、T01 は三次被覆からの明示的な飽和イデアルと記述できた。T02 の一次像は15次元の既約超曲面である。他の既知成分の一般族・定義イデアルを順に構成する。57件の台帳は全51階数四イデアルの完成ではない。

無限遠の追加探索は v07 の Q45 の22成分を入口にする。固定定数因子 v を持つ R0 では、旗の distinguished vector u と v は一般に異なる。G=x²(az+bw) だけへ制限してはいけない。一般的な旗は G=ell³ modx³ の射影である。滑らかな incidence conic の部分と、その退化型を分けて全てを監査する。退化型では旗以外の無限遠部分が現れる計算があり、全無限遠が旗だけという予想を前提にしない。

各課題は普遍的証明または適用開集合を明示した証明を完了条件とする。原稿の統合・保存は研究停止の理由にしない。

## 優先0: 基準状態を保持する

**対象:** 原資料、今回の監査結果、数え方。

最初に `STATE.md`、`CLAIMS.json`、`COMPONENTS.json` と必要な旧版だけを読む。リポジトリのセットアップ手順で `legacy/` を展開し、保存済み検証記録を確認する。環境・依存版・入力ハッシュが変わった場合、または該当の証明コードを変更した場合にだけ必要な検証を再実行する。旧版 runner の既定出力が原本を上書きし得るものは、一時コピーまたは別出力先を使う。

**完了条件:** 57の下界、階数二の8/6、正規因子の依存鎖、S61とQ45の差について取り違えず、現行の入力と再現記録が特定できる。保存済みの検証を毎回すべてやり直すことは要求しない。

## 優先1: 既知成分を $S_{61}$ 上の一つのカタログに移す

**目的:** 既知の51個以上の階数四成分と、全体の方程式を比較できるようにする。

**読む資料:** `legacy/v02/p4_poisson_completeness_v2.tex` §3、同版 `data/symplectic_slice_reduced_equations.json` と `symplectic_slice_certificate.json`、v02 `baseline/` のMPS40カタログ、v3–v8 の追加11成分の本文。

**着手:** まず S61 の復元写像・正規化写像と134方程式の座標規約を記録する。MPS40の内部IDと一次文献の図式はv27で全照合し、`COMPONENTS.json` は個別展開済み。この対応を使って未構築の一般族・S61座標へ進む。次に追加11族を同じ正規化へ移す。容易な族から始め、全51個を一度に処理したことにしない。

**各成分で残すもの:** 原資料のID、一般パラメータ族または軌道代表、正規化の適用開条件、S61への写像、一般階数・次元・退化因子の型、全空間での成分性の証拠、他成分との区別。

**完了条件:** 既知成分の全リストが一意な内部IDを持ち、S61上でどの部分を与えるかが検証できる。有限個のサンプルが方程式を満たすことだけで完了としない。各族の一般パラメータについて恒等式または証明を持つ。

## 優先2: 16個の重み1方向を含む大域還元を作る

**目的:** 45変数断面の22成分分類から、61変数全体の分類へ進む。

**着手:** 重み1の16変数の残存群作用と、重み4・5方程式への現れ方を明示する。一般の非零重み1部分を扱える正規化、階数条件、共変量による分岐を探す。まず一つの稠密開部分を完全に扱う目標を立て、除いた閉部分を明記する。

**大域化の候補:** 実際の退化因子を「正規」「被約非正規」「非被約」に分け、既約因子の次数、重複度、特異集合、線形対称性を組み合わせる。次数分割だけでは尽くしたことにならない。正規部分では優先3の監査状態を明示して既存の認識定理を使う。

**正当な成果の形:** 次のいずれかを得る。

- 明示した開部分の全解が既知族に含まれることを証明し、残余を具体的な閉集合へ狭める。
- 新しい既約族が全空間の成分であり、既知成分と異なることを証明する。
- 全体を有限個の明示した問題へ帰着する構造定理を証明する。

**第一段階の全完了条件:** 既知族の素イデアルを $\mathfrak p_i$ としたとき
$\sqrt{I_{61}}=\bigcap_i\mathfrak p_i$ を証明するか、それと同等な全点の被覆・成分性・相互区別の証明を与える。式の数、計算時間、得られたサンプル数は完了条件ではない。

**計算上の注意:** 有限体の非零 minor は標数零での階数下界に使える。有限体の標本だけから標数零の既約分解や網羅性を結論しない。重み、消去、saturation、局所化で問題を小さくした理由を記録する。

## 優先3: 正規な実際の退化因子の定理 — 統合・監査完了

現行の [91頁独立原稿](../papers/normal_actual_divisor_recognition_ja.tex) と [証明・組版記録](../reports/NORMAL_ACTUAL_DIVISOR_PAPER_20261103.json) が以下の指定完了条件を満たす。旧版の数値証明書は変更のない検証済み記録を再利用した。以下は元の課題仕様として保持する。

**対象:** v22 `thm:normal` と、その v14–v21 への依存鎖。

**着手:** `audits/2026-10-30/geometry_strategy.txt` から依存関係を展開する。まず inactive line の排除と、E8側で残る八つの四次支持型の尽くしを読む。特定の支持 $C=L$ の仮定を、可約支持内の直線成分に誤用していないかを確認する。次に各証明書の入力が幾何学的な全場合を覆っているかを確認する。

**完了条件:** v14のGorenstein Jacobian判定からv22の結論まで、必要な仮定をすべて追える現行原稿があり、各外部定理の仮定が確認され、有限の排除計算は独立に検証できる。Pymのpurely elliptic分類との差をIntroductionで明記する。複素閉点の分類と、v14の全スキーム認識を混同しない。

**問題が見つかった場合:** 当該主張を `CLAIMS.json` で要修正にし、影響する依存先を明記する。無関係なv24結果や57の下界まで一括して撤回しない。逆に未監査の鎖を既刊定理と同じ確実性で使わない。

## 優先4: 主要な収録成果を公表可能な単位へまとめる

**第一候補:** 階数二 locus の8成分と、全空間の一般階数二成分6個。v14の独自の安定分岐を中心に、非安定・因子的零点の被覆を接続する。

**第二候補:** SNC locusを避ける成分と一般非被約構造。とくに被約な退化因子を持つ6追加成分について、MPSの40成分との区別とtoric degenerationの問いへの意味を説明する。

**着手:** 原稿を旧版の年代順に連結せず、定理・背景・証明・計算証拠の順に再構成する。文献上の新規性を一次資料で調べ、既刊結果、収録独自証明、新しい観察を区別する。日本語原稿はLuaLaTeX。

**完了条件:** 読者が再帰的に旧版を開かなくても主定理と仮定を読め、補助計算の必要性と検証方法を理解できる原稿と最小限の再現資料がある。原稿整理は数学的完全分類の達成とは別に記録する。

## 副課題L1: 交差次元の指定段階は解決

[共通27次元族](notes/AMBIENT_INTERSECTION_DIMENSION27_20261103_v62.md)と独立raw係数レビューにより、dim_Q(CA∩CB)=27、Mはその真閉部分で交差成分ではないと証明した。Qを通る交叉成分は唯一Nと解決した。個々のCA/CB/Nの全完成環もAcとして決定し、実交叉の非被約性・埋込素idealの存在を証明した。Q1での実intersection完成ideal・nilradical・全annihilatorと正確な二つのassociated primesも決定した。全射影被約交叉もN27∪R19と解決した。全域scheme交叉と大域primary構造は維持する。

## 副課題L2: $Q$ の有限な同変モデルと共鳴

**現在:** κ=1/4の追加方向の全次数積分に加え、全350式を保持する35変数の正確な有限同変模型を構成・独立検証した。指定の完了条件は満たした。全枝・交叉・非被約構造の解析を後続課題として続ける。重み6閉切断は厳密に二枝・一般長4と1だが、全35変数の分解ではない。

**従来の既知:** 横断34変数、独立な二次障害53本。$\kappa=1/4$ では横断接空間が35次元となる。二次障害の消滅は全高次の延長可能性ではない。

**着手:** `audits/2026-10-30/local_audit.txt` のトーラス同変・重み付き消去案を検討する。これはv24で実装・検証済みの完成局所環ではない。補空間・局所化する行列式・消去順序をまず証明し、元のJacobi式への代入で検証する。共鳴点は一般パラメータの分母を無条件に特殊化せず、別のsliceとして扱う。

**完了条件:** 実際の有限な残余イデアルとその正当性、または共鳴追加方向の積分可能性・障害について一つの明確な定理を得る。零点スキームの $\kappa=-1$ での変化と混同しない。

**優先度を上げる条件:** 新しい全成分を与える、残る大域成分を制御する定理につながる、または独立した概念的価値の高い変形定理になる場合。端点を通る全成分の尽くしを、そのまま全 $\mathcal P_4$ の尽くしと呼ばない。

## 保留するが維持する最終課題

- 全特殊点の $\mathrm{PGL}_5$ 軌道分類。
- 全成分交叉、埋込み成分、非被約スキーム構造。
- 二次・三次・四次超曲面上のPoisson構造と周囲への延長の全体的な対応。

これらは第一段階の進展に応じて再優先化する。完了していないものを黙って目標から消さない。

## 成果ごとの更新

各まとまった進展後に、何を証明したか、何を検算したか、何が残るかを `STATE.md` と機械可読の二つのJSONへ反映する。反例・訂正が出た場合も証拠と影響範囲を残す。時間だけを理由に成分数や確実性を上げない。
