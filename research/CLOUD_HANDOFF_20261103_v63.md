# v63 Cloud 継続記録 — 2026年11月3日

研究継続中。これは停止指示ではない。全未解決課題を完了したとは主張しない。

開始SHAは `4c3f623fb698f9adc6077015b25e2325b11072af`。
v62保存は `e187b40de68bce534391f406f616940cddd045a7`、
`9e82f3c5120380028990cdb58b07318866c21f4b`、
`66a80ab53dfcf0ff842c67140320d8da0226039d`。全て通常push済み。
このv63保存のSHAはGit履歴で識別する。

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
