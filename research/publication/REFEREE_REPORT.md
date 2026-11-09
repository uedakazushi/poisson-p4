# 公表前の独立査読 — 現行結論

2026-11-09（Asia/Tokyo）。構築段階の147成分という結論は査読の仮定にしない。
本報告の「確認」は原証明を独立に追った範囲を意味し、全有限証明書の再実行を意味しない。
**数学的判定：指定した高リスクの証明義務を充足し、公表候補版として可。**
A1–A3/B1/B2は以下の通常の証明・正確な一次文献・有限証明書の射程で閉じている。
65構文パスの読解と実入力の修復も完了した。外部査読の採択や新規性・優先権の確定を意味しない。
現行全bundleの隔離検証・保存の実行状態は `CURRENT.md` と配布receiptで別に確認する。

## A1 — 被約固定五次の線形源：一般論を独立確認、三例を全空間で再計算

対象：`papers/normalization_residue_reconstruction.tex` の
`lem:quintic-ambient-lift`、`thm:reduced-fibre-linear` と直後の二系。
仮定は複素P4、被約五次F、正則な大域分子、closedなα/F。

独立導出：Λ²ΩDの商写像の双対と
0→Ω²→Λ²T→Λ²T|D→0、H⁰Ω²=H¹Ω²=0から任意の二重導分のliftは一意。
JacobiatorがD上で零ならH⁰(Λ³T(−5))=H⁰Ω¹=0からambientでも零。
各素因子の一般平滑点で単極closed二形式をdz/z∧η+γと書くとdDη=0、
Πは接し、Π²はその素因子で割れる。squarefree性と次数五からΠ²/2=qF。
三次元のFrobenius条件、被約Dへの延長、先のliftから、q=0も含め全V_FがPoisson。
線形性のpolarizationが互換性を与える。rank四はq≠0だけである。
逆形式と同じ線形双対性が全rank四fiberの逆を与える。代数閉拡大でも同じ証明が成り立つ。

最初に疑った点：特異・交叉因子、非零核だがq恒等零、相対kernelのrank jump、
因子像の包含から成分包含への飛躍。既存150係数adapterを使わず、C5の
**350個の三次二形式係数**にEuler水平性とFdα−dF∧α=0を課した別実装
`code/referee/check_reduced_fixed_divisor.py`でQQ上の全核を計算。
`reports/referee_reduced_fixed_divisor.json`に全基底とqを保存した。
五座標H積は(dim V,rank q)=(6,6)、滑らかな二次×楕円三次coneは(1,0)、
特異既約free平面五次は(10,4)。全基底自身と全基底対の和でJacobiとPf=qFを恒等検査した。
有限三例は一般定理の証明の代用ではない。

境界：q≠0のβとのpencilは全q=0点を同じ固定因子のrank四閉包へ送る。
q恒等零の例にはrank四fiberが無い。一定rankの**被約**底では可逆minorと
Schur補行列の消滅から核がvector bundleになる。この底の被約性を系へ明記した。
反例候補tx=0⊂A³(t,x,z)では全fiberが既約でも{x=0}と{t=0}は別成分で、
底像が包含されても全空間の成分包含は成り立たない。
同じ底像を支配する二成分は幾何学的一般fiberの既約性により一つなので、像の**一致禁止**は妥当。
非被約h⁵の既存例もJac(y,z,w)=−2を直接確認でき、仮定を除けない。
依存先：正規化の全降下、rank二境界、標識忘却と極大像の重複排除。

## A2 — 階数二支持：安定分岐と外部分岐の論証を独立に読了

対象：`papers/rank_two_support_reconstruction.tex`。
仮定は飽和反射的rank二接層F⊂TP4、零因子次数e、det F=O(e)、d+e=2。
安定分岐ではe=0、H⁰F=0、canonical cone liftはdiv=0でprimitive Plücker係数を持つ。

独立に追った核心：Euler商からF(1)の五切断と線形Casimirの核を復元。
一般切断の余核はI_Z(2)で、depth補題によりZはunmixed。局所自由・LCMを仮定しない。
二quadricsのCIと安定性からdeg Z=2,3,4。次数四はCIのKoszul分解、
次数三は残余平面とGorenstein双対によりHilbert--Burch分解を得る。
前者の二線形CasimirはRat113へ、後者の四切断は線形一形式のcurl pencilを使いRat122へ戻る。
次数二ではrank零surfaceがあれば全Plückerがそのidealの二乗に入り、
binary二次写像P1→Gr(2,5)の分裂からH⁰F≠0となる矛盾。
従ってBertiniとunmixed性を正しく使えてZは二平面の被約和。
非退化性から両平面の交わりは一点であり、この非LCM点も除外していない。

二平面のExt¹は四一次式の六Koszul syzygyを四行で割った二次元。
反射性が両拡大成分の非零性を強制し、二ブロックの尺度変更で正規化できる。
完全五切断の十quadricsはSym²C4*の基底で、OGr(2,5)=P3の全像を与える。
元の五切断の像Wを完全切断空間と混同せず、rank W≤3、5、4の二軌道に分ける。
div=0の最後の線形系はrank5と非等方核を排除し、等方核の全八係数では
dt∧dQ1∧dQ2=−ιPΩが明示されてRat122になる。この最後の表示は通常の多項式恒等式として確認した。

一次資料を照合：Corrêa–Muniz arXiv:2207.12880v5 本文Theorem4.1の安定分岐は
代数的可積分性までであり、Rat113/122全認識は上の内部証明が担う。
同Theoremの非安定分岐とProposition3.9の余次元二固定加法作用は本文の使用範囲に合う。
Loray–Pereira–Touzet arXiv:1112.3871v3 Theorem6.2（印刷17–18頁）の
degree1・任意余次元の二分岐をe=1へ適用する。引用の印刷添字を機械代入していない。
toricは互換pencil、Aff32は全b,t多項式変形でrank四閉包へ送る。
この閉包の置換はrank四の有限被覆を先に得た後で行うため循環しない。
六成分のambient上下界は全Jacobi微分と族微分の非零minorによる。
元六行列の前構築段階の再実行は既存証拠であり、今回新たに独立再計算したと数えない。
依存先：全rank二点の支持、六純rank二成分、最終有限閉被覆。

## A3 — 最後の五像：全源、proper条件と非包含の論理を確認

対象：主稿 `source:last-five-explicit`、`prop:five-last-full-components`、
`papers/maximal_image_synthesis.tex`。
全源はGとDの全PGL軌道、三つのA6磁気源の全PGL飽和。後三者の六係数を残して、
`check_universal_sources.py`で全Jacobiとactual五次を別実装で確認した。
Pfは順に−4hF4、非零定数×hxCi、h²C3で、磁気に依存しない。
free平面五次のCramer/Saito源では低次syzygyが無く、全混合列は定数GL2 frame、
残りが全R6となるので、標本軌道に置き換える飛躍は無い。
実shearの余核は第一二で(I,J)、第三で(I,J)の別の二次元商。
単元で一parameter化する前のA6全体を保持し、零磁気の個別外性は主張しない。

T36：quadric pencilのGrassmannianからの写像にbase点は無く、零を含むconeは閉。
点の零次・一次jetと二次jetの条件を課し、properな点標識を忘れるので全閉包へ延びる。
対象の全order二零点で非零horizontal三次一形式が既約cubic/quarticを不変にする。
そのcurveがquadric比のlevelなら次数二式で割れる矛盾。T38の零surface条件も
普遍零schemeのproper射とfiber次元の上半連続性で全閉包に有効である。

T37：主稿のHilbert--Burch分解からbaseはintegral・degree5。
bF=g1g2−q³、CI(F,g1−κg2)のdegree15にbaseが少なくとも二重に入り、
一般の葉曲面はdegree≤5。whole境界ではactual Dがsmoothなrank二点をsectionで標識し、
一枚のgeneric葉のHilbert polynomialを固定してproper Hilbert射を使う。
flat商のtorsionfree性が包含とHamiltonian不変性を延長し、char0のradical/minimal-prime保存により
標識点がembedded pointだけへ逃げることを排除する。重複因子の被約台だけへは適用しない。

別導出 `check_leaf_obstructions.py`：Ga9では二一般線形式でx,zを1,y,w(y)の線形結合にし、
Casimirのy⁹係数−1/36から葉degree9。TRIPLEでは
y=−(x³+b1x+b0)/(3x+b2−3/2)をCasimirへ代入し、分母と互いに素なdegree6式を得た。
三磁気源のannihilatorを独立に縮約し、留数は
(J,I/2,−I/2)、((I+J)/2,(J−I)/2)、−3(λ+1)。有限coverでも非零留数は消えず、
一般磁気では代数葉が存在しない。GのH上にはdegree3葉があるという反例候補も調べたが、
単純quartic上のdegree9がT37の必要条件に反するため、非包含を覆さない。

prior142比較：6 rank二、18低次元は閉rank/次元、同次元親は包含⇒等号を使って
稠密actual配置・内在的cofactorで比較。高次元親にはproper factor像、T21の二H接触旗、
T36/T38、上のT37、全S61グラフ又は全周囲接空間上界を使う。
Ga9 cuspidal24のmarked DVR議論では有限正規化のspecial fiberを最初からflat/smoothとせず、
CM・integrality・算術種数零からP1族を得る。DへのD5変形は相異なる四重根と三重根を
次数五式が持てないこととz=0の全場合分けで除く。
全IDのpartitionと根拠は主稿および `GA9_LINEAR_W_PRIOR142_WHOLE_COMPONENT_EXCLUSIONS.md`、
`TWO_LINES_PLANE_CUBIC_PRIOR142_AND_COMPONENTHOOD_SEPARATE_REVIEW.md`、
`CONSTANT_MODULAR_DOUBLE_PRIOR142_LE23_DISTINCTNESS.md` にあり、表の件数自体を証明としない。
有限全被覆で未知成分dim≤23を先に得てから五像を極大化する矢印は正しい。
依存先：五成分の存在・相違、141+6という最終極大閉像。

## B1 — 全51一般源：普遍恒等式と保存演算列を区別して確認

`reports/referee_universal_sources.json`：T01–T38と追加十一族は正規化前の全有理式、
T39/T40は保存DAG自体を全16入力で記号復号し、150係数をadapter無しで再構成して
全projective JacobiをQQ上（DPMはr²−3剰余）で検査した。全51で恒等零。
parameterだけの分母を払って多項式計算するので、分母零での正則性は主張しない。

`reports/referee_serialized_maps.json`：残49DAGは読了した元constructorを再実行し、
**214,525節点、全入力・61出力・正規化metadataが完全一致**。
これはserialize再現であってconstructorの独立実装ではない。
constructorの一般証明は、M P(M⁻¹x) Mᵀ、FCFᵀ=J、
Sp4大cell [[A,AS],[TA,TAS+A⁻ᵀ]]、Euler24商、定数可逆応答行列によるperspective四方向の除去。
いずれも開上の実射影座標変更で、Poisson性を保つ。
MPSでは対角torusが一般基点を(1:1:1:1:1)へ動かしperiod/edge源を保つ。
追加十一では一般基点四係数とscaleを保持するので、同じ議論をtorus無しで適用できる。
全symplectic frameの大cellを残したため一般PGL飽和を覆い、slice像は稠密である。

源はaffine空間の非空開。forced-zeroはperiod線形核に既に組み込み、
既存の正確な評価点で全inverse分母が非零であることは非空性の証拠として再利用する。
その標本を普遍Jacobiの証拠にはしない。DPMのrは一つの複素定数であり超越parameterでない。
反対埋込みでも同じmarked3A2 pairの全Saito源、磁気商177±100√3≠0の同一軌道へ戻る
（主稿 `prop:marked-quartic-A0-DPM-positive`）。根を別成分と数えない。
T39/T40は全A16グラフなので別のSp4飽和を省略したという問題は無い。

一次資料 MPS arXiv:2010.08692v2 のTheorem1.4は**SNC locusと交わる40成分**のみ。
Theorem6.16はX≤2上の局所成分、Corollary6.23の条件2はExample6.24で射影空間上成立し、
同Example末尾がn=2で全40の大域相違を確認する。この入力と実period/edge cocycleの対応で
成分名を認識し、像の次元一致だけでMPS成分と断定していない。
依存先：カタログの一般族、T39/T40を使う非包含、S61表示。

## B2 — 全点被覆から極大像までの必須接続を確認

再読した最小の接続は次の通り。元の有限認識を再利用する節点を明記する。

|入口|必要な全点定理と依存の使い方|
|---|---|
|actual既約normal五次|別稿 `thm:classification`。正規性だけをPymのpurely elliptic性へ置換せず、不活性・有理支持排除の全有限証明を入力として使う|
|actual既約nonnormal五次|主稿 `prop:irreducible-nonnormal-all-pointwise`。S2と正規化商、Σδ deg S=6−gからdegree1–6、rank0/2へ。planeは全35/30/9源、symbol零・scalar・leading零も別枝。nonplaneの全表の修正版認識を再利用|
|Q2+C3|共有稿 `thm:quadric-cofactor-reconstruct` と退化二次の全31代表。radicalのtrace零とSylvester余核を残し、各全支持証明へ戻す|
|Hを持つ被約因子|主稿の全因子分割表と `prop:HQQ-all-pointwise-1125`。simple HでPf最高零を仮定せず、全normal場・plane磁気を保持|
|非被約因子|重複H無しはHQ²だけ。重複Hは `thm:marked-affine` のdeg Pf≤3から全Jordan/translation余核へ。定数modularのdouble/tripleとunimodular源も含む|
|最後の残項|P*の真completed T19曲線、TRIPLEのproduct Poissonからの全e二次族とactual Ga9 affine逆。1−ex/3を反転せず、e=3,12も保持|
|rank二|A2の全支持へ、上のrank四有限閉和を得た**後**で閉包項を代入|
|極大性|既存142の全周囲minor/真正形式切片等とA3。有限既約閉和の補題は最後にだけ使う|

原稿の全pointwise命題の合成に循環は見つからなかった。共通復元の線形源から
各有限認識・positive membershipを飛ばす矢印も無い。
ただし本査読で223 primary claimsの全局所有限排除を別実装で再証明したとはしない。
全周囲行列・全係数恒等式等の収録証明書を、上の限定された節点で再利用する。
今回新たに正確再計算したのはA1/A3/B1の各範囲である。

最終依存監査は3539論理ファイル、16840記録hash、未解決の記録入力0。
静的走査が出した65構文パス・206組を、実際のASTと出典SHAに照らして読解した。
実入力8、書き出し先51、native生成出力4、事前生成入力1、archive locator1であり、未分類は0。
このうちv21の二runnerは旧bundleに無かったので現行依存へ追加した。
旧八支持runnerの文献PDFは算術に不要なのに必須readだったため任意のprovenanceへ修正し、
PDF無しの隔離実行で元の全算術・支持被覆fieldとの一致を確認した。
既存証明書が要求する元runner一版を追加保持し、凍結入力は実数で97件となった。

一つの本当の事前生成入力はweight6のfresh比較reportである。既収録producerを隔離実行し、
時間field以外の全fieldが旧reportと一致した。四つのnative出力はrunnerが生成後に読むものであり、
事前の係数入力ではない。native jobの実行には所定のCASが必要で、今回全てを再実行したとはしない。
追加したstdlib checkerでは、八つのtree模型の全350係数、Pfaffian、原始関数等249項目を再計算した。
これは旧3229個の局所有限排除を全再実行したとの主張ではない。

`static_path_review.json` と `check_static_paths.py` がこの有限読解を固定し、出典が変われば判定を無効にする。
RUNTIME/CHECKPOINTという名前だけで欠落を免除する規則は撤去した。
二つの既知の失敗実行metadataだけは、出典と診断hash及び置換された係数証明書を照合して除外する。
未知のruntime名readを出力と誤認しない反例テストも行った。
任意の生成コードの動的全readを網羅したとは主張しないが、検出された必須依存は閉じている。
exact入力materializerの黙示skipも例外へ訂正し、既存二bindingの成功と欠落bindingの拒否を確認済み。

## C1/C2 — 再現・編集と証拠の射程

C1では保存時receiptのSHA256と一致するローカル五volumeを用い、
3579ファイルと二相対リンクの隔離復元、全SHA検査、四PDF各二回組版を実施した。
全四PDFの抽出テキストは同梱PDFと一致した。
Libraryの五volumeの保存参照・MIME・サイズは取得済みで、保存記録は別receiptにある。
**Libraryからのバイナリ再取得はユーザーの指示により検証義務から除外した。**
今後は現行Gitまたは既存ローカル入力の実復元と数学的証明の査読に集中する。
65パスを修復した現行v74では、補足overlayに依存しない全bundleを再生成する。
実行・保存の範囲は `CURRENT.md` と `referee_runtime_input_validation.json`、配布receiptに記録する。

C2では四稿の主定理・導入・依存順を照合。主稿のMPS40の射程と安定rank二の内部復元を明記し、
一般S61の五点定理を主分類とは独立な補論と明記した。証明を削ってページ数を減らす編集はしていない。
Pym arXiv:1507.05668 Theorem1.3/Definition4.2/Theorem4.5、Kaledin math/0310173
Theorem1.5/Corollary1.6、Gualtieri–Pym 1203.4293v1 Proposition3.11も一次本文で照合した。
Pym総説のPMC入口は今回captchaだったため、上の原論文を直接確認した。
新規性・優先権の全面確定、外部投稿・公開は行っていない。

**推奨：公表候補版を確定し、検証した現行全bundleと四稿を保存する。**
数学的な未解決補題又は必要な算術入力の欠落は、指定した査読範囲では残っていない。
既存有限認識を再利用した範囲と独立再計算した範囲を保ったまま、本査読を閉じる。
全147計算の独立自動認証、全PGL軌道・全scheme厚み、文献上の優先権の確定は結論に加えない。
新しい証拠が得られた場合は、その依存先だけを再開する。外部投稿・公開は別の明示指示を要する。
