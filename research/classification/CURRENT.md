# 現在地 — 分類に必要な情報だけ

2026-11-07 / v64-classification-first。既存数学基準はa00a541b、現行計画はPLAN.json。
**下界58、一般rank4は52。全成分数の上界・一般五点尽くし・全分類は未証明。**
一般rank2の全空間成分は六族に尽くされる。現在の58成分の成分性・相違の依存監査は完了。

## 現在使える主線

| 義務 | 確認と残る穴 |
|---|---|
| C0 | S61→primitive betaの橋と優越切片をQQで確認。55消去・Bianchiは原証明に依存 |
| G_direction | 完了。交換子による階数三枝排除と全K16への橋を独立再検証。これは旧20261103の方向graphにも含まれる結果 |
| G | 未解決。全134式の未知残余Bの空性が必要。次は正規被約退化因子と可解unimodular零点のincidenceを調べる |
| E0 | 既知49proper閉像を区別。一般五点の底分母はG28を因子に持つ必要があるが十分ではない |
| E_zero | 旧Q45_22・MPS・既知成分定理を入力にbeta=0の全被約台を既知52閉包で被覆 |
| E_poly | SNC最高部を持つ二次多項式代表の全Poisson低次項を被覆。全Q閉包の逆像は主張しない |
| E_quadric | 完了。実退化因子H+2Q、rank Q>=4の全点を既知閉包で被覆。rank3横断配置も余剰逆形式を含めて被覆 |
| E_quadric3 | 未解決。rank3の頂点非包含・接触配置。全19/25次元の逆形式とPfaffian条件を保持 |
| A_rank2 | 完了。全安定・非安定・因子付き支持、CM/LPT、六族とrank4閉包への接続 |
| A | 完了。既知58の成分性・次元・相違を監査。追加11族の真の横断切片とNEW25八初期式による支持上界を接続 |
| E/F | 他の退化因子型、全例外底・全分類は未完了。未知特殊成分を既知49像の和に限定しない |
| M | 支持被覆と成分性監査を同じ主TeXへ統合。分類主定理未証明を維持 |

## 二重二次曲面の全点被覆

DOUBLE_QUADRIC_SUPPORT.md、TANGENT_DOUBLE_QUADRIC_MPS_COVERAGE.mdと独立証明書を読む。
全90係数の閉逆形式をQQで確認し、Pfaffian条件を別に全て課している。
rank5では既存の一軌道認識を接続。rank4では頂点回避配置は空、非接触配置はbeta零へ帰着。
接触配置は全点でisotropic五変数源とanisotropic一軌道に分かれる。
前者はT18の実源への恒等式・軌道極限・全周囲接次元23からT18へ入る。
後者は全パラメータの実projective Poisson pencilと二つの軌道極限でMPS40の和に入る。
新成分を数えたのではなく、指定層の全支持被覆を閉じた。

rank3横断配置Q=1+u²+v²は閉形式20次元で、affine原始形式19次元以外の一クラスも必要。
RANK_THREE_QUADRIC_SUPPORT.mdの全点変形と独立な全140×90整数行列検証により、
全点をbeta零飽和閉包で覆う。残る二配置はQ=u²+v²+z²とQ=w+u²。
後者の六追加クラスu dr−r duを落とさない。

## Gの現在の数学と中止した経路

NORMAL_DIRECTION_INCIDENCE.mdの全134式と同値なF5+R表示は保持する。
交換子の非零minor21866とcross19216、全45中心と45×3frameの橋は独立確認済み。
多項式左逆Uからdelta=det(I+UH)を作ればA=E5×Bが直接得られるが、B=0は従わない。
有限48/144被覆による忠実chartもunit未証明。全底分母の補集合はEに残る。

直接GBは五回未完了で停止中。normal-only展開も四次式前に331008項へ密化した。
低次構造だけの経路は、同じ可逆交換子・二次cross・平面五点・有限性を持ちながら
第六の孤立点が残る反例により限界を確認した。旧compact graphや境界次元還元を
新しい尽くしとして繰り返さない。

次の入口はDOMINANT_COMPONENT_GEOMETRY.md。
未知のbeta優越全空間成分の射影次元は25,26,27に限られる。
一般退化因子が正規被約なら、Pymの定理により、T39以外ではunimodular零点を持つ。
そこで退化因子の重複度は3以上、三次初項は一次Lie–Poisson構造のCasimirである。
非可解4次元一次Lie代数はsl2⊕Cであり、Whiteheadによる全形式線形化でrank4と矛盾する。
従って残る一次部は可解unimodularである。この局所排除と外部定理の適用は独立確認済み。

未証明のG-normalは、T39を除く正規被約・unimodular零点incidenceのS61一次像が
次元15以下であること。特異点のLie一次部とrank4基点のprimitive betaを混同しない。
この分岐が閉じても非被約・被約非正規退化因子と特殊betaは残る。

## 依存と検証の射程

KNOWN_COMPONENT_AUDIT.mdが58下界の入口。追加11族は全周囲の真の形式切片が零次元なら十分で、
精密なnil指数や全普遍jetは追加前提にしない。NEW25では八つの真正初期式と三次恒等式・
Sylvester行列式で横断次元<=3を再証明した。元の全350式模型と成分証明書は既存入力である。
CM/LPT/MPS/Pymの外部定理自体、全旧証明書・全134生成を今回再証明・再実行したとはしない。

新しい逆形式の全行列、実Poisson変形、軌道重み、必要な接空間minorを独立に検証した。
共通double-quadric runnerの初回producer全bytesはCLASSIFICATION_DOUBLE_QUADRIC_SOURCE_20261107.jsonに保存。
古いproducer/hashと証明書を変更せず、現在runnerにはanisotropic modeを追加している。
計画checkerと制御テストのPASSは全分類の数学的認証ではない。
Libraryへのアップロード、完全依存bundle生成は未実行。B18精密構造・旧minor farmは休止。

主原稿PDFは現TeXから16頁で再組版し、最終二passに未解決参照・overfullなし。
1/13/15頁を目視確認。build報告はCLASSIFICATION_MAIN_MANUSCRIPT_BUILD_20261107.json。
計画制御23テスト、strict check、変更codeの構文、diff検査を通過した。
