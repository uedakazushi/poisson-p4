# 現在地 — 分類に必要な情報だけ

2026-11-06 / v64-classification-first。既存数学基準はa00a541b、方針更新はe1d36387。
現行計画はPLAN.json。旧STATEは証拠の入口であり実行キューではない。
**下界58、一般rank4は52。上界・一般五点尽くし・全分類は未証明。**

## 今回確認した主線

| 義務 | 確認と残る穴 |
|---|---|
| C0 | 一巡を完了。元S61→primitive betaの橋det36/625、標準六切片det8428695240、疎六源det137438953472/996451875をQQで再構成・独立検算。generic5にはSp4とGmが必要。55消去・Bianchiは原証明に依存 |
| G | 既存splitterを仮定すると、超越tauのell=Σtau^j n_jによる一つの忠実チャートで全残余の空性を検出できる。元134式+1-u ell、任意で259法線式の単位証明が正確な目標。単位所属は未証明。別途全K16上で固定etaの42軌道spanを独立確認（minor27110）。K5の同名形式とは異なり輸送は未証明 |
| E0 | 準備を完了。52写像の全hash・源inverse・像次元・補集合を区別。一般五点開D(f)なら既知49閉像上でf=0、特にG28がfを割る。既存成分・S61対応・G28定理に依存 |
| E_zero | 旧Q45_22・MPS・既知成分定理を入力にbeta=0の全被約台を既知52閉包で被覆。22族の全包含、10SNC極限・2全9パラメータ恒等式・5稠密軌道を別実装で確認 |
| E | 一般証明の底分母f_Gは未確定。G28、判別式、beta零、未知fiber局所化の残部を保持。既知49像の和が全特殊底とはしない |
| A | 部分監査。rank≤2支持を六全空間成分とclosure(rank4)で覆えば十分。六族の相違・58下界の依存を整理。全独立監査は未完了 |
| M | 主原稿に座標橋、忠実チャート、必要底因子、弱いrank2接続を統合。分類主定理の未証明表示を維持 |

証明・独立査読は同フォルダのC0、GENERIC_NORMAL_CHART、EXCEPTIONAL_BASE、KNOWN_COMPONENT_AUDIT。
小行列はreports/CLASSIFICATION_SLICE_BRIDGE_20261106.json。
E0のbeta零制限は全134行のうち二次70行、零64行、未知数45個を保持し、別のQQ代入で全行一致を確認した。
NEW25全1350祖先節点とDPM全1263節点も独立照合。DPMのsqrt3は自由変数でなく一つの複素埋込みを固定した係数である。

## Gの実計算と次の一手

固定保存中心・GF(32003)(tau)の全134式+逆元式をSymPy F5bで診断した。
55 CPU秒・1.5GiBで未完了（SIGXCPU）、unit・反例・一般被覆の証拠は得ていない。
reports/CLASSIFICATION_GENERIC_CHART_PROBE_20261106/に元入力hash・実行producer bytes・停止状態を保存した。
再生成可能な二つの圧縮診断入力は自動送信審査の対象外とし、work/regenerable_diagnostic_inputsへローカル保全した。数学の証明書ではない。
tauは超越のままだがbetaと標数は固定しているため、成功してもQQへの移送は別義務だった。

次は同じ点・次数・wrapperの反復ではなく、native exact backendか有限syzygyの構造を選ぶ。
通常のexec sandboxではproxy接続が拒否されたが、標準のnetwork.enabled追加権限で公開Debianパッケージに接続できた。
Singular4.4.1と必要ライブラリをwork/native-cas/rootへ展開し、optimized moduleの読込みとstd smokeを確認した。
ソースとhashはreports/CLASSIFICATION_NATIVE_ENV_20261106.json。M2は未導入。システムインストール・proxy迂回は行っていない。
固定normal0のGF32003診断も305wall秒で未完了。直接SymPyGBを二回続けてもunitは得られず、このbackendの反復を停止した。
native stdも315wall秒で未完了。全9257係数のQQ恒等式で元ideal内と確認したG6（固定点では790項）を加えてslimgbを使ったが、4GiBでは33秒、同じ式の12GiB評価でも約65秒で実メモリ上限に達した。
最大16GiB環境で12GiBと2GiB以上のreserveを事前確認した。いずれもbasis・unit・nonunitは未決。NATIVE_FIXED_FIBER_SUPPORT.mdを参照し、同じ直接GBを自動反復しない。
同じ固定T01のbetaではG28=671334843179656461002307382635849/9039207968、法32003で25057と正確に確認した。既知NEW25像だけは外れる。他の特殊像や一般性は未証明。
別の構造的還元では全K16上の固定eta一つの135式へ移せるが、K5へ群変数を落とす橋は未証明。
K5+tauを保持するか、全K16で係数を扱う。どちらのunitも未証明でありGをdoneにしない。

## 維持する区別

* full134の一般根基包含と、選んだ切断・minor・jetは異なる。
* K[n42]上の有限5生成はK上有限を意味しない。
* actual42のdelta_E次数28と特殊底NEW25の像G28は別物。
* 一つの忠実チャートはK(tau)上の全点が同じ開に入るという意味ではない。
* 源の分母、底の分母、未知qの逆元を区別し、全補集合を残す。
* rank2の支持被覆は既存の安定接層証明と一次文献に依存する。

B18 H0/T2/nil指数、H絶対既約性、旧軸minor farmは休止のまま。
独立の分類への不可欠な依存がなければ再開しない。

## 検証と保存の射程

計画制御の23テスト、strict check、job hash、構文、diff検査は通過。
これらのPASSは分類の証明ではない。原134式の全再生成や58全証明の再実行はしていない。
不足していた日本語TeX二パッケージをwork/native-tex/rootへ展開し、主原稿PDFを現TeXから8頁で再生成した。
2回の完全組版で参照を解決し、最終passに未解決参照・overfullはなく、1頁と7頁を目視した。
再現環境・package/PDF/source hashはreports/CLASSIFICATION_MAIN_MANUSCRIPT_BUILD_20261106.json。
Libraryへのアップロードや完全依存bundle生成は今回未実行。
