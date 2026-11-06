# 現在地 — 分類に必要な情報だけ

2026-11-07 / v64-classification-first。既存数学基準はa00a541b、方針更新はe1d36387。
現行計画はPLAN.json。旧STATEは証拠の入口であり実行キューではない。
**下界58、一般rank4は52。全成分数の上界・一般五点尽くし・全分類は未証明。**
一般rank2の全空間成分は、監査した支持証明と既存成分証明を用いて六族に尽くされる。

## 現在使える主線

| 義務 | 確認と残る穴 |
|---|---|
| C0 | S61→primitive betaの橋det36/625、六切片の優越性をQQで再構成。55消去・Bianchiは原証明に依存 |
| G_direction | 完了。全K16の有理模型の方向incidenceでrank[A d]=4が全v≠0に成立。交換子det21866 mod32003、全45中心と45×3frameを独立再構成 |
| G | 未解決。rank3枝は排除でき、残るのはrank[A b]=rank[A b d]=4上の全F5。有限48/144被覆で追加超越変数なしの忠実chartとGalois根基復元が可能だが、unit未証明 |
| E0 | 準備完了。既知49proper閉像を区別。一般五点の底分母はG28を因子に持つ必要があるが十分ではない |
| E_zero | 旧Q45_22・MPS・既知成分定理を入力にbeta=0の全被約台を既知52閉包で被覆 |
| E_poly | 完了。SNC最高部を持つ二次多項式代表の全Poisson低次項を被覆。十正常形のkappa非零部分と実際のPGL飽和・閉包。全Q閉包の逆像は主張しない |
| A_rank2 | 完了。安定証明の全読解・独立査読、CM/LPT一次本文、spinor二実装と全Aff32変形の再現で全rank2支持を六族とrank4閉包へ接続 |
| A | rank4追加11族・NEW25の成分性と相違の全依存監査が残る |
| E/F | 全例外底・全分類は未完了。未知特殊成分を既知49像の和に限定しない |
| M | 現主原稿へ有限被覆、rank3除外・直接splitter、rank2支持、低次項SNC被覆を統合。主分類未証明を維持 |

## Gで新たに閉じた部分

NORMAL_DIRECTION_INCIDENCE.mdと独立reviewを読む。元134式と同値なF5と
R=X(z⊗n)+D1n+D2(n,n)を使う。L X=I、T_h=−(L D1)_hについて
[T0,T1]の42角行列式が非零である。共通固有ベクトルが存在しないのでrank3枝を排除する。
正当な全K16有理模型との特殊化は普遍三平面と明示中心で照合済みであり、
固定fiberの空性を一般空性へ移したものではない。

さらにM_h=z_hI−T_h、C=[T0,T1]とすると
U=C^(-1)[−M1,M0,0]Lはz次数1の多項式左逆である。
H=(1/2)Jac_n D2、delta=det(I+UH)からdelta n=0、delta≡1 mod(n)を得て、
delta²=delta、(n)=(1−delta)A、A=E5×Bが直接従う。
旧quintic表示・精密nil構造はこのsplitter証明に不要。**B=0はなお未証明。**

FINITE_SEPARABLE_NORMAL_CHART.mdはK5の有限48/144被覆を使う。
42独立係数によるell一つでB→B_K'[1/ell]が忠実平坦となり、係数体の超越次数は5のまま。
全259式の係数引戻しを別実装で検査した。全分岐補集合と元分母はEへ残す。
旧全K16固定etaをK5の同名形式へコピーすることは依然として許されない。

## 次の数学と計算限界

次はrank4方向の全F5が解を持たないこと、または未知点・族を実際に示す。
Cramer表示の全minor chartと無限遠・除外方向を保持する。
有限被覆chartや未展開の小行列積は使えるが、表示を作るだけでは支持問題は閉じない。

保存中心の直接GBはSymPy二回、native std、G6/slimgb4GiB/12GiBの五回とも未完了。
同じ直接GBの反復は停止中。今回のnormal-onlyコスト検査では元29,301項に対し、
五つの四次式を入れる前に331,008項となったため、追加GBを実行しなかった。
Singular4.4.1はwork/native-cas/root、環境のpackage/hashは既存NATIVE_ENV報告にある。
同じ固定betaでG28≠0をQQで確認済みだが、他の特殊像の除外や一般性は別問題である。

## 依存と検証の射程

rank2支持の証明はRANK_TWO_SUPPORT.md、低次項被覆はLOWER_DEGREE_SNC_COVERAGE.md。
CM/LPT/MPSの外部定理自体を再証明したとはしない。全旧証明書・全134生成を再実行していない。
新しい行列、全normal/center対応、全低次項重み、必要なrank2恒等式を正確に確認した。
計画checkerと制御テストのPASSは数学の全分類認証ではない。
主原稿PDFは現TeXから12頁で再組版した。2passで未解決参照・overfullなし、1/5/9頁を目視確認。
source/PDF hashと環境はCLASSIFICATION_MAIN_MANUSCRIPT_BUILD_20261107.json。
計画制御の23テスト、strict check、構文とdiff検査を通過した。
Libraryへのアップロード、完全依存bundle生成は未実行。
B18 H0/T2/nil指数、H絶対既約性、旧minor farmは休止のまま。
