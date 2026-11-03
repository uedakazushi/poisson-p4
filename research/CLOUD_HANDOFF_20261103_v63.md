# v63 Cloud 継続記録 — 2026年11月3日

研究継続中。これは停止指示ではない。全未解決課題を完了したとは主張しない。

開始SHAは `4c3f623fb698f9adc6077015b25e2325b11072af`。
v62保存は `e187b40de68bce534391f406f616940cddd045a7`、
`9e82f3c5120380028990cdb58b07318866c21f4b`、
`66a80ab53dfcf0ff842c67140320d8da0226039d`。全て通常push済み。
このv63保存のSHAはGit履歴で識別する。

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

## 未完了と次手

一般アフィン未知曲面・曲線の除去、一般五点の尽くし、特殊一次像上の未知成分、
既知全成分の定義prime、全大域根基一致、成分交叉と非被約構造、全軌道分類を維持する。
細分化した次手は `research/NEXT_TASKS.md`。下界57、上界null。

進行中の追加査読：基底Zのpinched-surface完成環、個々のCA/CB/Nの全pinched完成環、
rank≤2上の被約交叉。producerの証拠があっても独立監査が終わるまでは本保存の定理へ加算しない。
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
