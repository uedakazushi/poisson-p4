# 現在地 — 公表前査読の結果と残る再現義務

2026-11-09（Asia/Tokyo）/ v73-prepublication-review

**査読全体は未完。** 147成分（一般rank四141、rank二6）の分類稿について、
高リスク定理の独立導出・全パラメータ検算・最小依存の読解と原稿訂正を行った。
数学的反例は今回見つかっていない。以前はLibrary上のZIPのバイナリ再取得をC1の未達要件としていたが、これはユーザーの指示により査読条件から除外する。
公表準備の最終判定は独立査読の数学的射程と現行入力による再現証拠から行い、未判定とする。詳細は `research/publication/REFEREE_REPORT.md`。

## 今回確認したこと

|対象|今回の根拠と範囲|
|---|---|
|A1 固定被約五次|Ω²とΩ¹の消滅からlift/Jacobiを独立導出。350個の三次二形式係数から三つの全核を別実装で計算し、全基底対の互換性・Pf=qFを検証。dim V=1でもq恒等零の例を含む。相対系へ底の被約性を明記|
|A2 rank二|飽和反射層・五切断・非LCM二平面拡大・spinor・非安定/因子付きの全枝を読解。八係数のJacobian恒等式を独立再計算。CM/LPTの一次本文の仮定と比較|
|A3 最後の五像|全五族のJacobi/Pfと全磁気六係数を別実装で確認。proper flags、whole T37の標識付きHilbert極限、prior142の比較論理を読解。Ga9/TRIPLEの葉degree9/6と三磁気族の非零留数を別導出。零磁気の個別非所属へ外挿しない|
|B1 全51一般族|全パラメータprojective Jacobiを独立計算。T39/T40は保存DAGの全16入力で直接計算。残49演算列は元constructor再生成の214525節点・全入出力・metadataが一致。一般正規化・稠密性とMPS40の正確な射程を別に読解|
|B2 全支持の合成|全因子/Jordan/translation余核→有限閉和→rank二閉包代入→極大像の順を確認。元の有限認識・成分性証明書は再利用し、全223主張を今回独立再証明したとはしない|
|C2 原稿|MPS既知40と総数147の違い、安定rank二の内部復元、独立S61補論の位置を導入へ明記。四稿二回組版、未解決/重複参照・box問題0。主稿418頁・構造稿41頁・別稿91頁/118頁。修正ページを画像確認|

最小の新再現資産は `code/referee/` の四スクリプトと対応する `reports/referee_*.json`。
一般論・恒等計算・元producerによるserialize再現・標本だけの検査を区別している。
font/scaling/hyperrefの既存警告は残る。

## 再現ツールの訂正

依存監査の3536論理ファイル・16826記録hash・記録入力不一致/欠落0を再確認した。
未解決構文パス65個を自動的に「出力等」と断定していた記述を改め、
runtime read解析の完了を保証しないことを機械出力にも明記した。
exact入力materializerは欠落bindingを黙って飛ばさず例外にする。
既存の正確二bindingの準備成功と、欠落bindingの拒否を確認した。
これは分類の数学的反例又は必須入力65個の欠落という意味ではない。

## C1の現在有効な再現基準

**Libraryの実ZIPバイナリ再取得は、ユーザーの明示指示により必要条件から除外された。**
保存済みのローカル五volumeが保存時のSHA256と一致することを確認し、隔離復元で3579ファイルと相対リンク2個を検査した。
隔離treeで四PDF各二回組版を完了し、抽出テキストが同梱PDFと一致した。
この実行範囲は `reports/referee_reproduction_validation.json` にある。
Libraryへの保存の事実は `reports/referee_library_persistence.json` の保存receiptと参照の読み戻しで別途扱う。
Libraryバイナリの再ダウンロードを試行・要求しない。
査読合格の最終判定では、必要な主定理の証明・runtime入力の現実的な依存と今回の独立検証の射程を評価する。
既存の四PDF・五分割基礎bundle（source `a5f82342870b12583da0d0ca40bb5ba96c742b0e`）は保持する。
今回の査読補足は別に保存し、基礎bundle又は査読合格版への置換とはしない。
査読source commitは `e94cd15c4f1bad3d63d40929be32a2c1b55f5301`、main先端で照合済み。
その28ファイル補足を隔離treeへ適用し、四検算の全JSON一致を確認した。
改訂PDF二つと補足ZIPはLibraryへ追加し、Page sequence6で全十二参照（既存九＋追加三）を照合した。
主稿のPDF本文readbackは抽出テキスト1MiB上限で拒否されたが、uploadのdigest/accessは確認済み。
生成・隔離再実行・保存のreceiptは `reports/referee_delivery_validation.json` と
`reports/referee_library_persistence.json`。この保存receipt自体はsourceの次のcommitへ保存する。
無人の継続は約束しない。

資料Page：https://chatgpt.com/space/page_98c9971dc7c48191a2c4b7624fd99fa0
外部投稿・新たな外部公開、既存Libraryファイルの削除は行っていない。
