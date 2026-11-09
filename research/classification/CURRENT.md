# 現在地 — 公表候補版の検証・保存完了

2026-11-09 / v74-publication-candidate

**数学的査読の判定：指定した高リスク義務を充足し、公表候補版として可。**
複素P4の射影Poisson空間（PGL5商前）の被約成分は、現行稿で147成分
（一般rank四141、rank二6）と分類する。新たな数学的反例・未解決の必須補題は残っていない。
再利用した有限認識と独立再計算の範囲は `research/publication/REFEREE_REPORT.md` に明記した。
外部査読の採択、新規性・優先権の確定、全PGL軌道・全scheme厚みの分類は主張しない。

## 査読の根拠

|対象|確認した義務|
|---|---|
|A1 固定被約五次|消滅定理からlift/Jacobiを独立導出。350三次二形式係数の全核を三例で別計算し、全基底対の互換性・Pf=qFを確認。q恒等零・rank jumpと因子像一致/包含を区別。相対系の底の被約性を明記|
|A2 rank二|飽和反射層、五切断、非LCM二平面拡大、spinorからRat113/122へ。非安定/因子付きの全枝とCM/LPT一次資料、六成分とrank四境界の区別。八係数Jacobian恒等式を独立計算|
|A3 最後五像|全五族・全磁気六係数、proper flags、whole T37の標識付きHilbert極限、全prior142非包含。葉degree9/6と三磁気源の留数を別導出。零磁気の個別外性へ拡張しない|
|B1 全51一般族|全パラメータJacobiを別実装で確認。T39/T40は保存DAGの全16入力で直接検算。残49DAGの214525節点と全入出力を再生成。一般正規化・稠密性・MPS40の正確な射程を読解|
|B2 全支持|全因子/Jordan/translation余核→有限閉和→rank二閉包代入→極大像の非循環性。元有限認識・成分性証明書を再利用し、223主張の全再証明とはしない|
|C2 原稿|既知MPS40と147総数、安定rank二の内部復元、独立S61補論を導入で区別。主稿418頁・構造稿41頁・別稿91頁/118頁。既改訂箇所の画像確認と四稿二回組版済み|

## C1と実入力の閉鎖

Library ZIPのバイナリ再取得は指示により行わない。
既存ローカル五巻は保存SHAと一致し、3579ファイルと二相対リンクを隔離復元、
四稿二回組版と同梱PDFとの全抽出テキスト一致を確認した。
その隔離treeに査読補足を適用し、四検算の全JSON一致も確認済み。

さらに、65構文パス・206組を出典SHAとASTで読解した。
実入力8、書き出し先51、native生成出力4、事前生成入力1、archive locator1で未分類0。
旧bundleに欠けていたv21の二runnerを必須依存へ追加した。
八支持runnerの不要な文献PDF依存を任意化し、PDF無しの隔離実行で元の算術・被覆結果と一致した。
元runnerの正確な一版を追加保存し、凍結入力は97件。
weight6事前生成入力も隔離再生成して時間以外の全fieldが一致した。
追加したstdlib checkerで八模型の全Jacobi/Pfaffian等249項目を確認した。

最終の記録入力監査は3539論理ファイル・16840記録hash、欠落/不一致0、未分類構文パス0。
任意の生成コードに対する完全な動的read解析ではなく、報告された依存の有限読解である。
RUNTIMEという名前だけの免除は撤去し、二つの既知失敗metadataだけを正確な出典/hashと置換証明書で除外した。
根拠は `static_path_review.json` と `reports/referee_runtime_input_validation.json`。

## 配布の実行状態

全source・必要入力・査読報告を統合したv74全bundleを生成・検証し、五巻をLibraryへ保存した。
配布元はmainへ保存済みの `fdec6f8bd7eb01072093c05008f11932ed628046`、
bundle IDは `b3818c83dc41a73ae53cac9248c2e581a59c2be203b002934648984c1a4bd6e9`。
五巻だけから3599ファイル・二相対リンクを新しい隔離treeへ復元し、全chunk/file SHA256と論理参照を確認した。
最終復元treeでも206判定、八模型249項目、PDF無しの八支持runnerを再確認した。
四稿二回組版の全TeX入力hashと最終bundleの入力hashが一致する。
旧五巻と補足のoverlayは不要。

Libraryの五巻は保存先SHA256・サイズが生成物と一致した。
現行四PDFは同一hashの既存保存物を再利用し、現行九ファイルの参照・MIME・サイズを照合した。
資料Pageはsequence 7、全17参照で、更新前の12参照をすべて保持している。
生成・復元・限定再実行の証拠は `reports/referee_runtime_input_validation.json`、
保存・Page読戻しのreceiptは `reports/publication_candidate_delivery.json`。
配布元commitと、保存完了を記録する後続commitを分け、receiptの自己参照を避ける。
現在の目標に必須の未完作業はない。新しい証拠又は具体的な誤りがあるときだけ関連する依存を再開する。
Library ZIPのバイナリ再取得、既存Libraryファイルの削除、外部投稿・公開は行っていない。
資料Page：https://chatgpt.com/space/page_98c9971dc7c48191a2c4b7624fd99fa0
