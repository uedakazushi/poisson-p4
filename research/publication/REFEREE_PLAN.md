# Poisson P4 完全分類：公表前の独立査読の現在計画

2026-11-09 / v74-publication-candidate

以下は本査読で適用した検査基準である。結論は REFEREE_REPORT.md、配布工程の状態は
classification/CURRENT.md に記す。全bundle検証と保存も完了した。
以下は適用済みの基準であり、既済みの導出を繰り返すための未完作業一覧ではない。

## 0. 判定：次の仕事は新分類ではなく、最も脆い証明への反証試験

複素 P^4 の射影Poissonスキームの被約既約成分（PGL5商前）について、
現行主稿は147成分（一般rank4:141、一般rank2:6）を主張する。
主稿418頁、構造稿41頁、正規因子稿91頁、一次還元稿118頁、147族の全カタログ、
五分割の証明・再現bundle、既存の最終条件達成報告がある。
Git mainの完成宣言を、査読の証拠そのものとしては使わない。

**研究方針の結論：今から新しい成分数の探索、抽象化の追加、テスト項目増加を優先するのは得策でない。**
優先すべきは、分類全体への依存度が高く、かつ仮定の取り違えや
「標本検算→全パラメータ」という飛躍が起こりやすい箇所の独立確認である。
以下をその証明依存順で進める。反例があれば直ちに依存する結論を見直す。
確証できないときは未確認を具体的に残す。原稿・図式・成分数を守るために仮定を増やさない。

既存の FINAL_TARGET.md の完成宣言は過去の構築段階の成果であり、
本査読は別の「投稿準備品質」フェーズである。
完成宣言を真とも偽とも初めから決めず、検査可能な証拠で結論する。

## 1. 優先順位A：再利用される主張に対する反証志向の査読

**A1：被約実退化因子上の全fiberが線形開であるという定理。**
対象 papers/normalization_residue_reconstruction.tex の
lem:quintic-ambient-lift、thm:reduced-fibre-linear、
cor:reduced-fibre-boundary、cor:reduced-divisor-component-images。
固定した被約五次式Fに対する V_F=ker(F dα−dF∧α) が
互換Poisson双ベクトルからなり、全rank4のfiberがP(V_F)\V(q_F)となることを
原稿の論証を前提にせず導出してみる。非正規・交叉・singularな因子、
q_F恒等零、rank jump、generic fibreの幾何学的既約性、
divisor imageの「一致禁止」と「包含は含意しない」の区別を一つずつ検査する。
Jacobiが零であることとPfaffianがcFであることを混同しない。
小さな反例候補を独立コードで探索し、見つからなければ一般証明を記す。
反例探索の成功・不成功だけを定理の証明と呼ばない。

**A2：一般rank2の安定分岐の全Jacobian復元。**
対象 papers/rank_two_support_reconstruction.tex の
prop:ranktwo、thm:stable-rank-two-actual-rational、五切断・延長・spinor等の補題、
prop:rank2-support、六成分のambient componenthood。
「foliationが代数的可積分である」ことから
Rat(1,1,3)又はRat(1,2,2)という具体的全源へ進む全論理を追う。
飽和・反射性・det twist・零因子なし・stable/nonstableの分岐、
Aff32とtoricのrank4閉包、閉包項をrank4被覆の後に戻す非循環性を確認する。
Corrêa–Muniz（arXiv:2207.12880）とLoray–Pereira–Touzet（arXiv:1112.3871）の
一次本文の仮定・結論を直接比較し、外部定理が保証する以上を引用しない。

**A3：最後の五つの23次元成分と既存142成分の非包含。**
対象 papers/maximal_image_synthesis.tex、
papers/component_family_catalogue.tex の source:last-five-explicit、
研究証明 GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md と
その依存する全proper-flag・T37葉曲面・磁気不変量。
五族の全パラメータJacobi・actual因子を独立実装で再構成する。
それらと既存142成分の非包含を、単なる次元、一般退化因子の違い、
ソース一般点の帰属だけで処理していないことを確認する。
zero/unit・分母・異なる磁気frame・特異な旗を全域で保持し、
同次元の他候補への包含を漏らしていないか、明示的な退化族を試して反証を探す。
有限閉被覆の抽象補題だけでは成分性は得られない。実際の全source
及び非包含が独立の入力であることを示す。

A1–A3の各項目について、「正確な主張／使用した仮定／原証明箇所／
自力で導いた検査／最初に疑った失敗点／結果／依存する主定理」を短く記す。
未確認なら具体的な未解決補題を残す。
一つの入口に局所的な証明穴があれば、単なる問題提起で止まらず、
正しい補題・修正・反例に応じて主稿と全依存を直す。

## 2. 優先順位B：列挙した族が真に「一般族」であるか

**B1：T01–T40と追加十一族の全パラメータ源。**
対象 papers/component_family_catalogue.tex、code/data/component_manuscript_catalogue.json、
code/data/known51_rank4_s61_coordinate_catalogue_20261103_v62.json、
code/build_component_catalogue.py、code/s61_model.py と現行51演算列。
現行の「正確な標本評価PASS」が保証するのは座標規約と一標本のJacobiであって、
任意パラメータのJacobi、mapの正則定義域、開集合の非空性、
既約性、像の稠密性ではない。
MPSの40成分の原定理が主張する射程を一次本文で確認し、
smoothing diagram・period・full symplectic frame・S61正規化の
どの一般定理が全族の存在と稠密性を保証するかを、具体的な定理番号・式で追う。
証明が見つからなければ、対象となる実族について必要な全パラメータの式を
独立に証明または計算証明書で補強する。全51族で巨大な記号展開を
機械的に実行することを最初の方策にしない。
非零分母、forced-zero、√3の選択、特殊パラメータ・rank落ちを確認する。

**B2：全点被覆と147個の極大閉像への接続。**
A/B/C/rank2の全源→既存閉像→非包含という証明の
「最小の必須依存節点」を特定して再読する。
reports/required_proof_dependency_closure.jsonの223 primary claimや
凍結入力（現行97件）を「数学的査読済み」と読み替えない。
特に最後の五成分の存在、141 rank4の全support、
安定rank2の6成分、一般rank4閉包をrank2へ使う順序を検査する。
全147×147の交叉表や全係数証明書の再実行は必要ない。
欠落を見つけた場合は、その枝の元全sourceへ戻り、
条件付き結果・未証明の推論・完全な証明を区別する。

## 3. 優先順位C：クリーンな再現と公表原稿の読み直し

**C1：再現性の検査。**
現在のリポジトリ内の再現資産、または既にローカルにある五つのZIP volumeを使い、
UNPACK.pyの指示通り隔離ディレクトリへ復元し、
manifest、全SHA256、相対参照、必要な93/96等の番号・本数を
本当の入力ファイルから再確認する（本数を推測で一致させない）。
**LibraryからZIPのバイナリを取得・ダウンロードすることは要求せず、試行もしない。**
Libraryについては既にある保存receiptとファイル参照・メタデータを確認すればよい。
実データの復元性はリポジトリの実ファイルと既存のローカル五巻から検証し、
その証拠の射程を明記する。
四原稿のLuaLaTeXを再組版し、引用・参照・数式番号・欠けた図表を確認する。
独立再実行を選ぶ際は、根拠となる証明書とproducerを初めから同じものを
無批判に呼ばず、独立の簡単な導出または別実装と照合する。
必要な有限証明を支える資産だけ残し、歴史的な旧版・旧ログ・キャッシュは除く。
ただし原証明の再現に必須な凍結版は旧形式でも保持する。
ソース/文献PDFを不必要に大量収録しない。
すべての再実行完了を約束するのではなく、主定理に必要な証明依存の閉鎖を確認する。

**C2：投稿前の編集上の査読。**
主稿418頁・構造稿41頁・別稿91頁/118頁の全依存関係を読者の順に確認し、
本論文の中心定理とMPS40/Pym/Corrêa–Muniz等の既知結果との違いをIntroductionに具体化する。
重要な関連文献は一次資料で再照合し、書誌と適用範囲を推測しない。
主稿中の証明を単に別冊に移して短縮したと主張せず、必要な証明は
主稿または同時配布の自己完結した別論文に保持する。
一冊が明らかに読みにくければ論理的に完結する複数論文へ分割してよいが、
ページ数削減だけで分割しない。貢献・仮定・既知部分・新規性未確認を明示する。
原稿の内容と本文参照を変えず、README/RESEARCH_POLICY/NEXT_TASKS/STATEの
完成／未完文言が矛盾している点を現行事実へ同期する。
論文の外部公開・投稿は指示がない限り実施しない。

## 4. 実装手順と評価基準

最初に FINAL_TARGET.md、CURRENT.md、全四つの原稿の目次と主定理、
必要な証明・主要検証を読む。まだ元の主定理の証明を読まないうちに
大量の新しいテストや新しい研究計画を書かない。
まずA1を独立に導出し、反例を探す。次にA2、A3、B1、B2、C1、C2と進む。
証明依存が変わる場合のみ優先順序を変更する。

新しい作業の成果は、現在有効な **独立査読の結論** と
実際に必要だった定理・証明・原稿修正とする。
毎回「PASS」「形式的な別証明」「監査表」の数を増やすことを成果としない。
証明が十分な箇所は簡潔に理由と所在を記し、次の高リスク節点に移る。
失敗した方法・未完成な補題は必要な障害だけノートに残し、全作業実況は残さない。
対象を無限に広げず、主要な依存をすべて通過したら公表候補版を確定する。

**最終成果物**：
- research/publication/REFEREE_REPORT.md：高リスク項目の合否・根拠・未解決事項・推奨の簡潔な現行版（作業ログではない）。
- 必要な箇所を訂正した現行LuaLaTeX論文・共有TeX・研究ノートと、実行した独立検算の最小再現資産。
- 公表用の独立した証明依存関係とIntroductionの先行研究比較、再生成可能な現行bundle。
- mainのnonforce pushと先端確認、Libraryの追加承認不要の保存。既存Libraryファイルを削除しない。

査読の**合格**には、A1–A3およびB1/B2の証明上重要な義務が
原稿中の通常の証明・正確な一次文献・検証可能な証明書によって閉じていること、
C1/C2の再現・読者確認が済んでいることが必要。
一部が未証明なら合格を宣言せず、修正を進める。
数学的誤りを発見したら、その依存先・成分数・主定理・Introductionの範囲まで見直す。
実行限界など外的理由の中断では、可能な作業を先に行い、未達箇所・再現可能な阻害・最後のcommitを
CURRENT.mdへ記して未完と報告する。無人で将来続くとは約束しない。

## 5. 文献の一次確認先と既存資料

- B. Pym, Constructions and classifications of projective Poisson varieties,
  Lett. Math. Phys. 108 (2018), 573–632, DOI 10.1007/s11005-017-0984-5;
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5818673/ 。
  2018年のfull classification未解決という記述は2026年の既知の完全分類不在の証明ではない。
- M. Matviichuk, B. Pym, T. Schedler, A local Torelli theorem for log symplectic manifolds,
  arXiv:2010.08692; https://arxiv.org/abs/2010.08692 。
  Theorem 1.4、Theorem 6.16、Corollary 6.23、Example 6.24とSNC-meeting40の射程を確認する。
- M. Corrêa, A. Muniz, Holomorphic foliations of degree two and arbitrary dimension,
  arXiv:2207.12880; https://arxiv.org/abs/2207.12880 。
  結論を「特定Jacobian構成の全源」へ無断で強めない。
- Kaledinの正規化定理、Gualtieri–PymのPoisson module、
  Loray–Pereira–Touzetのfoliation定理は、現行主稿の引用と原論文で正確に確認する。
- 既存の主要証明：papers/normalization_residue_reconstruction.tex、
  papers/rank_two_support_reconstruction.tex、papers/maximal_image_synthesis.tex、
  papers/component_family_catalogue.tex、
  research/classification/GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md。
- 現行検証：reports/final_manuscript_validation.json、
  reports/component_catalogue_checks.json、reports/dependency_reconstruction_checks.json、
  reports/final_bundle_validation.json。これらは数学的な独立認証の代用ではない。

現在の分類稿と有限証明書を無断で破棄したり、旧版をアーカイブとして増殖させたりしない。
