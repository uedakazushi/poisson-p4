# 現在地と最終到達条件

2026-11-09 / v72-normalization-residue-integration

**全工程は未完。中間保存後もFINAL_TARGETの未達条件へ連続して進む。**
最終仕様ad2800acと数学9830a71fはa205a79dで統合しmainへの保存を確認した。
現行稿・台帳は147成分（rank4:141、rank2:6）を主張するが、今回その全証明を独立認証したとはしない。

## 達成証拠と未達条件

| 条件 | 現在の証拠と、なお必要な作業 |
|---|---|
| M1 全点被覆 | 主稿prop:irreducible-nonnormal-all-pointwiseの全次数/rank表を保持。thm:normalization-residue-source、prop:conductor-pole-residue、thm:curve-residue-sourceが全降下・留数とdegree5/6を共通化。R2全支持をsec:rank-two-reconstructionへ完全統合。安定な五切断の全証明と因子付き全枝を実sourceへ戻した |
| M2 極大性・相違 | 旧有限比較を保持。cor:reduced-divisor-component-imagesにより一般因子被約な別成分は同じ閉因子像を持てない。maximal_image_synthesis.texで全周囲・真正切片・有限閉和の三機構と全必要比較を合成。計算証明書の依存閉鎖はR1で続行 |
| M3 具体的記述 | 台帳と案内表あり。component_family_catalogue.texに全147 ID、明示source/全源参照・開条件・次元・階数・実因子/foliation・成分性参照を統合。51有理演算列を正確評価して全射影Jacobiを確認（一般証明の代用ではない） |
| M4 依存閉鎖 | Kaledin math/0310173v1 Theorem1.5/Corollary1.6を一次全文で確認。一般源の証明・37記号検算を追加。最終主定理の必要証明書依存閉包は未解消 |
| C1/C2 概念的統合 | A/B/planeに加え正規化・導手pole bound・全B/A線形降下・scroll全留数kernelを統合。prop:partial-conductor-familyが二つの実T40曲線と023族の重複を置換。全10係数determinant境界はprop:plane-fixed-divisor-pencil。R2の全証明とtoric/Aff32の閉包項を接続済み。全極大閉像の最終依存比較と本文の全ID対応は未完 |
| P1 公表原稿 | 現行主稿417頁・構造稿41頁の中間版を二回組版、追加と置換箇所をPDF画像で確認。原稿全体の理解順の再構成、主要一次引用・全一般元表の最終確認は未完 |
| R1 再現一式 | reports/normalization_residue_checks.jsonに新37項目、既存共通29項目を再現。一般証明の代用ではない。必須旧証明書のbundle依存閉包は未収録で、最終条件未達 |
| S1 保存 | 仕様・数学の統合commit a205a79dをmainで確認。R2統合5ec71f37もmainで確認。最終bundle/Library保存は未達。既存Libraryファイルは削除していない |

検証の正確な射程とPDF hashはreports/normalization_residue_validation.json。
被約fixed-D fiberはthm:reduced-fibre-linearにより線形開で幾何学的既約だが、非被約h5には反例がある。
導手次数だけでは降りないcuspのda/aと、二位極da/a2が許される実例、node-cusp相対導手を保持する。
全ID・成分数をこの新定理だけで再認証したとはしない。優先権は未確認。

## 次の具体操作

R2の完全証明、M3全IDカタログ、M2の極大閉像の共通論理を本文へ統合した。
reports/component_catalogue_integration_validation.jsonに本文417頁、二回組版、参照/box検査、実画像確認、51正確評価の射程を固定。
最終依存の223 primary claimから直接754ファイルを追い、typed closureは3362ファイル約661MB。
この機械走査には不要なprovenance/outputのリンクもあり、全件が証明入力との主張はしない。
元コードhashの不一致43種のうち、8種をGit、34種を保存証明書の埋込sourceから正確回収。
残る一つは更新済み数学ノートの旧hashで、旧未採択proof、出力ログ、歴史tree保存記録を正確な数学依存と区別する。
必須work入力rank3-mixedはreportsの同一bytesをhashで回収した。
次はこの依存解決を再現スクリプトと現行manifestへ固定し、必須code/input/certificateの欠落を解消する。
bundleが旧ZIP/未同梱入力依存を解消するまでR1完了とはしない。その後mainとLibraryへ保存確認する。
数学の完成又は全147の独立再検証を、この中間作業だけから報告しない。
外的阻害を終了理由として認定しておらず、未達条件へ続行する。
