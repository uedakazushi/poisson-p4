# 現在地と最終到達条件

2026-11-09 / v72-normalization-residue-integration

**全工程は未完。中間保存後もFINAL_TARGETの未達条件へ連続して進む。**
最終仕様と継続命令の更新ad2800acを、数学チェックポイント9830a71fと統合している。
現行稿・台帳は147成分（rank4:141、rank2:6）を主張するが、今回その全証明を独立認証したとはしない。

## 達成証拠と未達条件

| 条件 | 現在の証拠と、なお必要な作業 |
|---|---|
| M1 全点被覆 | 主稿prop:irreducible-nonnormal-all-pointwiseの全次数/rank表を保持。thm:normalization-residue-source、prop:conductor-pole-residue、thm:curve-residue-sourceが全降下・留数とdegree5/6を共通化。R2全支持の共通接続をこれから確認する |
| M2 極大性・相違 | 旧有限比較を保持。cor:reduced-divisor-component-imagesにより一般因子被約な別成分は同じ閉因子像を持てない。全147の必要比較・特殊底・非被約入口の依存確認は未完 |
| M3 具体的記述 | 台帳と案内表あり。全IDの一般元・条件・次元・階数・幾何・証明と最終本文の一対一整合は未確認 |
| M4 依存閉鎖 | Kaledin math/0310173v1 Theorem1.5/Corollary1.6を一次全文で確認。一般源の証明・37記号検算を追加。最終主定理の必要証明書依存閉包は未解消 |
| C1/C2 概念的統合 | A/B/planeに加え正規化・導手pole bound・全B/A線形降下・scroll全留数kernelを統合。prop:partial-conductor-familyが二つの実T40曲線と023族の重複を置換。全10係数determinant境界はprop:plane-fixed-divisor-pencil。R2と全極大閉像の最終接続は未完 |
| P1 公表原稿 | 現行主稿385頁・構造稿31頁の中間版を二回組版、追加と置換箇所をPDF画像で確認。原稿全体の理解順の再構成、主要一次引用・全一般元表の最終確認は未完 |
| R1 再現一式 | reports/normalization_residue_checks.jsonに新37項目、既存共通29項目を再現。一般証明の代用ではない。必須旧証明書のbundle依存閉包は未収録で、最終条件未達 |
| S1 保存 | 数学のローカル中間commit9830a71f、最終仕様remote ad2800acを非forceで統合中。最終bundle/Library保存は未達。既存Libraryファイルは削除していない |

検証の正確な射程とPDF hashはreports/normalization_residue_validation.json。
被約fixed-D fiberはthm:reduced-fibre-linearにより線形開で幾何学的既約だが、非被約h5には反例がある。
導手次数だけでは降りないcuspのda/aと、二位極da/a2が許される実例、node-cusp相対導手を保持する。
全ID・成分数をこの新定理だけで再認証したとはしない。優先権は未確認。

## 次の具体操作

RANK_TWO_SUPPORT.md、KNOWN_COMPONENT_AUDIT.mdとstable_rank2原証明を読み、
prop:ranktwoの飽和foliation・零因子から全支持へ接続する。
六つの純rank2成分と、toric外積/Aff32等のrank4閉包項を区別し、
全rank4被覆後に閉包項を戻す論理順を主稿に統合する。
その後M2/M3、主稿再構成、必須証明書の全依存・bundle・main/Libraryへ直ちに進む。
現在、数学・実行上の外的阻害を終了理由として認定していない。
