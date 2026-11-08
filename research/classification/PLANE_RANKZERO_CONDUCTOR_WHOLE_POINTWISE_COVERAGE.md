# rank零conductor plane：全first symbol・全位数の全点被覆

2026-11-08、PRIMARY E。統合担当の純合成証明。個々の全source証明と別著者receiptを使い、G次元上界から全点帰属を推測しない。新計算は無い。独立合成reviewは別receipt。

**命題。** 複素P4上のrank四Poisson tensorのactual Dが被約既約quinticで、rank零のconductor planeを持つなら、その点は既に採録された閉成分の有限和に入る。全位数、全first symbol、全higher係数を含む。

## 全入口

`PLANE_RANK_ZERO_WHOLE_SUPPORT_NONDOMINANCE.md` §§1–4 の全source分岐だけを再利用する。そのG結論をpointwise結論へ流用しない。`REDUCED_DEGENERACY_ZERO_SURFACE_SEPARATE_REVIEW.md` により被約D上のrank零曲面のconormalは可換で、全mのweight零35・一30・二9 sourceが有効である。traceless coHiggsの共通saturationはO(k), k=0又は1。k一/rankM一はglobal Vへ、k零もglobal Vへ戻る。rankM二以上、global V、traceless零の分岐は全ての零補集合を尽くす。

| 全source分岐 | 全点の結論 | 全源の証明入口 |
|---|---|---|
| O(1),rankM二以上 | 既知閉和又は空 | `PLANE_COHIGGS_RANK_AT_LEAST_TWO_WHOLE_POINTWISE_COVERAGE.md` と同POINTWISE_SEPARATE_REVIEW。全kernel/nonkernel、全Borel/nondegenerate、全τ、全normal位数 |
| global V、normalSS、T,V独立 | CartanはMPS、他は空又はessential三八成分 | `PLANE_RANK_ZERO_M2_COMPLETION_SEPARATE_REVIEW.md` のCartan全period/SNC、`PLANE_NORMAL_SS_TWO_PARAMETER_WHOLE_EXCLUSION.md`、`PLANE_NORMAL_SS_MIXED_THREE_PARAMETER_WHOLE_FACTOR.md` と `PLANE_NORMAL_SS_WHOLE_SEPARATE_REVIEW.md`、`PLANE_GLOBAL_LOWPARAMETER_POINTWISE_COVERAGE.md` とNIL_LOW21別review |
| global V、normalnil b非零 | 空 | `PLANE_NORMAL_NIL_AFFINE_WHOLE_ACTUAL_FACTOR.md` と `PLANE_NORMAL_NIL_LOW21_POINTWISE_SEPARATE_REVIEW.md`、全ρ及びhighm |
| global V、normalnil b零、T,V独立でleading determinant非零 | 空又はessential三八成分 | `PLANE_NORMAL_NIL_ZERO_B_MULTIPARAMETER_WHOLE_FACTOR.md` と `PLANE_NIL_B0_MULTIPARAMETER_WHOLE_FACTOR_SEPARATE_REVIEW.md`、全mixed/Cartanと全補集合。残低21は全点定理 |
| global V、T,V独立だがdet(E,T,V)零 | rank二以下、指定枝は空 | `PLANE_RANK_ZERO_COMMON_IMAGE_WHOLE_RANK_TWO.md` と `PLANE_COMMON_IMAGE_ROOT_WHOLE_REVIEW.md`。全四sourceと全P2 |
| global V、T=tV | essential三八成分又は空 | `PLANE_DEPENDENT_GLOBAL_VECTOR_POINTWISE_COVERAGE.md` とH二source独立review。全六sourceと全Ore八module |
| traceless零、T非零 | 空 | `PLANE_SCALAR_T_NONZERO_WHOLE_EXCLUSION.md` と `PLANE_SCALAR_T_NONZERO_ROOT_SEPARATE_REVIEW.md`。全15源/全壁 |
| traceless零、T零、K非零 | rank二以下、指定枝は空 | `PLANE_RANK_ZERO_ZERO_TRACE_BINARY_SEPARATE_REVIEW.md`。binary平方も含む |
| weight零symbol全零 | T36 | `PLANE_ZERO_SYMBOL_WHOLE_T36_ACTUAL_CURVE.md`、全39一軌道定理と別review。全higher/全cocycleの実逆 |

normalSS又はnormalnil b零でdet(E,T,V)が非零なら位数は二であり、有限first-symbol分類のgeneric二/三parameterと全low21で尽くす。独立でdet零なら全可換二場の分類によりcommon-imageへ戻る。dependentなら全m theoremへ戻る。この分岐で高位数をm二sourceのgeneric性から処理していない。

scalar、coHiggsとzero-symbolを低21へ暗黙に含めず、各全m定理へ分離した。全有限normalizer逆、六actual shearが作るinhomogeneous P2、FULL9、全代数的壁と元domainの零/非零は各別receiptに保存済みである。

## 閉和への結論

各表の全pointwise定理は元Poisson sourceを必要十分又はそれを含む必要sourceへ戻し、その**全実解**を既知閉成分又は空へ送る。従って有限表の和も既知閉和に入る。essential三源の全点八成分被覆とMPSのSNC定理を使い、像24以下だけの旧証明をpositive membershipに置換していない。

新成分加算零。rank二planeは別の全pointwise定理、非plane conductorは別の全pointwise定理で扱う。可約/非被約D及び全E/総上界の完成を本命題から結論しない。
