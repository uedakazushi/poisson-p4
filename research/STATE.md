# 研究状態 — 正確147成分

2026-11-08。方針 `RESEARCH_POLICY.md`、active queue `research/classification/PLAN.json`。

主分類の数学的義務は完了した。複素P4上のPoisson parameter spaceの被約支持（PGL5商前）は、**正確に147個の相異なる既約成分**の和である。一般rank四141、rank二6。台帳1158主張、未処理primarysource・未同定候補・未確認の分母補集合は零。
各成分の全空間成分性・相互相違・次元・一般rank・一般元のactual構成・退化因子幾何を台帳と同じ主原稿へ接続した。全PGL軌道分類、全schemeの非被約厚み、正確nil指数/H0/T2等secondaryは非主張。固定点の接次元からsmoothness/唯一局所枝を推測しない。

最終証明は `GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md`、最終manifestは `inputs/global_exact147_component_classification_final_proof_manifest_20261108.json`。全支持の無条件別作者最終査読は `FINAL_147_COMPONENT_SUPPORT_EXHAUSTION_COMPLETION_ADDENDUM.md`、最小5依存manifest `inputs/final147_component_support_exhaustion_completion_manifest_20261108.json` にPASSで凍結済み。
全caseの既142＋七閉候補被覆→unimodular P*の真T19実曲線で六→Ga9/plane三/constantDOUBLEの新五23全成分でknown147＋constantTRIPLE一22軌道→最後のTRIPLE actualPoisson族で既Ga9constant23へ。
wholeT37の新必要条件は「actualDのmultiplicity一component上でD自身smoothかつrank二」の点を次数≤5不変曲面で覆うという最弱範囲だけ。Ga9葉次数9、TRIPLE葉次数6、plane三とDOUBLEの非代数葉で親を排除した。重複因子のreducedsupportだけsmoothな点へ適用しない。
全rank四閉包を有限閉和へ入れた後に、rank二全支持 theoremのclosure項を代入し、循環は無い。全17rankone theoremはmarkedmultipleH限定。最後TRIPLE族のgenericfiberはsimpleHなので別のGa9完全十源へliteralactual正逆で入れる。

今回の新五成分はいずれも23次元・一般rank四：
- `R4_GA9_LINEAR_W_E6_CUSP_TANGENT_23`（E6平面四次と尖点接超平面、一般orbit23/商0）
- `R4_PLANE_TWO_LINES_CUBIC_NODE_FLEX_TANGENT_THROUGH_FLEX_23`（節点三次の変曲接線と、その変曲点を通る通常接線、orbit22/商1）
- `R4_PLANE_TWO_LINES_CUBIC_CUSP_FLEX_SECANT_THROUGH_FLEX_23`（尖点三次の変曲接線と変曲点割線、orbit22/商1）
- `R4_PLANE_TWO_LINES_CUBIC_CUSP_CUSP_TANGENT_CUSP_SECANT_23`（尖点接線と尖点割線、orbit22/商1）
- `R4_CONSTANT_MODULAR_DOUBLE_LINE_CUBIC_23`（二重Hとcuspidal cubic、orbit23/商0）
全source境界を閉像として保持する。特殊磁気点と他成分の交叉無し、generic scheme厚み、局所smoothnessは追加主張しない。最後のTRIPLEは既成分の境界なので加算零。

独立exact再生：P*→T19同曲線は全150/350/126/両逆/固定符号PASS（2.827秒）、TRIPLE→Ga9同族は元積Poissonから全150/350/126/両逆/actualframe/Phi_uPASS（0.315秒）。別作者は元whole-source仮定と全c1十parameter/六括弧を純査読PASS。証明読解、算術再生、plancheck/gateを区別し、後二者だけで数学を証明したとはしない。
旧失敗・限定FALSE・quadratic余核零・単独ansatz陰性はimmutableのまま保存。旧149/148上界は弱い正しい過去証明であり、現行未完キューではない。全古いsource/code/report/job/manifestは上書きしない。

この再開区間の開始mainは `2b7c071f168d0709bb28e1cb436ec137d91b4200`（旧142/1125/349頁）。Cloud開始は `4c3f623fb698f9adc6077015b25e2325b11072af`、数学baselineは `a00a541b835b2bff18bafcba545ecbd4bae7175b`。
現行日本語主稿は**147成分/1158主張・368頁**へ統合完了。最終二回LuaLaTeX・参照/引用/box0・新規警告0、editor23頁及びroot8頁の可視独立確認PASS。最小再現入口は `REPRODUCE_EXACT147.md`。
TeX SHA `6062102738663ba55d6797b16049d05d2513e80ab19459e894a3a056816b6b1a`、PDF SHA `784539e12a0ce9f4ceb548f4c783ac9c2b9b93618fb4f39b544a9f421fc35339`。root最終組版/受理receiptは `reports/CLASSIFICATION_MAIN_MANUSCRIPT_COMPONENT147_CLAIM1158_BUILD_20261108.json`。
数学・主稿・再現資料・独立確認のprimaryは全て完了。通常main commit/nonforcepushの実状態は `git rev-parse HEAD` と `git ls-remote origin refs/heads/main` の一致で確認する。本文で自分自身の将来commit SHAを先取りしない。

既Cloud .venv/Macaulay2 1.24.11/SymPy/Fraction/LuaLaTeXを再利用。引継ぎ元1.26.05との必要QQ互換・軽量動作確認は既完了、setup全体の再実行不要。Fedorchuk一次PDFを直接再読、MathNet shell再接続HTTP403はreceiptに区別したが、必要formalversal入力は直接読めてprimary阻害無し。
復旧は `research/CLOUD_RECOVERY_20261107.md`。旧二legacy optional runtime/gateや一historical未commit review lossは現採録の必要証明に使わない。現primaryの接続・権限・資源阻害無し。
rootのみ台帳/Git、editorのみ主稿/frontier。旧local移行一時停止は既に再開済みの研究を止める指示ではない。研究完了後にsecondaryを無断再開しない。将来は具体的な証明欠陥又は新しいユーザー指示を入口にする。
