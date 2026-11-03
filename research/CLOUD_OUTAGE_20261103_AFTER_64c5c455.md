# Cloud 接続障害と未 push 成果の再開記録 — 2026-11-03

この記録は研究の完了宣言ではない。Cloud 実行環境が unavailable / offline となり、実行・ローカル保存・通常の git push ができなくなったための障害引継ぎである。GitHub 接続だけは利用できたため、この記録のみ GitHub Contents API による追加コミットで保存した。研究コードを通常 push できたという意味ではない。

## 確実に GitHub 保存済みの状態

- 全 Cloud 研究の開始 SHA: `4c3f623fb698f9adc6077015b25e2325b11072af`。
- 最後の通常 commit / push: **`64c5c455ac642cbeacbf59fc63fd638e80102712`** (checkpoint J)。push 成功と main の一致を切断前に確認した。
- 直前の I: `b5545eeadc928bb5a0cb11cb2d9906f7113eedc2`。H: `6d82c007d497c5e6edbfdd006274d3133bc405fd`。
- J の CLAIMS は364件。新しい全 Poisson 成分の加算なし。下界57、上界未確定。
- J audit: `reports/V63_CHECKPOINT_J_AUDIT_20261103.json`。1319件の再帰的な現行入力 hash、全 claim 参照、v6 archive の49 member を確認済み。
- J integrity: `reports/20261103T105143.729205Z-integrity/verification_summary.json` PASS。
- `research/STATE.md`、`research/NEXT_TASKS.md`、`research/CLOUD_HANDOFF_20261103_v63.md` と v63 source index は J まで同期済み。
- v62 の凍結主原稿・旧 source index は変更していない。

J の主な確定結果は、全被約 R∩N=B18、全 R の Sing=Nonnormal=NonCM=B17∪C15、C15の全 pinched 完成環、generic B17 二枝環、実 CA/CB 交叉の minimal 一般長 N1/R2、実 Rat113 A20 prime と全 R37 の明示 a-torsion、全真四次初期式の被覆、実 blowup による二つの15平面の一般方向の除去である。全35模型の次元上界14/12や全枝の尽くしとは区別する。全 Qa の T02/T03 帰属も J に保存済みで、J 時点の独立確認済み点上成分下界は T01/T02/T03/CA/CB の5つ。

## 障害の観測と実際の制約

対象の既存実行環境:
`ccarenv_b64_Y2NhcmVudl9mNmRhOGJlZjMwMzg4MTkxODhkNDhhNzdkODQxOWQ2ZQ`。

ユーザーの environment_context は unavailable。exec は「No usable execution environment is available」で失敗した。その前には exec-server transport disconnected / recovery timed out25s も観測された。Cloud Environment の status は desired_phase=running、observed_phase=running、connectivity=offline、failure=null。したがって停止原因そのものは確定していない。権限拒否、ネットワーク allowlist 不足、セットアップ破損だとは結論しない。利用可能なツールに、この既存 runtime を restart / reconnect する操作はなかった。

GitHub connector で main の STATE.md を読めたので、この障害記録だけは remote に保存できた。Cloud のローカルファイルは切断後には読めていない。以下は切断直前に作成完了を確認したものの一覧であり、復旧後に同じ workspace 上で保持と hash を確認する必要がある。

## 未 push だが、producer と独立監査まで完了した成果

次の成果は **J の GitHub 保存に含まれていない**。ローカルで凍結済みと確認したもの。

1. **generic B16 の全完成像。**
   - `code/ambient_R19_B16_completed_image_20261103_v63.py`
   - `research/notes/AMBIENT_R19_B16_COMPLETED_IMAGE_20261103_v63.md`
   - `reports/AMBIENT_R19_B16_COMPLETED_IMAGE_20261103_v63/report.json`
   - 独立 review: `code/ambient_R_B16_review_20261103_v63.py`、同 stem の REVIEW note / JSON、および REVIEW directory の `review.m2`, `review.log`。
   - generic open B16 orbit で Rhat=C[[u16,t,r,c,p²,pr,pt]]、normalizer は r²−25t²−2cp=0 の quadric。depth18、edim22、multiplicity3。全 B16 / 実交叉厚みの尽くしではない。

2. **extended toric arcs の四成分。**
   - `code/q_mps_extended_toric_arcs_20261103_v63.py`
   - `code/data/q_mps_extended_toric_arc_witnesses_20261103_v63.json`
   - note / reports: `Q_MPS_EXTENDED_TORIC_ARCS_20261103_v63`
   - 独立 review: `code/q_extended_arc_review_20261103_v63.py` と同 stem の REVIEW note / JSON。
   - 全有限 Qa が T02/T03/T16/T17 に入ることを確認。J の結果と合わせ、独立監査済み下界は少なくとも7既存成分。全150係数、divergence/Jacobi、period、密度と closedness を確認済み。

3. **非平方の有限次数 unit lift。**
   - `code/s61_binary_nonsquare_bounded_unit_lift_20261103_v62.py` と二つの generic native utility。
   - note / reports: `S61_BINARY_NONSQUARE_BOUNDED_UNIT_LIFT_20261103_v62`
   - root の独立 review: `code/s61_bounded_unit_lift_review_20261103_v62.py`, `code/nmod_selected_rank_review_20261103_v62.cpp`。
   - `reports/S61_BINARY_NONSQUARE_BOUNDED_UNIT_LIFT_20261103_v62_REVIEW/` に report.json, rank.log, selected-matrix.in (約4MB) を保存済み。
   - 154項の次数3 scalar aJ witness を別の Fraction 実装で検査。第22列の witness 係数は非零有理定数であり、置換行列式1ではなく定数単元である。
   - 強い isolated target の次数≤5係数行列11302×2577を独立再構築し、2577正方小行列の full rank を別素数32003の独自 FLINT rank で確認。QQでこの有界解がないことのみを証明。
   - 13×45縮約の全残余 SHA を fresh producer と比較済み。全 core rank、自由 coker、四次65被覆、周囲除去は未証明。

## producer 凍結済みだが独立再検証が未完了の成果

以下は producer が全証拠を作成・確認して凍結したが、独立 coefficient/minor replay の PASS をまだ保存していない。未解決を完了扱いしない。

- `Q_MPS_IMPLICIT_TORIC_ARCS_20261103_v63`:
  `code/q_mps_implicit_toric_arcs_20261103_v63.py`、対応 witness JSON / note / reports。
  17個のQQ source witnessから T04–T15,T18–T21,T28 の全 Qa 帰属を主張する packet。producer PASS6.68s。証明橋は読了・概念監査済みだが、独立した全17係数/minor再計算は未完了。
- `Q_MPS_HENSEL_TORIC_ARCS_20261103_v63`:
  `code/q_mps_hensel_toric_arcs_20261103_v63.py`、対応 witness JSON / note / reports。
  16個の正の p-smooth source witnessから T22–T27,T29–T38 の全 Qa 帰属を主張する packet。producer PASS1.43s。positive Hensel/etale existence の橋を使い、mod p 空性からの持上げではない。独立再計算は未完了。
- 上二 packet を含む producer の結論は、全 Qa∈T01–T38,CA,CB、点上で少なくとも40既存成分である。**この40を独立再検証済みと記録しない。** 成分総数の下界57自体は変わらない。
- reviewer が開始した fresh outputs:
  `work/v62-review/Q-MPS-implicit-toric-arcs-replay` と
  `work/v62-review/Q-MPS-hensel-toric-arcs-replay`。
  exec sessions42825/25938は切断後の完了を観測できていない。復旧後まず実出力を読む。
- `Q_CA_CB_N_SOURCE_SMITH_20261103_v63`:
  `code/q_CA_CB_N_source_smith_20261103_v63.py` と note / reports。全QQ[a] source image Smith。CA rank28、N27は全a、CB28はa非零、a0で27。target tangent rankと混同しない。
- `Q_CA_CB_N_ALL_NONZERO_COMPLETED_MODELS_20261103_v63`:
  `code/q_CA_CB_N_all_nonzero_completed_models_20261103_v63.py` と note / reports。
  全物理 a≠0 で CA/CB=A8、N=A7 の旧完成型が維持されることを主張。固定βの Qa は自由な k0=9a/2 だけを動かし、missing odd idealにk0はない。物理 a0 の完成型は主張しない。
- `Q_CA_CB_ACTUAL_INTERSECTION_TRANSLATION_20261103_v63`:
  `code/q_CA_CB_actual_intersection_translation_20261103_v63.py` と note / reports。
  共通 source shift k0+=λ,k2−=pλ,k5+=p²λ/3 が target の u38 だけを変える。旧218式と45 colon式が該当v4に依存しないことを確認し、旧 actualD の平方零nilradical、Ass2、full depth20等を全物理a≠0へtransportする packet。独立監査未着手。
- 上の CA/CB/N 三 packet の reviewer は証明を一部読んだ段階で切断された。新しい REVIEW code/note/report は未作成、PASS なし。
- `Q_GLOBAL_WEIGHT6_TWO_AXES_20261103_v63`:
  `code/q_global_weight6_two_axes_20261103_v63.py`、note / reports (certificate.json.gz 約624KB)。
  a単元の閉三方向切断の全350式を
  (W1,W2)*(W1,W2,eta)+((9a−8)(9a+8)eta²)
  と同定する producer PASS12.31s。全350↔6係数、primary分解、square-zero A-free3 nilmoduleを検査。負軸は既知T40と全150係数一致、比1440/41。正の+8/9軸は既存成分未同定。親全環の新しい非平坦性や新成分は主張しない。独立監査未完了。

## まだ公開 continuation packet にできていない途中作業

### J=0 非平方 coker

`work/nonsquare-J0-pid-v63/` に保存。公開 J0 checkpoint は切断時点で未作成。

- new code: `code/s61_nonsquare_J0_pid_prepare_20261103_v63.cpp` と同系 verify / smith / rowgcd Python。
- safe inputは `work/nonsquare-v62/cubic27-v6/after21` (36×70)。ad57 basisと、aFBCだけを反転する21 pivotから来る。Jを反転した27×61をJ0へ代入してはいけない。
- J0-matrix.txt: 3.8MB、108530項、SHA256
  `4e8ed3e3c8eb89b614b15ab9a11c91ff51b90e5ef6cca826a27043eb6167e683`。
- d=(6e+4)/(5a)の全2520成分をQQで独立照合。aの冪だけをrowごとに払った。unit setは a,F,Bs,C,3e+2、B=2Bs/(5a)。
- 既知単元13追加pivotで23×57、37.9MBの J0-residual.txt を作成した。これは coker 消滅の証明ではない。
- `QQa-e-smith/`, `QQe-a-smith/` はそれぞれ180CPU/約190wallでtimeout、input readyのみ。Smith / ann / unit出力なし。
- 初期 `QQe-a-rowgcd/` は Singular の rational-function polynomial division / extgcd backend failure。receiptは arithmetic_errorへ訂正済み。ログの NONUNIT_FACTOR を数学的結果として使ってはいけない。
- `rowgcd-rational-smoke-v3/` は独自 VerifiedDivRem/VerifiedExtGCD fixture と replay がPASS。
- `QQe-a-rowgcd-verified/` の実計算は120CPU/約130wallでtimeout。11 stepは完了確認、12番目の witness は proof.sing (4261B) に保存されたが、最後のmatrix update / replayは未確認。
- `checkpoint-replay.sing` を開始した後に接続断。FIRST read `checkpoint-replay.log` / JSON、`resume-matrix.sing`。12step replay成功を推測しない。
- 有限欠損、実非単元曲線、全余核消滅、全四次被覆はいずれも未証明。失敗と成功した部分恒等式を分けて公開保存する。

### V8 lookahead (未archive)

`work/nonsquare-v62/product-field-v7/low-image-D3-K10.*`。

- degree3 multiplierでtarget aJを保ち、他行像degree≤10とするexact QQ証人。3779×921係数行列、897-square Dixon候補を全QQ式で確認。
- redundant columnをappendしてunit pivot:13×46、167410項。これは invertible replacement と区別する。
- k12の場合は true replacementで13×45、186299項。
- k8の4915×921係数行列はmodular full columnrank921、QQでこの有界解なし。k9未検査。
- LLLはconditioning変更だけで新しい全点結果なし。
- V8 archive / manifest の作成コマンドは接続断で開始できなかった。既存immutable v6はJにGitHub保存済み。

### 未開始・途中の認識実験

- +8/9軸のT01 regular-nilpotent認識は未開始。CA/CB外という説明も未監査・未昇格。
- T39/T40 extended arc sourceのspecial-a障害は `work/q-mps-periods-v63/cycle-source-small.{py,log,json}` に途中データ。all120-label cofactor follow-throughは未完了で、最後の作成コマンドが実行されたかも未確認。既存JのQa/T39,T40 graph結果を上書きしない。

## 復旧後に必要な具体的操作

1. **同じ既存Cloud workspaceを復旧・再接続する。** 新しい空workspaceやcloneで未push分を上書きしない。runtime接続の復旧が必要であり、成功済みセットアップの再実行は不要。
2. `cd /workspace/poisson-p4`、`source /workspace/cloud-setup/activate.sh`。
   旧sourceをimportする再現は `PYTHONDONTWRITEBYTECODE=1` を使う。
3. `git status --short --branch` と `git rev-parse HEAD` を確認してから `git fetch origin main`。この障害記録はAPI追加コミットなのでローカルmainよりremoteが進んでいる。diffを確認し、未pushファイルを保持して `git merge --ff-only origin/main` で取り込む。reset / clean / force-pushはしない。
4. 全frozen reportのinput/evidence SHAを確認。上記の未観測replay結果を読む。J0とV8の再現可能な途中状態をcurated archive / 日本語ノートとして先に保存する。
5. B16 / extended4 / bounded-unit-liftのaccepted証拠を保全し、QQimplicit17 / Hensel16 / CA/CB/N三packet / weight6二軸の独立監査を終える。未完了のreviewをPASS扱いしない。
6. `research/CLAIMS.json`, COMPONENTS, STATE, NEXT_TASKS, v63 handoff, `papers/RESEARCH_SOURCE_INDEX_v63.json` を同期。存在する未push証拠だけを選んでstageし、`git diff --cached --check` と `.venv/bin/python scripts/verify.py --suite integrity` を行う。
7. 通常の `scripts/checkpoint.py --message ... --push` で main を保存する。その後も未解決研究を継続する。ユーザーの既存のmain保存承認は継続しているので再承認を求める必要はない。
8. 初回から存在するuntracked `reports/20261102T142937.117200Z-all/` は無関係であり、巻き込まない。不要な旧native ELFもstageしない。

全一般アフィンfiberの五点尽くし、成分総数の上界、残る非平方cubic/quartic支持、全局所枝・全大域primary構造は未解決である。今回の中断理由は実行環境の接続断であり、研究を完了したからではない。
