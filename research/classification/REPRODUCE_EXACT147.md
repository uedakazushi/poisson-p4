# 正確147成分分類：証明と最小再現入口

2026-11-08 / 147 components / 1158 claims。

全支持・全成分性・相互相違・一般元の定理は `GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md`。無条件の別作者最終査読は `FINAL_147_COMPONENT_SUPPORT_EXHAUSTION_COMPLETION_ADDENDUM.md`。最終manifestは `inputs/global_exact147_component_classification_final_proof_manifest_20261108.json` と `inputs/final147_component_support_exhaustion_completion_manifest_20261108.json`。主稿は `papers/poisson_p4_classification_ja.tex`、PDFは同stem。全147の元source・証明・証明書・引用入口は `research/COMPONENTS.json` と `research/CLAIMS.json`。

## 証明の最短依存

1. 元full150/Euler/350及びfaithfulS61全134の同値、rank2全支持、既各成分証明を既監査で接続。
2. `KNOWN142_SEVEN_CANDIDATE_WHOLE_POISSON_COVER_AUDIT.md` とroot separate review：全factor partition/全modularJordan-coker/lower/zero-unit/非飽和thetaから既142＋七閉候補の全支持有限閉被覆。
3. P*の真T19曲線で全h5を既知成分へ、候補一つを除く。`UNIMODULAR_H5_T19_AND_GLOBAL_SIX_CANDIDATE_COVER_SEPARATE_REVIEW.md`。
4. `SIX_CANDIDATE_PROPER_PARENT_FLAGS_SEPARATE_REVIEW.md`、`T37_SIX_TARGET_LEAF_OBSTRUCTIONS_SEPARATE_REVIEW.md`、全prior142各projection、全候補相互非包含から五つの新23全成分。`FIVE_CANDIDATES_FULL_COMPONENTS_AND_ONE_ORBIT_FRONTIER.md` と別作者catalogue監査。
5. 最後constantTRIPLEの全actualPoisson族が既Ga9constant23へ。`CONSTANT_MODULAR_TRIPLE_GA9_FULL_SOURCE_SEPARATE_REVIEW.md` とroot全係数再生。全rank四閉包を有限和へ入れてからrank二closure項を戻し、exact147。

この順序に循環は無い。原sourceのlocalized unitと全zero補集合はそれぞれのimmutable note/manifestへ保持。後続のより強い結果を、旧noteの歴史的『未完』へ遡って書き換えない。

## 環境と軽量検査

既導入Cloud環境を再利用する。M2 1.24.11（引継ぎ元1.26.05との必要QQ互換確認済み）、`.venv`、Fraction/SymPy、LuaLaTeX。setup全体や全legacy実験の反復は不要。

```sh
python3 scripts/classification.py check
python3 -m unittest discover -s tests -p test_classification.py
```

これらはキュー/制御コードを検査し、数学定理を自動的に証明するコマンドではない。元run前のgate PASSは各契約/receiptへ保存済み。現在Eはdoneなので、新規run用gateはTask is not active/readyと正しく拒否する。最終保存時の契約hash/予算検査は、完了taskの読取専用historicalcontextで行い、実キューを再活性化しない。

## 最後の全係数族を再生する

rootの独立replayは新wrapperを必要とせず、固定caller JSONに全sourceを保存している。

- `inputs/constant_modular_triple_ga9_full_family_root_review_caller_20261108.json`：元積Poissonから全六、aux両逆/det、canonical150、原350、actualPf126、e0、Ga9 actualframe/Phi_u。
- `inputs/unimodular_pole_t19_full_same_curve_root_review_caller_20261108.json`：元completedT19から同じ固定curve全150/350/126、両逆、固定符号と旧universal reportの全係数一致。

元reportとcallerはexclusiveであり上書きしない。例えばrepo rootで次を実行すると、**出力先だけ**を未使用workへ変えた同じ数学計算を再生する。inputhashは全て元のまま、元callerのsourcehashも先に検証する。再生reportには実際に走ったsource全bytesが保存されるので、I/O差分を検査出来る。

```python
from pathlib import Path
import hashlib, json, subprocess, sys
caller = Path('research/classification/inputs/constant_modular_triple_ga9_full_family_root_review_caller_20261108.json')
data = json.loads(caller.read_text())
source = data['inline_source']
assert hashlib.sha256(source.encode()).hexdigest() == data['inline_source_sha256']
old = "OUT=ROOT/'reports/CLASSIFICATION_CONSTANT_MODULAR_TRIPLE_GA9_FULL_FAMILY_ROOT_REVIEW_20261108.json'"
new = "OUT=ROOT/'work/replay147/triple_ga9.json'"
assert source.count(old) == 1
work = Path('work/replay147')
work.mkdir(parents=True, exist_ok=True)
script = work / 'triple_ga9.py'
assert not script.exists() and not (work / 'triple_ga9.json').exists()
script.write_text(source.replace(old, new, 1))
subprocess.run([sys.executable, str(script)], check=True)
```

`.venv/bin/python` でこのPythonを実行する。root同族の実測は0.315秒/RSS50072KiB、上限60秒/1024MiB。P*同曲線は2.827秒/RSS54360KiB、同予算。再生の時刻/RSS/埋込I/O pathは元reportとbyte一致しないが、全恒等式と全係数順は固定される。既存exact証明書の入力SHAを数学変更のために改訂しない。

## 歴史的入力と証明範囲

旧台帳のprojectionは保存されたGitblob又は当時のcommitに固定してある。current147台帳で旧142/126/97 sourcehashを置換しない。普通のGit cloneはその履歴も保持する。途中の未commit台帳三Gitblobは `inputs/exact147_historical_input_blob_archive_manifest_20261108.json` と隣接exactbytes snapshotへ保存したので、同じGitblobが通常pushで到達可能となる。旧manifestは変更しない。新五metadataの別作者projectionは `inputs/five_new_components_catalogue_audited_projection_20261108.json`。

旧point/profile/高次nil/H0/T2等secondary実験は最終分類の再実行義務ではない。一件の旧未commit review消失と二optional legacy runtime/gate制約は別receiptに保存済みで、現採録の必要証明には依存しない。全scheme被約性、精密厚み、全orbit分類、文献上の新規性はこの再現手順又は主定理から追加主張しない。

最終二回組版・参照/引用/box・可視確認・control tests・inputbinding検査は `CLASSIFICATION_MAIN_MANUSCRIPT_COMPONENT147_CLAIM1158_*BUILD_20261108.json` に別保存する。GitHubへの通常main pushの実状態は `git rev-parse HEAD` と `git ls-remote origin refs/heads/main` で照合する。数学・組版・保存は別の状態である。
