# GitHub main だけからの再現

## PDFs are build outputs (2026-11-10)

**GitHub tracks sources, not generated PDFs.** All 7 formerly tracked PDF files have a corresponding tracked TeX source. They have been removed from Git tracking, not deleted from history. The mathematical text, formulae and proofs are retained in the TeX sources; generated PDF output is ignored. Build from the root of the repository:

```bash
latexmk papers/poisson_p4_classification_ja.tex
bash scripts/build_papers_latexmk.sh
```

The repository `.latexmkrc` always selects LuaLaTeX, including for `latexmk -pdf`. Outputs appear in `work/latexmk/` and never require ChatGPT Library. The formerly tracked PDF/source pairs are:

| Previously tracked PDF | Canonical source |
|---|---|
| `papers/ca_cb_local_normalization_ja.pdf` | `papers/ca_cb_local_normalization_ja.tex` |
| `papers/normal_actual_divisor_recognition_ja.pdf` | `papers/normal_actual_divisor_recognition_ja.tex` |
| `papers/poisson_divisor_symmetries_ja.pdf` | `papers/poisson_divisor_symmetries_ja.tex` |
| `papers/poisson_p4_classification_ja.pdf` | `papers/poisson_p4_classification_ja.tex` |
| `papers/positive_weight6_component_ja.pdf` | `papers/positive_weight6_component_ja.tex` |
| `papers/s61_first_jet_reduction_ja.pdf` | `papers/s61_first_jet_reduction_ja.tex` |
| `research/notes/p4_poisson_status_strategy_ja_20261030_v1.pdf` | `research/notes/p4_poisson_status_strategy_ja_20261030_v1.tex` |

Prior delivery/build receipts (e.g. `reports/final_manuscript_validation.json`, `reports/CA_CB_LOCAL_NORMALIZATION_PAPER_20261103_v63/REPORT.md`, `reports/MANUSCRIPT_SYNC_20261103_v62.json`) document the prior PDFs but are **not requirements** to retrieve them. A historical script that expects `papers/foo.pdf` can first run `latexmk -outdir=papers papers/foo.tex` (generated output is ignored). PDF-only mathematical content must never be discarded merely to shrink Git. The source-only current bundle generator includes all standalone TeX, the old evaluation note TeX, mathematical proofs and necessary proof input, but does not require prebuilt PDF files.



2026-11-10。**GitHub mainが唯一の正本**。過去のLibrary ZIP、会話、未追跡のwork/は再現・再開の前提にしない。GitHub Actionsは使用しない。

## 取得と環境

```bash
git clone https://github.com/uedakazushi/poisson-p4.git
cd poisson-p4
git switch main
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
source .venv/bin/activate
```

Python3、SymPy 1.14.0、NumPy 2.3.5、mpmath 1.3.0（root requirements.txt）、LuaLaTeX（ltjsarticle、luatexja、HaranoAji、amsmath、mathtools、geometry、hyperref、xurl、needspace等）、Popplerのpdfinfo/pdftotextが必要。一部の既存分類の有限証明書を完全再生成する場合には、各producerに指定されたMacaulay2/Singular/FLINT等のCASも手動で用意する。Actionsや自動CI環境は使わない。

## 完全等方Log点の局所次元・長鎖成分同定・非CM性（v11）

`code/subvariety_incidence/check_local_dimension28.py` はSymPy 1.14.0の有理数演算で、全150一次Jacobi方向を消す二つのSchouten汎関数、30接方向の独立性、二次式 `9*t21*t29`・`-4*t25*t27`、Log(1,1,1,2)族の微分階数28を示す小行列式 `-57177414000000`、rank4曲線を検証する。証明書は `research/subvariety_incidence/certificates/local_dimension28.json`、本文は `papers/quadric_log_component_intersection.tex`。局所次元28は147成分完全性を使わずに得られる。

**新しい11番目のproducer** `code/subvariety_incidence/check_rank4_long_nonconstant_membership.py` は同次五変数の全Jacobi、アフィン六Poisson括弧、modular場のJordan型、`a != 0`の実アフィン線形（従ってPGL5）座標変換を `QQ(a)` 上で検証する。変換先は従来の `RANK_ONE_LONG_NONCONSTANT_COMPONENT_ENTRANCE.md` の非定数g全sourceに一致し、残余二次式の半Hessian determinantは `-1/(32*a**4)`。計算のJSONは `research/subvariety_incidence/certificates/rank4_long_nonconstant_membership.json`、日本語証明は `papers/rank_four_long_component_membership.tex`。既存の25次元rank4成分性および28次元rank2成分性の原証明を入力として、pの局所環が非等次元・非Cohen--Macaulay・非Gorenstein・非局所完全交叉と結論する。ambient非被約性、全接錐分解、他のrank4成分の有無は未証明。

全11 producerは `bash scripts/verify_github_source.sh incidence` で手動再生成・保存済JSON照合する。GitHub Actionsは使わない。生成PDFはgitignoreされ、Library取得は再現の前提でない。

## 新しい部分多様体incidence：完全な手動再生成

```bash
bash scripts/verify_github_source.sh incidence
```

GitHubの `papers/poisson_subvariety_incidence_ja.tex`、`papers/quadric_log_extension.tex`、`papers/quadric_log_orthogonality.tex` と `code/subvariety_incidence/*.py` を用いる。
`research/subvariety_incidence/certificates/*.json` （現行11種）をexactに再計算して同一内容を照合し、一般剛性の二つの非零minor（70875/88、4443337729552/7688671875）、rank60/syzygy90、101項目の明示恒等式を検証。第3のproducerはquadricのPfaffian二例と一般Ratのstrict so(5) stabilizerを消すrank10 minor 27980515917168640を検証。追加のLog証明書には直交性wallの全パラメータJacobi計算、5つのrank2 H例（c=0・等方bを含む）、さらにrank1 Hの2つの正規形（c=0,1）のstrict stabilizerと全Jacobi/Pfaffian検算を収録する。全11 JSONを11 producerから再計算・照合する。追加のisotropic境界証明書はq|Wが階数1/0の三例のstrict stabilizer、全Jacobi理想(0)、(a(b-4))、(ac,c²)、Pfaffianを検査し、scheme非被約性と被約階数4枝を区別する。
LuaLaTeXを二回実行して現行の論文PDFを生成し、15頁以上であること、未解決参照やOverfullがないこと、抽出テキストが空でないことを検証する。抽出テキストのSHA256を出力し、PDFバイナリのbitwise再現を要求しない。
さらに9番目のproducer `code/subvariety_incidence/check_ambient_quadratic_obstruction.py` は、全150個の周囲斉次二次方向の射影Jacobi線形化と5係数の二次障害汎関数を独立に構成する。rank95／augmented rank96、障害9、接空間次元30を `research/subvariety_incidence/certificates/ambient_quadratic_obstruction.json` と照合する。これはambient非被約性の証明ではない。
PDFのバイナリバイト列の同一性までは主張しない。
出力はgitignore済みの `work/subvariety-incidence-check/` のみ。

原稿の一般del Pezzo剛性は一般の非空Zariski開で証明され、全特殊例まで証明したものではない。
ここでの検査は147成分完全分類の独立証明ではない。

## 既存147成分分類を含む手動preflight

```bash
bash scripts/verify_github_source.sh full
```

必須入力のhash監査、147成分カタログ、グローバル恒等式、分類主稿・構造稿・正規因子稿・一次還元稿の四稿の二回組版、
incidence研究の全検査を行う。旧223 primary claimの全てを改めて独立に証明するコマンドではない。
既存の限定計算の原入力・producer・文献は `research/reproduction/README.md`、
`research/publication/REFEREE_REPORT.md` と関連する証明・証明書に整理されている。

## 現行一式をGitHubだけからローカル生成（任意）

```bash
python scripts/make_structural_bundle.py --output work/current-github-source-bundle
```

追跡済みの原稿・証明書・入力を、SHA256 manifest付きの分割ZIPに再生成する。
このコマンドは**Libraryへの保存や取得を行わない**。出力を別ディレクトリに展開して
manifestと全ファイルのhashを検査できる。GitHub mainの実ファイルが正本であり、
ZIPの旧版・キャッシュ・入れ子ZIPを正本として扱わない。

## 将来の更新ルール（必須）

1. 研究方針→ノート→論文の順に読み、証明済み／条件付き／計算確認／未証明／既知を分離する。
2. 必要な**最新の**TeX・証明・証明書・コード・入力・再生成手順をすべてGitHub mainへcommitする。
   Library、会話、work/や旧ZIPだけに必要な資料を残さない。
3. 変更した数学・組版・依存ファイルの検証を手動実行して範囲を記録する。
4. 他者変更を保護し、期待したリモートheadを確認してmainにnon-force push、
   最後にリモート先端のSHAを読み戻す。ローカルcommitだけで「保存成功」とはしない。
5. **GitHub Actionsを作成・dispatch・実行しない。** hosted CIも依存しない。
6. **Libraryへの通常保存はしない。** 実際にGitHub pushが失敗した場合だけ、
   現行の自己完結bundleを緊急退避として保存してよい。GitHub未保存を正直に報告し、
   次にmainへの復旧を優先する。既存Libraryファイルを勝手に削除しない。
7. 過去の版、バックアップ、作業ログを現行一式に蓄積せず、今の研究内容だけで再開可能にする。

**今回の検証範囲**：長鎖成分同定の全六括弧・同次Jacobi・Jordan型・Jacobian・二次因子のQQ(a)恒等式は独立のSymPy 1.14.0局所実行でPASS。既存の局所次元28コードも再実行でPASS。一方、GitHubへの直接DNS接続ができず、fresh cloneからの全11 producer照合や主稿全体のLuaLaTeX組版はこの作業時に再実行していない。過去のv74の分類監査や従来PDF組版receiptを今回実行したものとして扱わない。GitHub mainの実SHAを読み戻して保存を確認する。
