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

## 新しい部分多様体incidence：完全な手動再生成

```bash
bash scripts/verify_github_source.sh incidence
```

GitHubの `papers/poisson_subvariety_incidence_ja.tex`、`papers/quadric_log_extension.tex`、`papers/quadric_log_orthogonality.tex` と `code/subvariety_incidence/*.py` を用いる。
`research/subvariety_incidence/certificates/*.json` をexactに再計算して同一内容を照合し、一般剛性の二つの非零minor（70875/88、4443337729552/7688671875）、rank60/syzygy90、101項目の明示恒等式を検証。第三のproducerはquadricのPfaffian二例と一般Ratのstrict so(5) stabilizerを消すrank10 minor 27980515917168640を検証。追加のLog証明書には直交性wallの全パラメータJacobi計算と、さらに5つの非退化二次形式（c=0、等方bを含む）のstrict stabilizer・Pfaffian検算を収録する。全6 JSONは対応する6 producerから再計算・照合する。
LuaLaTeXを二回実行して現行の論文PDFを生成し、15頁以上であること、未解決参照やOverfullがないこと、抽出テキストが空でないことを検証する。抽出テキストのSHA256を出力し、PDFバイナリのbitwise再現を要求しない。
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

**今回確認した範囲**：新しいincidence原稿と2本の正確な計算コードは元の局所実行資産と一致し、
証明書の再生成・LuaLaTeX11頁・PDF抽出テキスト一致を局所で確認した。
一方、本実行環境ではgithub.comのDNS接続ができずfresh `git clone`で147分類の
全preflightを今このセッションで改めて実行したとは主張しない。
既存v74の検査記録はその射程で再利用する。
