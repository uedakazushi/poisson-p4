# GitHub mainからの手動再現

2026-11-11 / v18。GitHub mainが唯一の正本。生成PDF、LaTeX中間生成物、旧ZIP、キャッシュをGit追跡せず、GitHub Actionsも利用しない。Library保存は実際のpush失敗時のみ緊急退避とする。

## 環境

Python 3、SymPy 1.14.0、LuaLaTeX（ltjsarticle、luatexja-preset、HaranoAji、amsmath、mathtools、geometry、hyperref、xurl、needspace）、Poppler（pdfinfo、pdftotext）を要する。VS Code LaTeX Workshopはrootの `.latexmkrc` と `.vscode/settings.json` を使用。追加の147成分分類検算にはrootの依存設定と各producerで指定するCAS・入力が必要。外部資料PDFやLibraryの古いbundleは現行incidence検証の入力には使わない。

## 最優先：Poisson不変部分多様体incidence

```bash
bash scripts/verify_github_source.sh incidence
```

`scripts/reproduce_subvariety_incidence.sh` は現在の **13 producerと13保存JSON** を照合し、`papers/poisson_subvariety_incidence_ja.tex` をLuaLaTeXで二回組版する。出力はGit非追跡の `work/` のみ。

最も新しい証明は `papers/del_pezzo_pfaffian_incidence.tex`。反標準零曲線Cに沿うPfaffianの二重消滅（必要五次空間41次元）、余法線Lie括弧κとPfaffianの法線一次ジェットとの積公式、κ=0とD∈I_S²の同値、可換型の必要五次空間15次元、非可換型の可換化集合がCに含まれること、滑らかなCで有限可換化スキーム長16という幾何学的主張は本文に証明する。

`code/subvariety_incidence/check_dp_pfaffian_conormal.py` と `check_dp_conormal_lie_pfaffian.py` はそれぞれHilbert空間のrankおよびPfaffianの局所一次ジェット・次数数値を正確に検査する。最新の保存証明書は `research/subvariety_incidence/certificates/dp_pfaffian_conormal.json` と `dp_conormal_lie_pfaffian.json`。これらの計算は幾何学的な全域的証明の代替でも、rank4延長存在の証明でもない。

```bash
python3 code/subvariety_incidence/check_dp_conormal_lie_pfaffian.py \
  --output work/dp_conormal_lie_pfaffian.json \
  --check research/subvariety_incidence/certificates/dp_conormal_lie_pfaffian.json
latexmk papers/poisson_subvariety_incidence_ja.tex
```

## 独立副研究：完全等方Log点

```bash
bash scripts/verify_github_source.sh tangent-cone
```

現行日本語稿は `papers/isotropic_log_tangent_cone_ja.tex`、その共有TeXは `papers/isotropic_nilpotent_pencil.tex`。全二次初期関係7本、真の非被約接錐、局所次元28、交叉次元23/24、冪零場・C*重みの独立producerと証明書を検査する。元の局所環の被約性や三次以上の完全初期イデアルは未解決であり、主incidence研究より優先しない。

## 全147成分分類の既存preflight

```bash
bash scripts/verify_github_source.sh full
```

分類の現行原稿・既存厳密計算資産、指定された高リスク内部査読の範囲は `research/classification/CURRENT.md`、`research/publication/REFEREE_REPORT.md`、`research/reproduction/README.md` を参照。fullスクリプトの成功のみで全分類の原証明・新規性・外部査読が完了したとは扱わない。

## 実行状態の記録方法

- 本文に証明があること、独立に今回検証したこと、SymPyで確認したこと、LuaLaTeXの現行稿を組版したこと、GitHub mainへpushしてリモート読戻ししたことは異なる状態として扱う。
- 数学的主張に誤りが判明したら関連する論文・方針・ノート・証明書を一貫して修正する。未実行の全producerや全稿の組版を実行済みと報告しない。
- 作業履歴や古いbundleをこの現行再現手順に蓄積しない。
