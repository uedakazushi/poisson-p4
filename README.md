# Poisson structures on projective four-space

**Source policy (2026-11-10):** Standalone papers and research notes are tracked as **TeX, not PDFs**. The 7 previously tracked PDFs have corresponding TeX sources and reproducible build instructions; use `.latexmkrc` with `latexmk papers/poisson_p4_classification_ja.tex` or `bash scripts/build_papers_latexmk.sh`. Generated PDFs are ignored in `work/latexmk/`. Primary proofs, exact calculation inputs and certificates remain tracked. No GitHub Actions; GitHub main is canonical; ChatGPT Library only on an actual failed push.


2026-11-10 / v10-incidence + v74-classification

## 最新の正本と再現方法（2026-11-10）

**GitHub mainが唯一の研究正本です。** 過去のLibrary ZIPや会話は必要ありません。
GitHub Actionsは一切使わず、手動スクリプトで原稿と証明書を再現します。
通常のLibrary保存は行わず、GitHub pushが実際に失敗した場合にだけ緊急避難保存します。

- [GitHubのみからの再現手順](REPRODUCE_GITHUB.md)
- [部分多様体incidence現行論文v4（Q3一般三成分と対数境界の直交性wall）](papers/poisson_subvariety_incidence_ja.tex)
- [対応する研究方針](research/subvariety_incidence/RESEARCH_POLICY.md)・[ノート](research/subvariety_incidence/NOTES.md)
- [incidenceの正確な計算証明書](research/subvariety_incidence/certificates)
- [完全等方Log点の局所次元28・rank2/rank4成分交叉](papers/quadric_log_component_intersection.tex)：rank30の接空間上の二つの互いに素な障害式と、28次元対数型源の非零minorで独立に証明。rank4成分IDとambient非被約性は未解決。
- ambient二次接錐：完全等方Log点では接空間が30次元、局所次元は29以下。固定Q延長三方向の二次障害は-9c²。現行カタログの36次元・39次元rank2成分はこの点を含まない（証明はpapers/quadric_ambient_quadratic_obstruction.tex）。ambient非被約性は未解決。
- 再現実行: `bash scripts/verify_github_source.sh incidence`
- 既存分類を含む検査: `bash scripts/verify_github_source.sh full`

この現在の方針が以下の過去の保存報告に優先します。

## 現在の状態と次の研究

複素P4の射影Poissonスキーム（PGL5商前）の被約成分の**完全分類を記した現行稿**がある。
分類の主張は147成分（rank4 141、rank2 6）。最新の完成報告はM1–M4、C1–C2、P1、R1、S1を満たすと記録し、
主稿418頁、構造稿41頁、正規因子稿91頁、一次還元稿118頁、五分割の再現bundleを保存した。
構築時の完成記録とは別に、以下の公表前査読を行った。

**指定した高リスクの査読義務を充足し、公表候補版として可と判定した。**
固定被約因子・安定rank2・最後五像の証明を追い、51一般族の全パラメータJacobiを別実装で確認した。
一般源の49演算列は214525節点を再生成して一致を確認し、導入と相対系・再現ツールを訂正した。
結果と保証範囲は [REFEREE_REPORT.md](research/publication/REFEREE_REPORT.md)。
LibraryからZIPバイナリを再取得する作業は査読条件から除外した。
保存SHAが一致するローカル五巻の隔離復元、全ファイルのhash確認、四稿の再組版は実施済みと記録されている。
65構文パス・206組の実read/writeを読解し、旧bundleに欠けていた二runnerを追加した。
八支持runnerの不要な文献PDF依存も修正し、隔離実行で元の算術結果と一致した。
最終依存監査は3539論理ファイル・16840記録hash、欠落・不一致・未分類パス0。
現行全bundleの生成・隔離検証・main/Library保存は完了した。
配布元fdec6f8、3599ファイル・二相対リンク、五巻と四PDFの保存記録は
[CURRENT.md](research/classification/CURRENT.md) と [配布receipt](reports/publication_candidate_delivery.json) にある。

再開の順序：
1. [研究方針](RESEARCH_POLICY.md)
2. [現行の査読計画](research/publication/REFEREE_PLAN.md)
3. [現在地・証明依存](research/classification/CURRENT.md)
4. [既定の最終分類定理](research/conceptual/FINAL_TARGET.md)

Codexへの実行指示は [AGENTS.md](AGENTS.md)。
数学本文は papers/poisson_p4_classification_ja.tex と共有TeX、および必要な別稿にある。
147族の明示源と証明参照は papers/component_family_catalogue.tex。
現行再現手順は research/reproduction/README.md。
Libraryの現行4 PDF・五分割bundleは
https://chatgpt.com/space/page_98c9971dc7c48191a2c4b7624fd99fa0 。

この判定は外部査読の採択や優先権の確定を意味せず、旧全計算の独立自動認証とも区別する。
外部投稿・公開は別の明示指示なしに実施しない。
