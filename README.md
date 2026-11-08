# Poisson structures on projective four-space

2026-11-09 / v67-global-reconstruction

## 現在の主作業

**全分類を「標識付き復元と持上げ可能な変形」によって再構成する。**
H+2Qの局所補題や検算の反復を最優先にはしない。
現行分類稿・台帳は147被約既約成分（generic rank4:141、rank2:6）を主張するが、今回全証明を独立検証したという意味ではない。
対象は非零Poisson構造の射影化、PGL5商前。全軌道や全非被約構造は対象外。

再開は [研究方針](RESEARCH_POLICY.md) → [全体計画](research/conceptual/GLOBAL_PLAN.md) →
[構造原稿](papers/poisson_divisor_symmetries_ja.tex) → [研究ノート](research/conceptual/NOTES.md)。
Codexへの具体的な指示は [AGENTS.md](AGENTS.md)。

## 全体を整理する入口

rank4は「超平面因子あり」「超平面なしの二次＋三次」「被約既約五次」に分かれる。
超平面を標識した全群は一意なaffine quadratic Poisson模型に入り、f=Pf P、Z=div Pを固定すれば
`div P=Z, P#df=fZ`（線形）と `Pf P=f != 0`（二次）で正確に復元できる。
共通の共鳴、相対逆形式、一次・二次持上げ、グラフ閉包、因子標識の重複を統合した。
rank2は飽和foliationと零因子として別に扱う。通常の証明を構造原稿に収録したが、新規性は未確認。

最初の実作業は、単純H+四次と重複Hの二つ以上の旧入口を共通命題へ置き換え、主稿に統合すること。
固定因子の分類だけで終わらず、実変形と極大閉像の判定まで接続する。
H+2Qの0/1/6余剰は必要時の局所課題として保持する。

## 保存済み資料と検証の射程

全分類主稿 `papers/poisson_p4_classification_ja.tex` と同名PDF、成分・主張台帳、元証明書は今回変更していない。
`research/classification/GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md` は既存分類の論理入口。
旧classification/PLAN.jsonは構築の完了記録で、現在の作業キューはGLOBAL_PLAN.md。

```sh
python code/check_global_structure.py --output work/global-structure.json
mkdir -p work/structure-build
lualatex -interaction=nonstopmode -halt-on-error -output-directory=work/structure-build papers/poisson_divisor_symmetries_ja.tex
lualatex -interaction=nonstopmode -halt-on-error -output-directory=work/structure-build papers/poisson_divisor_symmetries_ja.tex
```

恒等式検算は29項目、全分類の自動証明ではない。今回の構造研究bundleは自己完結した構造原稿・方針・ノートを含むが、全147成分の証明bundleの代替ではない。
