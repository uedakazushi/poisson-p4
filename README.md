# Poisson structures on projective four-space

2026-11-09 / v68-marked-reconstruction

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

単純H+四次と重複Hの標識・modular・共鳴の入口を共通証明へ置換し、主稿に統合した。
新しい原始関数式 `df + i_Z α = 0` はPfaffianで割らずにJacobiを判定し、次数低下とtranslation余核を保つ。
二単純H+三次の全模型でも適用を確認した。個別源の帰属・成分性・非包含の証明は保持する。
TRIPLE→Ga9では次数二以下を保つ補正を導き、指定因子族の二次持上げに補正が必要なことを示した。
可換場による二重Hの分裂族では補正なしで全次数を持ち上げられる。
次は超平面なしQ+Cのcanonical symmetryと全31代表の復元を統合する。

## 保存済み資料と検証の射程

全分類主稿 `papers/poisson_p4_classification_ja.tex` と構造原稿は、
`papers/marked_hyperplane_reconstruction.tex` と `papers/poisson_lifting_examples.tex` を共有する。
主稿の既存の全証明を保持し、共通部分だけを置換した。成分・主張台帳と旧証明書は未変更。
全147 IDの一般実因子・入口・原証明への案内は `research/conceptual/component_entrances.csv`。
一般元はA:113 / B:10 / C:18 / R2:6であり、成分閉包の分割ではない。
`research/classification/GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md` は既存分類の論理入口。
旧classification/PLAN.jsonは構築の完了記録で、現在の作業キューはGLOBAL_PLAN.md。

```sh
python code/check_global_structure.py --output work/global-structure.json
python code/check_marked_reconstruction.py --output work/marked-reconstruction.json
python code/map_component_entrances.py --output work/component_entrances.csv
mkdir -p work/structure-build
lualatex -interaction=nonstopmode -halt-on-error -output-directory=work/structure-build papers/poisson_divisor_symmetries_ja.tex
lualatex -interaction=nonstopmode -halt-on-error -output-directory=work/structure-build papers/poisson_divisor_symmetries_ja.tex
```

既存29項目に加え、共通恒等式と表示した全パラメータ族・正逆を検算した。全分類の自動証明ではない。
環境内にTeX一式が別展開されている場合は `scripts/build_structural_manuscripts.sh` を使う。
今回の検証範囲と組版結果は `reports/marked_reconstruction_validation.json`。
