# Poisson structures on projective four-space

2026-11-09 / v71-end-to-end-completion

## 最終目標と実行指示

**複素P4上のPoisson構造の被約既約成分を、共通の幾何的認識・復元・変形から導く完全分類論文として完成させ、再現可能な現行一式をmainとLibraryへ保存する。**
対象は非零構造を射影化したJacobiスキーム、PGL5商前。
現行の結論は147成分（一般rank4が141、rank2が6）。この数を仮定にせず、完全性・成分性・相違と一般元の記述を証明して導く。

再開順：
[研究方針](RESEARCH_POLICY.md) → [具体的な最終定理・終了条件](research/conceptual/FINAL_TARGET.md) →
[全工程の実行経路](research/conceptual/GLOBAL_PLAN.md) → [現在の未達条件](research/classification/CURRENT.md)。
Codexの継続命令は [AGENTS.md](AGENTS.md)。

**一つの補題・工程、検算成功、PDF、commit/pushを理由に終了しない。**
非plane導手、rank2、極大閉像、論文統合、文献・検証、bundle、保存までを一つの課題として自律的に続行する。
中間保存は行い、その後すぐ次へ進む。通常の終了はFINAL_TARGETの全条件達成時だけとする。
実際の実行制約による不可避な中断では未完を明記し、完成や無人の継続を主張しない。

## 現在有効な数学

構造原稿 `papers/poisson_divisor_symmetries_ja.tex` と分類主稿 `papers/poisson_p4_classification_ja.tex` は、
marked_hyperplane_reconstruction、poisson_lifting_examples、quadric_cubic_reconstruction、
conductor_conormal_reconstructionの共有TeXを使う。
Aの標識/modular、Bのcofactorと磁気留数、Cのplane余法線を統合済み。
非plane・R2・全極大閉像への概念的接続と最終原稿の完成は残る。
今回の指示変更では数学本文・成分台帳・旧証明書を変更していない。

検証入口：
```sh
python code/check_global_structure.py --output work/global-structure.json
python code/check_marked_reconstruction.py --output work/marked-reconstruction.json
python code/check_quadric_residue_reconstruction.py --output work/quadric-residue-checks.json
python code/check_conductor_conormal_reconstruction.py --output work/conductor-conormal-checks.json
bash scripts/build_structural_manuscripts.sh
```
これらの成功は全分類の自動証明ではない。必要な原証明・適用範囲・証明書も確認する。
全軌道・全交叉・全scheme厚み、任意次元への完全分類、論文の固定ページ数は最終目標に含めない。
