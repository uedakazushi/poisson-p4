# Poisson structures on projective four-space

2026-11-09 / v66-conceptual-structure

## 主目標：列挙の検算から、構造・復元・変形による理解へ

既存分類原稿・台帳は147成分（一般rank4が141、rank2が6）を主張する。PGL5商前の射影Jacobiスキームの被約既約成分が対象である。
現在は、個別計算を共通の構造定理へ置き換え、なぜその族・次元・境界が現れるかを説明することを主軸とする。検証はそのための補助線であり、監査報告の増加自体を目標にしない。

最初に [研究方針](RESEARCH_POLICY.md)、次に [構造整理ノート](research/conceptual/NOTES.md)、[現在地](research/classification/CURRENT.md)、[次の課題](research/NEXT_TASKS.md) を読む。Codexの入口は [AGENTS.md](AGENTS.md)。

## 二つの構造原理

新原稿： [Poisson退化因子が生成する可換対称性](papers/poisson_divisor_symmetries_ja.tex)。日本語LuaLaTeX、現在7頁。

1. P^(2n)の一般に非退化なPoisson構造について、退化因子の被約既約因子数rから、r−1次元の標準的な可換Poisson対称性代数を作る。非被約・非SNCも許す。P4では射影線形軌道の次元は25−r以下。
2. 固定した重複度付き反標準因子F上のrank4構造は、閉有理二形式の線形条件とPfaffian二次条件で正確に記述できる。

通常の証明を原稿に収録した。基本原理自体の新規性は未確認。これだけで147成分の全分類や全独立監査を完了したとは主張しない。
次の具体的な目標は、H+2Qの複数の線形核計算を一つの有界極・次数制御の補題へ置き換え、二因子型のJordan/共鳴の入口を統合すること。

## 既存の数学を保持する

完全分類主稿は `papers/poisson_p4_classification_ja.tex` と同名PDF。今回その368頁本文、成分・主張台帳、元証明書は変更していない。
`research/classification/GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md` と `REPRODUCE_EXACT147.md` は既存分類の入口。
`research/publication/REVIEW.md` は必要な独立監査の補助キュー、旧classification/PLAN.jsonは分類構築の完了記録である。
TRIPLE→Ga9の独立検算は `code/verify_final_bridge.py` に保持する。具体的理由なしに全旧実験やsetupを反復しない。

全PGL軌道、全非被約厚み、精密nil/H0/T2は現在の主目標に含めない。新原稿のPDFは配布・Library用に生成し、Gitには再生成可能なTeXを保存する。
