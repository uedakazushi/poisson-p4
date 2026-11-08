# Poisson structures on projective four-space

2026-11-09 / v65-publication-review

対象は複素P4上の非零Poisson双ベクトルを射影化したJacobiスキームの被約既約成分であり、PGL5商前の分類である。

## 現在地

現行分類原稿・台帳の結論は **147成分（一般階数4が141、階数2が6）**。
主稿は `papers/poisson_p4_classification_ja.tex` と同名PDFで、既存ビルド報告では368頁・1158採録主張。
この原稿の証明収録・リポジトリ内の監査報告と、新たな独立再検証の完了を区別する。
**現在の作業は公表前の重点監査と説明の整理であり、全証明の独立監査は未完了。**
今回、最後のconstant-TRIPLEからGa9への明示族について、旧producerを使わない記号計算で81項目を確認した。
全147成分の被覆・成分性をこの検算だけで認証したとはしない。成分台帳と分類主稿は変更していない。

再開は [研究方針](RESEARCH_POLICY.md) → [現在地](research/classification/CURRENT.md) →
[公表前監査](research/publication/REVIEW.md) → [次の課題](research/NEXT_TASKS.md) の順。
Codexの入口は [AGENTS.md](AGENTS.md)。分類構築の `research/classification/PLAN.json` は完了記録であり、現在の監査キューではない。

## 検証資料

- 主分類の論理入口：`research/classification/GLOBAL_EXACT147_COMPONENT_CLASSIFICATION_FINAL_PROOF.md`
- 従来の再現入口：`research/classification/REPRODUCE_EXACT147.md`
- 今回の証明ノート：`research/publication/verification_note_ja.tex`
- 独立検算：`code/verify_final_bridge.py`、`reports/final_bridge.json`

```sh
python code/verify_final_bridge.py --output work/publication-review/final_bridge.json
python -m unittest discover -s tests -p test_final_bridge.py
```

必要なのはPythonとSymPy。既設環境を再利用する。全legacy suiteやsetupの反復は不要。
`python scripts/classification.py check` は従来の計画・参照検査であり数学定理の検証ではない。

全PGL5軌道、Jacobiスキームの全非被約厚み、全交叉・局所環の決定は主分類に含まれない。
過去のREADME・STATE・原稿案にある下界や未完フラグを現在の結論と混同しない。
