# Poisson structures on projective four-space

**2026-11-07 · v64-classification-first**
**主目標：Poisson構造の分類を主定理とする論文。派生論文はsecondary。**

最初に [研究方針](RESEARCH_POLICY.md)、次に
[現在地](research/classification/CURRENT.md) と [実行計画](research/classification/PLAN.md) を読む。
Codexの入口は [AGENTS.md](AGENTS.md)。[次の課題](research/NEXT_TASKS.md) は現行キューだけを指す。

## 分類の意味と現在地

対象は複素射影四空間上の非零Poisson双ベクトルをスカラー倍で射影化したJacobiスキームで、
PGL5商を取る前の空間である。第一主定理は、その**被約既約成分の尽くしと一般元の幾何**。
特殊点を含む全PGL5軌道分類や全冪零構造とは区別する。

現行台帳の成分下界は **92（一般階数4は86）**。上界・一般一次ジェットの五点尽くし・
全成分の尽くしは未証明。全92成分を個別同定し、既約二次因子＋既約三次因子などの全支持を被覆した。

現在の最短候補は、一般ファイバーの未知**被約支持**を消し、得られた開集合の補集合を
実際の方程式で尽くすこと。B18の正確なnil指数などは、必要性を示さない限り休止する。

```sh
python3 scripts/classification.py check
python3 scripts/classification.py next
python3 -m unittest discover -s tests -p test_classification.py
```

初回環境構築だけは `docs/CODEX_SETUP.md` を参照する。既設環境でsetupを繰り返さない。
数学検証は変更した依存範囲に限定する。上のcheckは計画と参照の検査で、数学の証明検証ではない。

## 必要な資料

- [証明の経路と停止基準](research/classification/PROOF_ROUTE.md)
- [証拠・コードの入口](research/classification/SOURCES.json)
- [主論文の日本語LuaLaTeX原稿](papers/poisson_p4_classification_ja.tex)（分類主定理は未証明）
- 既存 `research/CLAIMS.json` / `research/COMPONENTS.json`：数学的主張・成分台帳。

旧 `research/STATE.md`・過去handoff・既存五原稿は証拠資料として残るが、旧優先順位や
過去の「未決」記載を現行方針として読み直さない。原証明・旧証明書を保持し、新たな全源と独立検証を現行原稿へ統合する。
