# Poisson structures on projective four-space

複素射影4空間のPoisson構造の完全分類を目指す研究リポジトリ。
**2026-10-31 · v27-s61：既知三成分の一次部分への優越性、トーリック復元、MPS40の個別照合。**
完全分類は未完了で、既知成分数の下界57は維持します。

T01・T39=q(5,1)・T40=q(5,2)から16次元のprimitive一次部分への射影が
すべて優越的で、一般に有限であることを確認しました。十分一般の一次部分には
少なくとも三つの相異なる孤立被約な延長があります。未知成分や追加の延長は排除していません。
トーリック開部分は16変数coframeのS5商へ還元し、MPS40の全図式を一次資料と照合しました。
[現在地](research/STATE.md)、[楕円型の証明](research/notes/ELLIPTIC_WEIGHT1_20261031.md)、
[再現コード](code/S61_README.md)、[日本語研究草稿](papers/s61_first_jet_reduction_ja.tex)を参照。

## 再開の入口

1. [研究継続方針](RESEARCH_POLICY.md) を最初に読む。
2. [現在地](research/STATE.md) と [次の課題](research/NEXT_TASKS.md) を読む。
3. [Codexでの起動と環境設定](docs/CODEX_SETUP.md) に従ってsetupと検証を行う。

Codexでこのリポジトリを開いたら、**「このリポジトリの研究を継続してください。」**
で再開できるように、[AGENTS.md](AGENTS.md) に研究・検証・保存の指示を置いています。
通常の指示ファイル自体はバックグラウンド実行を開始しません。継続目標機能が利用
できる環境での使い方もsetup文書に記載しています。

## 研究の現在地

- 第一段階の目標は、射影化されたJacobiスキームの**被約既約成分の尽くしと一般元の幾何**。
  PGL5商を取る前の空間を扱います。全特殊軌道、交叉、冪零構造は別の達成基準です。
- v14には階数<=2部分の8成分の尽くしがあり、そのうち6個が一般階数2の全空間成分です。
  計算証明書は再現済みですが、すべての幾何学的議論の独立査読完了ではありません。
- 全成分数は **少なくとも57 = 40 MPS + 6 rank2 + 11追加rank4**。正確な総数は未確定です。
- 大域研究の主軸は **S61の61変数・134式**。16個の重み1方向を含む全模型を扱います。
- v22の正規退化因子定理には依存鎖の監査が残ります。v24の端点Qの局所成果は
  全体の尽くしとは区別します。

詳しい確実性・依存・証拠は [主張台帳](research/CLAIMS.json)、
[成分台帳](research/COMPONENTS.json)、[評価ノート](research/notes/README.md) にあります。
MPS40を含む既知57成分は個別台帳へ展開済みです。51階数四成分すべてのS61座標・イデアルは未構築です。

## 最短の再現手順

Python 3.11以上とvenvを用意し、リポジトリのルートで実行します。

```sh
bash scripts/setup.sh
.venv/bin/python scripts/doctor.py
.venv/bin/python scripts/verify.py --suite quick
.venv/bin/python scripts/verify.py --suite rank2
.venv/bin/python scripts/verify.py --suite integrity
```

Python依存は固定しています。日本語PDFの再ビルドにはLuaLaTeX、luatexja、HaranoAji等が
別途必要です。既存PDFは同梱されています。詳細は [setup文書](docs/CODEX_SETUP.md)。

## ファイル構成

| 場所 | 内容 |
|---|---|
| `RESEARCH_POLICY.md` | 恒常的な研究原則、目標、現行計画、保存許可 |
| `research/` | 現在地、主張・成分台帳、課題、文献、評価ノート |
| `papers/` | 公表用原稿の入口、元原稿の証明所在地、再編計画 |
| `sources/` | v02〜v24の全原資料を版ごとに保存したZIPとSHA-256台帳 |
| `legacy/` | setupが生成する原資料の展開先。Git管理外・原本変更禁止 |
| `audits/2026-10-30/` | 前回評価時の監査範囲と再検証結果 |
| `reports/` | この環境で実行した検証結果とsetupの記録 |
| `scripts/` | 展開、環境確認、検証、PDFビルド、bundle、checkpoint |
| `code/` | 今後の計算コード。原本の修正はここに新実装として作る |

原本の全ファイル内容は保存しています。ZIPコンテナは版ごとに再梱包しており、
元ZIP自体と同一のバイト列ではありません。[保存形式と検証](sources/README.md) を参照。
原資料には研究草稿・未完成の議論が含まれます。現在の評価は現行台帳を優先してください。

## 保存と継続

作業の区切りと、長い作業ではおおむね30分ごとに、意味のある変更をmainへ保存・push
する方針です。GitHub Actionsは使用しません。保存のために毎回ユーザーへ確認しません。
`scripts/make_bundle.py` は再開に必要な資料をまとめた配布ZIPを生成します。
Library保存は、その機能が実行環境にある場合に許可済みの操作として行い、実際の
保存結果を記録します。リポジトリの設定だけでLibrary接続が追加されるわけではありません。

既存の強い結果の論文化を進めながら、未処理領域の尽くしを優先します。
このリポジトリにあることをもって査読済み・完全分類済みとは扱いません。
