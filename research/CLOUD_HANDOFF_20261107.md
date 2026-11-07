# Cloud研究引継ぎ — 現行分類を継続する

2026-11-07。公開済みCloud、repository `/workspace/poisson-p4`、branch `main`。
これは移行先で継続中の研究であり、旧ローカル移行の一時停止文を研究停止命令として扱わない。

## 開始と環境

最初のCloud開始SHAは `4c3f623fb698f9adc6077015b25e2325b11072af`、
11月7日研究再開のbaselineは `a00a541b835b2bff18bafcba545ecbd4bae7175b`。
公開済み環境と `/workspace/cloud-setup/activate.sh`、既設 `.venv` を再利用した。
Macaulay2はCloud1.24.11、元環境1.26.05。使用するQQ行列kernel/既定GBの互換小例を確認済み。
新しい主計算はSymPy1.14と独立Fraction疎演算であり、1.26専用APIを仮定していない。
文献の接続・一次scopeは各専用receiptへ保持し、Saito原独語全文未取得を取得済みとしない。
現在の全八平面quintic尽くしは、局所quasihomogeneity外部判定に必須依存しない直接証明も独立確認した。

Singular実体は `work/nonsquare-v62/singular-runtime/prefix/usr/bin/Singular`。
旧 `work/native-cas` が無いことをソフト不在と誤認しない。セットアップ全体を理由なく再実行しない。
日本語LuaLaTeX、PDF可視確認、通常GitHub fetch/pushは動作済み。network追加権限は既設proxy/TLSを保持する。
`legacy/` はtracked原ZIPの展開入力。新しいmainへforce-addせず、必要memberのexact bytesを原ZIPから検証する。

## 数学と保存の境界

現行個別下界は111（rank四105、rank二6）、主張台帳931。
旧C45と新七平面quinticの全周囲23成分性・相違、全explicit source/cusp、全必要curve分類を別々に独立確認。
**本質三変数の既約退化五次は八つの23成分による全点被覆まで完了。全G/E・全上界は未証明。**
最新の義務・証明・未閉runは `classification/CURRENT.md`、選択順は `classification/PLAN.json`。
全m二rank零plane、scalar/zero symbol、common-image、全smoothquadric等は受理済み。
残るdependent高mの全例外、rank二cone/scroll/plane、canonical index二、Eのstable二源・markedquartic等を続ける。
D5+tangentはsource23/D101/二次余核零で止めた。これを新成分又は高次計算成功に数えない。

保存開始前の公開baseは `6ab28a327b64543ba0a7bf8cd302bea9def517e6`（104/916、187頁）。
このbundleは111/931、主稿201頁、最終2pass/主要頁可視確認と依存完全性を検証した保存単位。実際のGitHub保存状態は、
`git log`、`git status`、origin/mainと当該build receiptを確認する。
mainへの通常commit/pushはユーザーが明示許可済み。root統合担当だけがGit/大域台帳を編集する。
Libraryはupload receipt無し。生成、検証、GitHub保存、Library保存を区別する。

## 再開操作

`RESEARCH_POLICY.md` → `classification/CURRENT.md` → `classification/PLAN.md` → `NEXT_TASKS.md` を読む。
`.venv/bin/python scripts/classification.py check` と `next` を使い、非自明runは有限job契約をgate後に実行する。
checkerは数学を検証しない。既に閉じた結果を欠陥無しに再開しない。
二回非閉鎖なら支持幾何/actual family等へ変更し、同じGB/jet/minor列を延長しない。
既存FunctionDef/producer reports/certificatesを保ち、compatible拡張前の実行bytesと契約hashをsnapshotへ保存する。
節目ごとに主稿・CURRENT/PLAN/未完了課題を更新して通常pushし、その後も有効なprimary数学を続ける。
現時点でユーザーの追加許可が必要な環境阻害は無い。未完了の理由は記載された数学的義務である。
