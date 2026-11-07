# Cloud研究引継ぎ — 現行分類を継続する

2026-11-07。公開済みCloud、repository `/workspace/poisson-p4`、branch `main`。
これは移行先で継続中の研究であり、旧ローカル移行の一時停止文を研究停止命令として扱わない。

## 開始と環境

最初のCloud開始SHAは `4c3f623fb698f9adc6077015b25e2325b11072af`、
この継続区間の開始SHAは `37d56e2afb3fd467401d395e39f9c34f32bae7b3`。
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

現行個別下界124（rank四118、rank二6）、主張台帳967。
全essential三変数既約quinticは八23成分へ全被覆、μ六markedquarticは九新23成分を独立採録。
D5cusptangent全源は旧A0、3A2cusptangent全源は旧DPMへactual全source認識済み。
全rank零曲面支持のG非優越、全degree四/quadric rank二、plane m四全R等も完了。
**全G/E・全上界は未証明。** 残るrank二plane m二三、cubicScroll/cone、degree五canonicalindex二、
Eのstable二源、A5smoothhyperflex未知親、E6smoothflex未知親、他明示残支持を続ける。
全globalT39/T40のS61判定同値をconnectedSp4＋globalGa4切断で独立確認し、
既E6点を別Fraction再構成して両親の外と証明した。未知親/次元はそのまま残す。

直近の公開保存は `37d56e2afb3fd467401d395e39f9c34f32bae7b3`（111/931、201頁）。
本保存単位は124/967、主稿225頁。受理済み全G閉鎖と十三markedquartic新成分を本文へ統合し、
最終二回組版/全参照/主要頁可視及びpublic依存閉包を検証した。
実際のGitHub保存状態は `git log`、`git status`、origin/mainとbuild receiptで確認する。
最新の数学・残余は `classification/CURRENT.md` と `classification/PLAN.json` を使う。

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
