# Codexで研究を継続する

2026-11-03 v63研究追記：公開済みCloudを再利用中。再セットアップは不要。
M2 1.24.11の必要な既存計算は再現済み。QQのAlgorithm=>LinearAlgebraには
配列順序エラーがありdefault GBを用いる。新しい有限模型・原稿の再現入口は
[v63引継ぎ](../research/CLOUD_HANDOFF_20261103_v63.md)を参照。

2026-11-03の公開済みCloud再開では、既設の
`source /workspace/cloud-setup/activate.sh` と `.venv` を使う。
Macaulay2は1.24.11（移行元1.26.05）、FLINTは3.1で、全135式・有限行加群・QQ/GF拡大体・
原稿ビルドに必要な計算を再現済みである。成功したセットアップを再実行しない。
GroebnerBasisの列数には `numColumns gens G` を使う。
1.24.11のQQ加群で `gb(...,Algorithm=>LinearAlgebra)` は小例でもSIGSEGVを起こしたため、
正式証明書は既定のGBアルゴリズムを使う。QQ上のF4指定は既定法へ戻る場合があり、
F4で実行したと無条件に記さない。有理数の有限体化は分子と分母を別々に変換する。
再開・文献接続の実記録は `reports/CLOUD_RESUME_20261103.md` と
`reports/CLOUD_CONNECTIVITY_20261103.md` にある。

2026-10-30 / v25-setup。リポジトリの設定と、利用者のCodex側の実行環境を区別する。
このリポジトリは研究方針・原資料・依存・検証コマンドを用意する。
AGENTS.mdは実行中のエージェントへの指示であり、単独で常駐プロセスや定期実行を作らない。

## 1. Codexに渡す最初の指示

このリポジトリを選択して、次のように依頼できる。

> このリポジトリの研究を継続してください。最初にRESEARCH_POLICY.mdを読み、
> research/STATE.mdとresearch/NEXT_TASKS.mdから現在地を確認してください。
> 完全分類に向けた大域的な尽くしを優先し、意味のある成果の区切りで記録・検証・保存を
> 行い、続行可能なら次の課題へ進んでください。

短く「このリポジトリの研究を継続してください」でも、ルートのAGENTS.mdから同じ
方針へ到達する構成になっている。既存チャットの内容やLibraryの旧ファイルを渡し直す
必要はない。新しい優先順位があれば、その指示を追加する。

## 2. ローカル／デスクトップ／CLI

必要なものはGit、Python 3.11以上（検証時は3.12.14）、venvである。
macOS・Linux等、相対symlinkを利用できる環境を想定する。

```sh
git clone https://github.com/uedakazushi/poisson-p4.git
cd poisson-p4
bash scripts/setup.sh
.venv/bin/python scripts/verify.py --suite all
```

Pythonが別名の場合は、例えば `PYTHON=python3.12 bash scripts/setup.sh` とする。
2026-10-31はPython3.14.10でも全基準検証を再現した。この環境の既定Python3.15では
固定NumPy2.3.5を導入できなかったため、`PYTHON=python3.14`を明示した。
依存を勝手に更新せず、既存の不完全なvenvは確認してから別名へ退避する。
setupはリポジトリ内の`.venv`だけへ固定依存を入れ、版別原資料を`legacy/`へ展開する。
sudo、OSの変更、利用者のCodexの全体設定の変更は行わない。
以後はこのフォルダをCodexの主作業フォルダとして開く。
CLIが導入済みならリポジトリのルートで`codex`を起動する。

ネットワークが必要なのは初回の依存導入、文献の調査、GitHubへのfetch/push等である。
固定済みの計算証明書の検算は依存導入後オフラインで実行できる。
GitHubへの書込み認証は実行環境が提供する接続を使い、リポジトリへ認証情報を置かない。

## 3. Codex Cloud

既にこのリポジトリを使える公開済み環境がある場合は、それを選択する。
新しく作る場合は、公式の環境作成画面で `uedakazushi/poisson-p4` を選び、次を準備する。

- Python 3.11以上、Git、venv。
- Install scriptとしてリポジトリのルートで `bash scripts/setup.sh`。
- 起動確認として `.venv/bin/python scripts/doctor.py --strict`。
- 最初の検証として `.venv/bin/python scripts/verify.py --suite all`。
- 原稿を編集する環境には、日本語LuaLaTeX一式も導入する。
- 文献調査とpushに必要なネットワーク・接続は、その環境の設定に従って用意する。

2026-10-30に確認した公式手順は、`Work in > Cloud > Select environment > Create
environment` からリポジトリを選び、セットアップ確認後に環境をPublishする流れである。
`Settings > Codex Cloud > Environments` からも作成できる。古いCloud画面では
setup script欄等の名称が異なることがある。

リポジトリへのcommitだけでは、利用者のCodex Cloud環境の作成・Publishは実行されない。
このセットアップ作業の検証結果は `reports/SETUP_20261030_v25.md` に記録する。

一次資料（2026-10-30確認）：
[Cloud environments](https://learn.chatgpt.com/docs/environments/cloud-environments)。

## 4. 継続目標と自動実行

対応するCodex環境でGoalsを使う場合、いきなり不明確な完了条件だけを設定せず、
完全分類を上位目標とした上で、検証可能な研究単位を設定する。
最初の単位の例は「S61の定義式・復元写像・成分対応の根拠を現行入口へ統合し、
既知族の座標表示と再現手順を確立する」である。詳細な完了条件はNEXT_TASKSを使う。

対応するCLIでは `/goal` で目標を設定できる。利用できない場合も、AGENTS.mdの継続
方針を通常タスク内で適用し、終了時にSTATEとNEXT_TASKSを更新すれば再開できる。
Goalsはスレッド単位の継続目標であり、無制限の稼働や数学的な完了の保証ではない。
予算・時間・接続の制約に達したら、その事実を記録してcheckpointを作る。

定期実行の時刻や頻度はこのリポジトリでは設定していない。利用者が開始するCodexの
研究タスク内で継続するための指示を用意してある。必要なら別途、明確な実行条件と
資源上限を持つスケジューラへ接続できる。GitHub Actionsは使わない。

一次資料（2026-10-30確認）：
[AGENTS.md](https://learn.chatgpt.com/docs/agent-configuration/agents-md)、
[Using Goals in Codex](https://developers.openai.com/cookbook/examples/codex/using_goals_in_codex)、
[Non-interactive mode](https://learn.chatgpt.com/docs/non-interactive-mode)。

## 5. 検証と出力

| コマンド | 確認する範囲 |
|---|---|
| `python3 scripts/bootstrap_sources.py --check` | 23版のZIPと展開済み1885ファイルのSHA-256、版間リンク |
| `.venv/bin/python scripts/verify.py --suite quick` | v24の全1298項目。quickは一部の抜き取りではない |
| `.venv/bin/python scripts/verify.py --suite rank2` | v14の既存10検証群 |
| `.venv/bin/python scripts/verify.py --suite integrity` | 原資料ハッシュと、旧v24/v23 manifestの再現検証 |
| `.venv/bin/python scripts/verify.py --suite all` | 上のintegrity、quick、rank2 |

出力は既定で新しい `reports/<UTC日時>-<suite>/` に置かれる。
`--output reports/任意の新しい名前` でも指定できる。`--output` はディレクトリである。
元証明書や既存レポートを上書きしない。検証に失敗すれば非零の終了値を返す。
v14は`--jobs`で並列数を指定できる。

PASSが意味するのは、実行した恒等式・階数・証明書等の確認であり、全幾何学的な
場合分けの独立査読や、新規性の証明ではない。v3〜v8の全大規模証明書やv14〜v22の
全依存鎖を、この短いsuiteがすべて検証するわけではない。

既知の固定入力を毎回全部実行する必要はない。入力・実装・環境が変わった範囲と、
新しい主張の成立に必要な範囲を検証する。

## 6. 日本語LuaLaTeX

現在の評価ノートPDFとソースは `research/notes/` にある。
再ビルドは `bash scripts/build_notes.sh`、出力は `build/notes/`。
LuaLaTeX、ltjsarticle、luaotfload、HaranoAji、通常の数式・表・ハイパーリンク用packageが
必要である。OSのTeX環境はsetup.shでは変更しない。
Debian系なら管理された環境準備で `texlive-luatex`、`texlive-lang-japanese`、
`texlive-latex-extra`、`texlive-fonts-recommended` 等を用意する方法がある。
macOSでは日本語LuaLaTeXを含むTeX環境を使用する。特定OSのインストールを自動実行しない。

公表用原稿の本文は各版の原資料に保存してある。独立したA/B/Cの完成原稿を新たに
執筆した状態ではない。証明・PDF・ソースの入口は `papers/README.md`。

## 7. checkpointとbundle

まず状態文書とノートを更新し、diffを確認して意図したファイルだけをstageする。

```sh
git status --short
git diff --check
git add <今回保存するファイル>
.venv/bin/python scripts/checkpoint.py --message "research: 証明したことと残した課題" --push
```

checkpoint.pyはstage済みの変更だけをcommitする。main以外では停止し、他セッションの
変更を取り込んだり消したりしない。pushは通常のfast-forwardだけで、forceは使わない。
遠隔で変更が進んでいたら、fetchして差分を確認し、必要な統合後に再実行する。

```sh
.venv/bin/python scripts/make_bundle.py --version v25-setup
```

bundleは`bundles/poisson_p4_continuation_日付_版.zip`となる。原資料アーカイブを含み、
依存パッケージや再展開したlegacyは含まない。旧bundleを上書きしない。
このコマンドはGitHub checkoutで実行する。展開した配布bundleだけでもsetup・検証・
研究はできるが、版管理とpushを再開する際は元リポジトリへ適切に接続・統合する。
Libraryへの保存は別の認証済み機能が必要であり、このスクリプト自体はアップロードしない。
保存許可はRESEARCH_POLICYに記載済みなので、機能がある場合は再確認を求めず保存する。
