# Cloud接続中断からの復旧と非公開部分の再構成

2026-11-07。研究は継続する。分類全体のG/E・上界は未証明である。

## 確認した保存境界

中断直前、main の `8b46f29ddbb9f539abc7a026e492ee25bf8ad699` を通常pushし成功した。
個別124成分（rank四118・rank二6）、967主張、日本語主稿225頁はGitHubに保存済み。
PDF SHA256は `98ce53afff462fa85ec2be98809007c3ea9257843a1e8ad72b9b4ca756451a49`。

Cloudがunavailableとなり、exec-server接続不能が実際に発生した。
再接続後のworkspaceはHEAD `4c3f623fb698f9adc6077015b25e2325b11072af`、
research/classificationの保存済み内容も無い初期snapshotの状態だった。
原因の内部機構は未確認だが、単なる同一filesystemへの再接続ではないことを確認した。

rootは既存untrackedファイルを保護し、`git fetch origin main` と
`git merge --ff-only origin/main` で公開checkpoint8b46f29dを復元した。
復元後の124/967、PDF全hash、tracked差分無しを確認した。
接続復帰直後にworkerが作成していたcurl/erratum/partial-boundの三untracked noteも保持した。

## 中断後の再現規則

公開checkpoint後の未commit部分は新workspaceから読めなかった。
手証明・コードは会話に残る元tool argumentsから再構成し、既知のhashと一致するものだけをexact復元と呼ぶ。
計算reportには実測時間・実行source・exact inputsを含むため、再計算結果を消失reportと同じidentityで作り直さない。
fresh replayは新report名へ保存し、旧hash・旧測定値はhistorical metadataとして明記する。
独立算術を未実施の結果をPASSに変更しない。

主な非公開部分は、m二reducedR squarefreeの追加純証明、m三nonredR全19source、
smooth cubic scrollのdirectrix/root-line還元、degree五canonical index二の追加構造、
T37有限wallのquartic-field/全曲線閉包証明である。
各担当が手証明・入力を先に再保存し、必要な有限算術だけを新identityで再現する。
全成分数はこれだけで増やさず、台帳の受理はrootが独立証明の完了後に行う。

m三の消失旧producer reportは SHA256
`138f21ee8fb4c8966c90e3451617a17c06779fc5a9058becc3b5b34727ed933b`、
旧source SHA256は `95981c680721a562462a81d1f97489e87c612e7a8ac1384fe9bcf9ef374bf0b2`。
その数学はproducerで完了していたが、独立19source算術は中断時未開始。
旧runtime20.624秒/RSS123,960KiBをfresh replayの実測として流用しない。

## 環境と継続

公開済み `/workspace/cloud-setup/activate.sh`、.venv SymPy1.14、既設Macaulay2を再使用する。
セットアップ全体を再実行しない。旧work内のSingular最適化wrapperは復帰snapshotに無かった。
これは旧任意tool実体の消失であり、保存済みQQ/Fraction証明の無効化ではない。
必要なtracked ZIP memberだけを元hashで再展開し、既設CAS又はQQ算術で必要な計算を行う。

中断は研究停止指示として扱わない。CURRENT/PLANのprimary義務を保ち、
復旧済み主稿を根拠なく作り直さず、未完了の研究・独立確認・通常commit/pushを継続する。
