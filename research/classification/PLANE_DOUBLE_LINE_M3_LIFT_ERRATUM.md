# m三candidate手liftのaffine quotient符号訂正

2026-11-07。主義務G。旧手noteと初回finite reportのbytesは保存し、ここで正しいsourceを指定する。

`PLANE_DOUBLE_LINE_M3_CANDIDATE_LIFTS.md`及びfinite-connections note末尾の第四candidate不整合予測には誤りがある。
affine chart b=1でのnormal bracketの正しいquotient式は

`{z/b,w/b}=Pzw/b²+w Pbz/b³−z Pbw/b³`。

旧手展開では後二項の符号を逆にした。
global quadratic Kを自由な三未知数のまま実際のEuler-Pfaffianから再構成した
`CLASSIFICATION_PLANE_DOUBLE_LINE_M3_SQUAREFREE_REPLAY_20261107.json`は、この手予測を仮定していない。
全四candidateをbatchに保持し、正しい全30係数matrixからunique Kを求めた。

正しい全値は

| candidate connection | global raw K |
|---|---|
| C=a²、g=2 | 2z² |
| C=a²、g=3 | 0 |
| C=a²、g=4 | −2z² |
| C=b²、T=a²b² | 0 |

であり、第四もglobal Kを持つ。
前者三つのaffine normal Pfaffianは`−3t²[K−(6−2g)z²]`、第四は`−3t²K`となる。
旧式K=(6−g)z²及びK=6(a/b)z²は使用しない。
raw mixed columns、全root-cover/Riccati/collision分類及び全四candidateを保持する数学は変わらない。

第四のSF空性はK不整合からではなく、**完全weight一kernel上のactual D3像がa²z³-lineだけ**であることから従う。
g=3及び4も全D3像がrepeated-root cubicに限られる。
g=2の像にはsquarefree cubicがあるので、その枝は別のB≤9 source次元証明でG非優越へ戻す。

本erratumは誤った手予測の撤回であり、旧reportのfalse coarse上界を成功へ書き換えない。
全19sourceの別作者reviewは、保存されたactual K=2/0/−2/0を用いて再構成する。
Cloud resetで元の138f21ee… receiptは失われたため、fresh bounded replayを新identityで保存した。
元の実測値だけをmemory内で戻したhash比較は138f21ee…と一致し、他の全bytesの再現を確認した。
fresh reportには新しい実測値を保持し、旧runtimeを上書きしていない。
