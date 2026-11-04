# Bianchi quartic の11単元方向：独立監査

2026-11-04。対象は frozen
`S61_NONSQUARE_BIANCHI_QUARTIC_UNIT_PIVOTS_20261104_v63` だけである。
専用 verifier は producer を呼び出さず、既に独立受理された full175
Schouten operator、oracle、frame 共変性、旧 cubic row-basis、旧22方向を
hash で固定して使う。239列の行列、C122、元175恒等式の再展開、GB、新しい点、
新しい pivot 探索は行わない。

## 有限算術の確認範囲

元 q/t 斉次135式を local chart に戻し、全135式の J/B atoms と照合する。
元 oracle DAG の prefix、全 core atom の既知単元分母、全 node の依存順・
重みを検査する。19/23 nuisance coordinates は E,Z2,Z3,Z4 のすべてで零。
従ってこれらを weight4 block から除く操作は元 operator の削除ではない。
上流の full135（G6 を含む）恒等式と tG6 は保持されている。

保存済み15×33 block に対する11回の行操作を、QQ(a,e,d) で独立に追跡する。
pivot は順に
`−1,−1,−1,−1/10,−1/5,−a,−2/5,−1/2,−1/2,−1,−2`。
各 pivot は aFBC 開の既知単元だけを使い、残り4行は厳密に零、選択11列は
I11 である。保存した elementary transform と transformed block も全係数で
比較する。

full63 archive の full member 全体を SHA で照合し、実際の header
`63 91 57` を読む。writer と parser により第3整数57が target metadata で、
63×91 のサイズだけが検査されることを確認する。57を rank 証拠に使わない。
最初の実行失敗と一回の resume の記録も保存されたまま照合する。

補足5 cubic 行の native anchor、row normalization、実際の F_g(Dv) anchor を
独立に比較する。whole-row identity の根拠は、受理済みの全63行の native
source 再構成と `Tn=T34D, Un=ud U` の共変性である。この証明は元行ごとに
適用され、anchors はその既証 scalar を定めるだけである。新しく21行全体を
展開比較した、または一つの anchor だけから行恒等式を推測したとはしない。

全15 piece の multiplier に、weight1 の E、weight2 の E+Z2 と y=Dv を代入し、
6 quadratic units と全21 genuine cubic units を使う。選択34式は exact graph
上で零なので除ける。正重みにより Z3/Z4 以降はこの coefficient 抽出には
寄与しない。各354係数 DAG を復元し、35個の既受理単元関係を順に引く。
各削除 step と、残る全319係数を保存値と比較する。35関係には quartic
成分がない。最後に elementary transform と共通 D_H² の除去を行い、全11本の
319係数、およびその選択 quartic 成分が I11 であることを検査する。
これは有限 DAG による厳密な恒等式輸送であり、全239列の展開検算ではない。

## 結論と限界

上記の有限検査が通れば、新しい11本は **aFBC≠0 全域で点別に独立な actual319
生成行の関係**である。旧22本は quartic 成分が零なので、既受理の独立性開
`FULLU=a*d*F*B*C*J≠0` に制限すれば合計33本が独立になる。
旧22本を aFBC 全域へ拡張しない。80本 gate には **47本不足** する。

full80、C122 rank64、239点別下界、full quartic 全射、affine exhaustion は
未証明のままである。残余4本の quartic block が零であることも、それらの
他成分での独立性を否定しない。実行済み検査の成否・時間・全 input hashes は
専用 `..._REVIEW/report.json` に保存する。
