# Bianchi 残余4方向と合計37方向：独立監査

2026-11-04。対象は frozen `S61_NONSQUARE_BIANCHI_RESIDUAL_FOUR_20261104_v63`。
既受理の親15本の actual319 関係、その10710-node DAG、11個の quartic pivot
の独立レビューを使う。producer、pivot 探索、full175 の展開、239列行列、
C122、GB、新しい点や体の計算は実行しない。

親の15×15行変換の最後4行は、元 block の行0,4,13,14の単位ベクトルで
あることを確認する。従ってこの4本は既に証明した actual319 関係そのもの。
保存された全319 source definitions と合わせ、112 cubic 列と33 quartic 列を
独立に特定する。親の15係数から選んだ root と今回保存された root を全て比較する。

author の recursive cache や Oracle evaluator は呼び出さない。必要な scalar
node の依存閉包を先に取り、添字の昇順に **別の bottom-up QQ(a,e,d) 演算**で
評価する。全132 quartic 成分が零、全448 cubic 成分が保存値に一致し、その
119成分が非零であることを確認する。親 DAG は追加・変更しない。

保存済み4回の行操作を、探索せずに独立に追跡する。pivot は

    189a²/(2F), −756a/F, −36450a/(FB), −50000a/(FB)

であり、分子・分母を a,F,B,C の冪と非零 QQ 定数へ戻して検査する。
swap、各 elimination、4×4変換、最終4×112全係数、および選択列の I4 を
比較する。pivot の319列は176,177,237,181。対応する source は
`G18*v2,G18*v3,G86*v0,G20*v0`（各 G は既受理 row-unit 正規化済み）である。
さらに `scope-receipt` の明示的4×15係数が、実際にこの行変換と親4行の
合成であることを検査する。未展開の quadratic cofactors もこの合成により
一意に定まるので、それらを展開する必要はない。

新11本と新4本の選択15列は `(I11,*;0,I4)` となる。したがって **全15個の
direct weight4 pieces は aFBC≠0 全域で点別に独立**である。旧22方向は
quadratic source columns だけに支持を持つため、既受理の独立性開
`FULLU=a*d*F*B*C*J≠0` 上では三段の三角 block によって合計37本が独立。
80方向 gate には **43本不足**する。旧22の独立性を aFBC 全域へ広げない。

元 report の `old22_plus_direct15_rank_on_FULLU=null` は cubic 全零なら、という
未実行の代替分岐の field である。今回の cubic rank4 と明示的合成による37の
結論を否定しない。`scope-receipt` と本レビューはこの正の分岐を明示する。

全 weighted Bianchi multiples や全 syzygies を尽くしたとはしない。
full80、actualC122 rank64、全点239下界、full quartic 全射、affine exhaustion
は未証明のままである。専用 report は arithmetic、入力 hash、時間とこの範囲を記録する。
