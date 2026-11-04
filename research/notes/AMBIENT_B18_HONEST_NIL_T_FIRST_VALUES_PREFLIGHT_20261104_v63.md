# honest T_nil の小さい定数・一次係数：既存 source の完全性と計算上限

2026-11-04。本packetは preflight のみである。m0、X0/X1、rank、nil生成元や
T_nil長さを計算したとは主張しない。凍結済み128×972 honest operatorを変更しない。

係数環は actual CA-F9 の S=k[[zA0,...,zA6,rA]]、ell=5(sA−rA)=0。
切断変数rAとDVR parameter ellを区別する。全九 proper 座標の全八方向一次係数は
NATIVE128_FIRST phase4に、実u5=G8[:5]とN12の全一次係数は同phase5に存在する。
第三gamma sectorでのみ必要なrA二乗係数は、独立監査済み native S2/B11 phase3に
保存されている。B11が捨てるmonomialは今回のz一次・r一次・r二乗に寄与しない。
不足するκ方向は無い。

元のproper/source座標から実dR9、D8×9を評価する。C36はdual9と自由四trace基底の積であり、
nuA/nuBの一次値を保持する。tau平方を零に固定したまま係数だけ微分してはいけない。
最初の五pivotも、交叉上の有理rootではなく実u5を使う。実sourceのN12を単元逆行列で
反転し、A=N12^-1 injection Dを作る。旧B640の全一次値をRes4(A^t)で変換する。
最後140列はradical sourceの全35 syzygyを同じC36で制限する。旧P8はそのまま保持する。

旧640列のキャッシュはB部155 entries（定数4個、一次245個）、P部105 entries
（定数零、一次160個）という疎なものだった。native35の非零entryは70個。
従って新X0/X1出力は20,356個以下の非零QQ係数で構成できる。
定数消去、84×780のW0上限、単元逆行列、固定minorとlive bufferまで加えた保守的上限は
128,812 QQ slots、全出力32MiB、60秒・2GiBである。巨大823680積は列挙しない。

将来の値workerでrank(m0)=24を確認できれば、24本のliteral columns Eを選び
n=E(mE)^-1とできる。Vgamma=d1mからVgamma n=d1、mn=Iが正確に従う。
旧nのr三乗正則商を数値評価する必要はない。

r=rank X0のpivotを消去し、X0 c=0という定数組合せについてのみ一次値を使う。
このとき(mXc)_1=m0 X1 cであり、未計算のm1は消える。定数Koszul liftは一意で、
pair(i,j)の係数はcycleのj成分におけるzi係数から読み、全24×8等式で照合できる。
第三sectorの最後のregular generatorはrA二乗だが、pair(i,7)はzi係数により同様に一意。
従ってW0を得ればmu_S(T_nil)=104−r−rank W0となる。
これは最小生成数の式であり、総長、nilidealの全生成元、nil三乗や正確な指数の答えではない。

本preflightはsourceファイルの存在・schema・全方向の保存と疎supportだけを照合した。
新しいsource値、gauge/IFT、jet、Gaussian/rank workerは実行していない。
