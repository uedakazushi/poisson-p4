# 列17第四 pair 座標の一回限定実計算

2026-11-05。SOURCE_GATE 9e25c30c と literal-input V2 ab430bf4 を固定する。
worker は20 atom literal payloadと四つの既存κ多項式だけから49 sourceを評価し、
κ17=κ dot v17を計算する。元checkpoint/Q2/pair JSONは読み込まない。

係数は共通正整数分母と整数係数で保持し、指数 keyを出力項ごとには持たない。
各polyのbox/全次数による整数格子を暗黙添字にし、(a,e)行のoffsetを別途数える。
原 sparse inputだけはpacked keyとinteger係数を二slotとして計数する。
全denominator、row offset、poly header、49 nodeの数値metadata、整数temporary、
publication/readback phaseを含め75410 ownedを越えないことを実装で検査する。
548bit temporary、60秒/2GB、入力と全新成果物合計10.15MBを保持する。

各source atomのa shiftとQQ constant scaleはliteral sourceと項ごとに照合する。
全加減算でbackadd、全積で逆順全convolution再生を実行する。
source四値とκ17を別ファイルにstreamし、worker終了後、一ファイルずつQQ既約化・
再読込する。全四積は順序付きに足し、符号済κを二重に符号反転しない。

これはx4=κ17/Tだけの計算である。他三pair座標と十角block補正は保持し、
κ17の零/非零を最終一対八、全57、補集合の解決にしない。
