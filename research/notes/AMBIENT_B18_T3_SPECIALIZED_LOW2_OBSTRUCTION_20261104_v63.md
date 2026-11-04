# Actual CA平面へのspecialization後のLOW2：top12は全て二次まで零

2026-11-04。T3の非零性を低次negative certificateで調べる最初の試行。
同じphysical R0からactual CA zero7平面へのspecializationを行い、二変数(r,s)の
二次まででQ[3]の像とtop12を比較した。
specialization後の像に入らないなら元の像にも入らない、という負の向きだけを使う。
specialization後のmembershipや原点matrix零性から、元のmembershipを推論しない。

original sourceのjetは既保存のcombined105解を読み、二次成分を取り出した。
新しいIFT/Newton solveはしていない。
784 coordinate商と四つのg商の未計算微分を0には置かない。
各正則商を既知の原点値と、独立なorder≥1の形式的偏差に分けた。
g_CB=rA+sBのA,Bは既知のorder≥2を使い、g_CB tail自身はorder≥3なので二次では零。
cloneの未計算変数を残さず、上記商とg tailでDAGの依存を切ることを確認した。

各Q[3]列の二行を、(r,s)と形式的偏差の総order≤2の多項式として計算した。
偏差を任意の正order(r,s)-seriesに代入して得られる全jetを含むように、
各偏差monomialの可能な像を独立にvector-spaceへ追加した。
この手続きは実際の像を拡大するだけである。
Q[3]2034列の全許容多項式倍を二次まで含めた拡大spanのrankは3。

しかしtop12の二行は、拡大に用いた独立偏差を形式変数のまま保っても、全て二次まで零だった。
従ってこのLOW2試行は非零classを検出できない。
T3=0、actual membership成立、全Dのdepth20のいずれも意味しない。
同時に、同じspecializationによるnegative testの最初の可能orderは3以上だと分かった。

Source: `code/ambient_B18_T3_specialized_low2_obstruction_20261104_v63.py`。
`reports/AMBIENT_B18_T3_SPECIALIZED_LOW2_OBSTRUCTION_20261104_v63/`に、
全2034既知vectorと偏差の拡大vector、span basisと出所、top12の形式多項式・残差を保存。
一回0.723408秒、最大RSS142636KiB。新点、GB、高次jet、source再生成は無し。
自動的な次数増加はせず、このLOW2結果を独立に凍結した。
