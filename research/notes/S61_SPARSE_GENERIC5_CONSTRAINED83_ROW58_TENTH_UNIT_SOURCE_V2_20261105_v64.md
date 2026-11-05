# 第10段固定単元：全167×14と全RHSを保持するSOURCE

ROOT独立value review fff391b8による実pivot

p10=F10/(353894400000 H t0 F9)

を使う。F10は指定された317項のprimitive整数多項式であり、既約性やF9との互いに素は用いない。既存のH,t0,F9等のchartにD(F10)を追加する。F9を分母から消さず、F10=0を含むすべての旧・新補集合は未処理として残す。

現在の全168×15からrow58/col56を消去し、全167×14を得る。topは101×14、selected squareは14×14、残り87本のtop RHSと66本の零制約RHSをすべて保持する。S10[i,j]=S9[i,j]−lambda[i]S9[58,j]、lambda[i]=S9[i,56]/p10、h10[i]=h9[i]−lambda[i]h9[58]。全pivot行と取り除いたh9[58]を復元用に保持する。selected14本のlambdaは受理済みcol56 binary参照、他87本は完全なsource callbackで未評価、66本は係数が恒等的0でもRHS整合条件を落とさない。row58全体を旧raw行やrow35と同一視しない。

既知U=Q^4 R^2 E、D9=U F9に対する共通表現は

D10=U^2 F9 F10,
N10[i,j]=U F10 N9[i,j]−353894400000 H t0 N9[i,56]N9[58,j]

である。これは型付きの式で、他14列の係数を計算したことにはならない。新しいパラメータ因子の除去・GCD・数値点・新しい行選択を行わない。P290は既存4列に掛けた一回だけを参照する。

全168×15は単元block[p10]と全167×14の直和へ移り、cokernelは同型である。一方、零制約行66本を落としたtopだけなら、全cokerはtop cokerとK5^66の直和である。selected14平方を全cokerと同一視しない。x56の復元式、その前の9段、10個の制約単元、49段の古い操作も残る。

新しい行列係数workerは未実行で、全体の逆行列・exhaustion・chart補集合の被覆を主張しない。次の係数計算には、固定した一列に対する支持・全途中高さ・owned・RAW・実装custodyの新SOURCE gateが別途必要である。

V1のSOURCE作成では数値owned式の手書き比較値42056に誤記があり、式の実値41200とのassertで停止した。実際の行列係数計算ではない。V1 code・型付きoperator・failure bytesをそのまま保存し、V2は比較定数だけを41200へ訂正する。予算を引き上げず、数学式・入力・row/column・chartを変えていない。
