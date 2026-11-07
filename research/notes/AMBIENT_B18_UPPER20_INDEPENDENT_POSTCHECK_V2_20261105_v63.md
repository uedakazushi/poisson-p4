# 固定20件の独立signed postcheck

ROOTの固定一回実行許可の後、既知入力だけで20件を独立に照合した。全20件のupper二成分とQの六成分は零で、元primaryの各成分および固定rank零・I20核証明と一致した。新しい物理・gauge・高次function・別case・rank探索は行っていない。

primaryのdirect N転置Sc積とは異なり、先に90個のcovector-L行を作る。通常68行に、未制限の形式的raw補助源の点成分を運ぶ22行を加える。16個のraw補助源はあくまで同じpartial-N拡張での形式値であり、8/32/80の非零S0を消していない。その後に完全な点核制限を取り、15個の実際の核源と10個の既証明SOURCE零を照合した。varying monicと全first-crossを保持する。

OLDはdense B21多項式を逆順に積和し、raw21件の形式的部分値を比較する。FIRSTはL0Se＋LeS0を個別源のまま再構成し、全1,742係数積と合計を照合した。native617の定数座標は既に補正済みなので一回だけ使う。raw617＋20raw563はOLD/FIRST/外側accepted92項だけへ適用し、raw563の部分upperを実際の零関係とみなしていない。

実測は147 bit、数値所有上界40,266、約0.642秒であった。事前上限は354/1063 bit、57,977/60,000数値slot、RAW 4,000,000 byte、外部wall 60秒、AS 2GBであり、変更しなかった。V1の制御だけを改善したV2は、現在のL累積値・kernel制限・FIRST部分・各case四部分およびpending演算を保存する。旧V1コードとcustodyは未実行のまま不変である。

既存の584件のSOURCE零、独立検算済み554番、今回の20件を合わせて、605件のupper零が確認された。残るopaque12件は579/580/593/594/597/598/605/606/611/613/615/616で、24個のSECOND依存は未決のままである。HF二次の等号、T2全長、H0、nilpotenceについて新たな結論を出さない。

主結果は `reports/AMBIENT_B18_UPPER20_INDEPENDENT_POSTCHECK_V2_20261105_v63/report.json`、実行コードは `code/ambient_B18_upper20_independent_postcheck_v2_20261105_v63.py`。結果reportを変更せず、別のpublication-custody.jsonに全phase・ROOT事前監査・584/554の依存hashを固定した。
