# 有限Tor63像：全38 first relationの一次倍による追加cheap stop gate

2026-11-04。前packetの十五exact targetは一つも再計算せず、その保存405係数をそのまま使う。
ここでも全640 second carriesやregular quotientの未知二次微分は作らない。

## 1. Actual first関係を全て使う小さいsubset

既独立受理されたactual D0の全640 first columnsは81係数で保存されており、QQ rank38である。
新checkerはこのrankと全640列のbasis係数を再検証する。
独立38列に九つのF9変数を掛けた342列を、free C0九generatorの二次405係数空間へ入れる。
これはactual relationの実際の一次倍だけなので、full quadratic relation initialの正当なsubsetである。

十五saved targetをこのspanで減らし、original342列へのmembership係数とexact residualを全て保存する。
intrinsic28だけの前検査で残った七本も含め、同じsource-bound targetを用いる。

## 2. 結論の範囲

計算結果は342列のspan rank **261**で、十五targetは**全てこのspanに入った**。
全十五本のexact membership係数を保存した。従ってこの二次経路は以下のstop条件で停止する。

全targetがこのspanへ入れば、選んだactual一次relationの一次倍を引いて、十五q-wedge-q像はm³D0に入る。
十二eta-wedge-q像も既にorder3以上なので、有限Tor63像全体はm³D0に含まれる。
その場合quadratic nonmembership経路はinconclusiveとして停止する。
有限像が零であること、socleが無いこと、元Dのdepth20は従わない。

残余があっても、まだfull second cancellation carriesを含めていないため、実際の非零像やdepth19は従わない。
必要な次の論理gateは、actual全second initialへの非所属を保証する上界である。
このpacketからsecond-carry worker、高次jet、source5やliteral transformを自動的に起動しない。

Code `code/ambient_B18_Tor63_full_first_multiple_stop_20261104_v63.py`。
Report/folder `AMBIENT_B18_TOR63_FULL_FIRST_MULTIPLE_STOP_20261104_v63`。
