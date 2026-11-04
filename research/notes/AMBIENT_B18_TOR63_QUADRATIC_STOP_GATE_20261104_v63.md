# 有限Tor63像のcheap quadratic stop gate：intrinsic C0 Koszulだけとの比較

2026-11-04。前preflightの有限moduleの像が非零かどうかを、まず最も安い二次検査で調べる。
元D0の全640 second carriesを作らず、完全なintrinsic28 Koszulの一次式のmaxideal倍だけを使う。

## 1. Exactな十五quadratic target

六actual cycleの原点q_alpha(0)と、source-bound同R0のeven差a_A-a_Bを使う。
cached E4 source2をCAのfull F9座標へ正確に戻すので、二族の異なるF9を直接引かない。
even差のconstant/linearは全て零である。従って

```
target_alpha,beta,2
 =sum_i,j (q_alpha,Ai*q_beta,Bj-q_beta,Ai*q_alpha,Bj)(0)
       *(a_Aj-a_Bj)_2 * origin(image(f_i*tauA))
```

は、そのactual Tor2像のC0 free9における二次成分全体である。
even成分はnu_A*f_i*f_jからorder3以上なので二次へ寄与しない。
残り十二eta-wedge-q像もorder3以上なので、今回は十五targetだけで十分。

## 2. 完全intrinsic relationの正当な小span

C0 basisは1,z7*tauA,gA*tauAで、gAはorder3。
そのmodule kernelは八regular generatorsの全28 Koszulである。
一次式を直接作り、凍結済みactual D0 first matrixの最初28列と全係数照合する。
九変数をそれぞれ掛けた252列は、9*45=405行のQQ matrixになる。
十五targetを疎Gaussで減らし、全original252列に対するmembership係数または残余を保存する。
大きいfull640 second-carry matrixやopaque未知二次微分のmatrixは作らない。

## 3. Stop条件とscope

計算ではintrinsic spanのrankは**196**。十五targetのうち八本はそのspanに入ったが、
cycle pairs(1,2),(1,3),(1,4),(1,5),(2,5),(3,5),(4,5)の七本は入らなかった。
残余はzero7成分・zero7 monomialだけを使い、gA成分やr,sは含まない。
特にq_f63を含む五本は全てintrinsic spanで消えた。
全membership係数と七つのexact residualを保存した。

十五targetが全てこのspanに入れば、intrinsic actual relationの一次倍を引くことで全27 minimal代表はm³D0に入る。
従って有限Tor63像全体もm³D0に含まれる。
しかし有限socleが高いm-orderに隠れることは可能であり、**これは像が零という証明ではない**。
この場合quadratic nonmembership経路はinconclusiveとして停止し、自動的な高次jetや全second-carry jobへ進まない。

逆にintrinsic spanの外に出ても、full D0の追加second relationsがあり得るため、それだけで非零像やdepth19を主張しない。
必要なら別途、actual全second relationの上界を証明しなければならない。
新source jet、source5、literal九105-clone transform、GB、point、primeは使わない。

Code `code/ambient_B18_Tor63_quadratic_stop_gate_20261104_v63.py`。
Report/folder `AMBIENT_B18_TOR63_QUADRATIC_STOP_GATE_20261104_v63`。
