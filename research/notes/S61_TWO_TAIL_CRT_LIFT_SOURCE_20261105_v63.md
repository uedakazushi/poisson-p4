# 固定二tail：40素数 CRT候補と全係数検証の SOURCE

2026-11-05。新しいprime gcd、CRT、候補値、整数因子除算は未実行。
既存32003のgcd3899/94項は別方式peer9041…で受理された有限体診断であり、標数零因子ではない。
元の4500/2317整数とD197、seed、全依存をhashで固定する。

## 既存証明だけで得る上界

任意のprimitive標数零共通因子Hはe次数≤25,d次数≤14で、a646/e19/d1 Kronecker像は単射かつ非定数。
Gaussにより整数商があり、元二式の像のLCが32003で生きるのでH像のLCも生きる。
H像は検証済みgcd3899を割るためKrondegree≤3899、従ってa次数≤6。
a次数6なら19e+d≤23、特にe≤1、e=1ならd≤4。この上界だけで存在は主張しない。

94支持のdecoded boxは(6,9,4)、全次数11、W=2a+eの一意最大13は(6,1,4)。
元二式の最大Wは42/35。これは後述の条件付き整数検証を小さくする支持情報である。

## 固定40素数と内容除去

32003から32363までのSOURCE記載40素数を固定する。全てLC/gamma/Dがunitであることをscalar metadataだけで確認する。
既存32003 seedは再実行せず、残り39だけを将来のworkerに予定する。二整数LCの正のgcd gammaは226 bit。
全積Mは599 bit。これが真のgcdの係数上界を保証するとは言わない。

gcd次数3899ならmonic化しgammaを掛けて3900係数のCRTに加える。
より大きければその固定primeだけskip、小さければSTOPし、次数resetや追加primeを自動で行わない。
固定リスト終了後の候補は一回だけ。M>2gammaを要求し、balance後のleading係数がgammaと一致することを確認する。
全整数係数の正contentを必ず取り、exact除去・逆乗算してprimitive候補へ直す。
**gamma-scaled gcdをそのまま整数除算してはならない。**

候補はseed94項だけに強制せず、32003で消えた係数の復活も許す。ただし実候補の
a≤6/e≤9/d≤4/全次数≤11、一意W13=(6,1,4)という固定形状gateを全係数で検証する。
上限207項。gateに合わなければ候補frontierを保存して停止し、別familyを開始しない。
univariateでの正gcdにはspurious factorがあり得るので、decoded候補には次の三変数検証が必須。

## 偽候補でも有効な整数検証境界

整数除算はW優先・同率lex a/e/d。Hのtailは全てWを落とし、a≤6,d≤4なので、
各更新はa,dを増やさずWを落とす。成功した真の商のe次数や全次数を途中に仮定しない。

- N0：途中Rはa≤16,d≤18,W≤42（8721セル）、Qはa≤10,d≤14,W≤29（3300）。
- N2：途中Rはa≤12,d≤14,W≤35（4680）、Qはa≤6,d≤10,W≤22（1309）。

途中eは42/35に達してよく、整数部ではKroneckerの小radixを使わない。
W*2²⁶＋legacy(a*2²⁰+e*2¹⁰+d)という別の単射keyを使える。

実候補からc=|LC|、T=L1(tail)、BH=bitlenL1(H)、r=min{r≥0:T≤c2ʳ}を求める。
入力L1 bit Bf、m=29/22に対し
\[
B_q=\max(0,B_f-\lfloor\log_2c\rfloor+\lceil\log_2(m+1)\rceil+mr),\quad
B_R=\max(B_f,B_q+BH)+1.
\]
両BR≤8192の場合だけ整数仕事へ進む。これは全途中式とfresh掛け戻しの完了境界であり、
候補が間違っている場合も最初の不正LC/単項式除算で止まるまで有効。
600 bit候補・207項ならBH≤608で、r≤250はBR≤8162という十分条件になる。
これは新しい因子algorithm用の8192条件で、旧518境界を引き上げるものではない。

候補はprimitiveなので、QQで割るなら整数商になる。二本の元三変数式に対する
全係数 H*Q0=N0、H*Q2=N2を必ず確認する。片方だけでは共通因子を主張しない。

## 全所有量と状態保存

元6817 sparse整数recordは一回だけ保持する。field phaseは二uint vector18521＋CRT3900＋input13634＋metadata16000で52055。
field vectorsとdense CRTを解放してから整数R/Qへ進めば最大54090、提案60000、上限150k以内。
新modular product cap60億、bigint product cap210万。45秒soft/60秒hard・2 GBで、完了時間は未保証。

8192 bit coefficientのbinary record≤1040 bytes、候補≤96、入力≤48。
整数prefix上限12914464 bytesを16 MiBで予約し、既存cofactor出力≤5 MiB、CRT/candidate≤1 MiB、report≤1 MiB、
全準備history/入力/code≤3 MiBでunion26 MiB<32 MiB。
Euclid内側j、CRT途中係数i/oldM/newM、content cursor、整数active quotientとH項cursorを全保存し、
完成部分を繰返さない一つのbinary prefixとする。自動retry・prime追加・atomic二重全stateを行わない。

SOURCEは実行可能な条件付き計画であって、40primeでの復元や因子存在の保証ではない。
新worker/custodyの独立ROOT監査とONE GOが必要で、現時点に新値はない。
仮に恒等式が通ってもselected二tailの共通因子のみで、実bad Fittingや全FULLU57/239・補集合・既約性は結論しない。
