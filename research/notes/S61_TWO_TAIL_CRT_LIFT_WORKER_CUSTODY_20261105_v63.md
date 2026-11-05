# 固定40素数 CRT候補／全係数因子検証：未実行 worker custody

2026-11-05。SOURCE6fbe…を実装する準備のみ。
compile、help、実データを開かないsynthetic binary schema fixtureだけを実行する。
新しいprime gcd・CRT候補・整数因子値は未実行。

## 固定した一回のlane

既存32003 seedを再利用し、残り39primeは固定リストだけ。
inputdegree保持を確認し、gcd次数3899より大きければそのprimeだけskip、
小さければSTOPしてresetや別候補を開始しない。
gamma-scaled3900係数をCRTし、全整数contentをexact除去・逆乗算してprimitiveにする。
M>2gammaと先頭gamma一致も要求する。
固定40primeで復元できるという保証ではなく、最終の元三変数式への全恒等式だけが因子証明になる。

実候補の形状207上限と一意W13=(6,1,4)を検査する。
実LC・tailL1から条件付き高さを計算し、両入力の全途中R境界≤8192の場合だけ整数除算へ進む。
途中eを元の32/25に制限せず、a,d,Wの偽候補にも有効な境界を用いる。
整数部はKronecker小radixを使わず、Wとlegacy座標を含む別keyを使う。

## 全状態と所有量

元6817 sparse recordは一回だけ保持。field二vectorは各prime前に古いstorageを解放する。
fieldcapacity18521＋CRT3900＋input13634＋metadata16000で52055。
denseCRTとfieldを解放後の整数部はR/Q/H/inputをすべてState内に置き、最大54090、hard60000。
最初の検証済みcofactor mapは保存後に解放し、二本目と重複保持しない。

`TTCRT001`は64cursor、14 scalar、全40primeledger、source paths、input/f/g/CRTのcapacityと全値、
H/R/Q mapをtyped little-endianで保存する。CRTのoldM/newMと更新途中i、
Euclidの内側jと係数、内容除去位置、整数active quotient・H項・未反映sumを区別する。
時刻STOPの保存中にclock guardを再発火させない。

stored coefficient8192 bitを守り、未commit加算scratchだけは8193 bitまであり得る。
これはState PENDINGに置き、prefix scalarに保存する。record padding1040 bytes以内で、旧518jobsは変更しない。
binaryのlen32を使うが予約済みrecord上限1040/48/96/80を超えない。
prefix読戻しはstreamでcanonical bytesのSHAを再計算し、二つ目の全prefixファイルを作らない。

synthetic fixtureは−3、600bit CRT値、8192bit map値と8193bit pendingを保存するIO検査だけ。
元のN0/N2/seedは開かず、gcd・CRT・因子算術を行わない。
全準備historyとfixtureを含むbase≤3 MiB、prefix16・cofactor5・CRT候補1・report1でSOURCE26 MiB、hard32 MiB。
60億modular／210万bigint product、45秒soft・58秒driver・60秒CPU・2 GB。

一回のvalue jobには新しいROOT source/code/whole auditとONE GOが必要。
仮に二つの全乗算戻しとseed reductionが通れば、既存のmodular quotient-coprime peerから
gcd maximalityを別lemmaで導けるが、実値・独立value peerまでは条件付きのままとする。
Hの既約性、実bad Fitting、補集合、全FULLU57/239の主張はない。
