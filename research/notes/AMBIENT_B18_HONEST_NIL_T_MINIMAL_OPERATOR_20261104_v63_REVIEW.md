# honest T_nil の全次数5×849とd3既知一次部分：独立監査

2026-11-04。既受理128×972の全actual sourceとrank24/18/81の固定証明を使い、
全次数の三段単元Schur縮約を監査する。新source値・jet・X2・m1は計算しない。

## 完全性とminimal target

Eをm0の選択24列とするとmEはoriginで単元であり、n=E(mE)^-1はmn=Iを満たす。
このsectionはE座標以外で文字通り零なので、最初の24section列のSchur消去のlower blockは
全次数で零である。残る104×948の上20行は元Xの補集合座標、下84行は元の形式函数h。
この新sectionも前の完全比較証明に適用でき、単なる定数行列への置換ではない。

次の18×18 origin minorは、この20座標に適合する保存固定minorへ全780 factorizationを
照合した。最後の81×81 minorは定数相殺後の保存W84から来る。各blockは非零QQ origin
determinantを持つのでSの全次数で単元となる。全60 typed matrix nodeと三段の
D−C U^-1 B、target projection PR−C U^-1 PIを確認し、自由非可換block恒等式で
単元消去がcokernelを保存することを確認した。

全972元関係の経路は重複なく、24+18+81列が単元対として消え、681 actual X列と全168 d3列が
849列として残る。全原点係数は既受理rank証明により零となり、残る五targetはminimal。
targetは元free44の42,43とsector0 wedge(0,7),(4,7),(5,7)。従って保存作用素は
同cokernelの**完全5×849形式函数表示**である。全849scalar級数を展開したという意味ではない。
T_nilの生成元を元nilidealの全生成元と同一視しない。

## d3だけの一次関係

d3は原点零で、最初の24sectionと次の18消去におけるupper targetは全次数で零。
従って最後の81消去後のd3一次項は、上二targetで零、下三targetでは
保存constant quotient projection3×84とliteral d3一次項の積だけで定まる。
全168列を独立に再構成し、固定15minorの両側逆・行列式と全列factorizationから
この既知部分のrankが正確に15であることを確認した。全一次関係空間とは主張しない。
単独の帰結はgr1次元≤40−15=25である。

## 未知項を消さない

X0のconstant residualではm2 X0は消えるが、二次入力m0 X2+m1 X1のm1項は消えない。
18pivotの一次消去にもh_pivot,0が現れるので、定数生成数の計算をそのまま一次層へ延長しない。
実old Bold0の残差を保存cacheから再構成すると九つが非零で、supportは旧48座標
[32,36,40,44]、span rankは固定4minorで正確に四。従ってadapter二次項を全て捨てる仮定は
認めない。既存B11が保存しない二次係数と、第三sectorのpure r³も未計算のまま保持する。

総長、全一次Loewy層、元nil生成元、nil三乗をこのpacketから結論しない。

再現: `python code/ambient_B18_honest_nil_T_minimal_operator_review_20261104_v63.py`。
