# constant-root m三の六残源：全Ore支持の幾何還元

2026-11-07、主義務G。producer手証明。全六のoperator/Ore source再構成は
`CLASSIFICATION_PLANE_NONREDUCED_M3_CONSTANT_ORE_20261107.json` に固定する。
単純な二次Jacobi coefficient rowspanによるlinear-factor判定は全源でnegativeだった。
これはfactor不存在の証明ではなく、以下では全reduced supportの幾何を使う。
独立review・残moduleは未了。全constant-root又は全m三完了はまだ主張しない。

## 1. 全sourceのOre橋

前十四sourceのP1全必要kernelと**FULL P2九**を戻す。
二段のEuler gaugeにより全tensorはc非依存の

`Π=Q(y)+∂c∧W(y), y=(a,b,z,w)`

へ入る。第一gaugeはL=(-∂cΠac,−∂cΠbc,0,0,0)。
残るL_c=ℓ(z,w)でW_zのz²及びW_wのw²を零にする。
全firstweight sourceがこのgaugeへ入ることを係数恒等式で検査し、
Qの全六bracket×十monomial coefficient matrixのcomplete rankとbasisを保存した。
Wは固定W0に、Wa/Wb各七、Wzのzw/w²、Wwのz²/zwの全十八係数を加えたものに含まれる。
QとWを独立に取ることはcovering supersetである。
元sourceの相関条件は落としても必要条件となるが、特にG3 lineを使う時には必ず元条件を戻す。

c次数を元projective Jacobiで比較して、全actual solutionには

`[Q,Q]=0, E_y∧[Q,W]=0`

が必要で、c非依存表示では逆も成り立つ。
これはglobal divergenceでstabilizerを選ぶshortcutではない。

## 2. 二EXACT_B whole lineはessential三

DOUBLE_LINE_E0_EXACT_Bの全Qは、全αで

`Q=∂a∧(F(z,w)∂b+X(z,w))`

となり、a非依存である。
DOUBLE_LINE_E2_EXACT_Bでは同様に

`Q=∂b∧(−F(z,w)∂a+u zw∂w)`

でb非依存である。generic fieldのbasisに分母が出るα値も元complete kernelを戻して確認し、
このconstant-factor性は各壁でも保たれる。

一般にQ=v∧Y、v定数非零、L_vQ=0ならPfQ零でactual Dはcに依存しない。
Dは符号を除きE_y∧W∧Qの係数。v微分し、module式E_y∧[Q,W]=0にvをwedgeするとv(D)=0。
従って非零既約Dは二独立方向c,vに依存せず、本質三変数の既受理八成分へ入る。
全ambient像上界二十三をglobalに使い、plane/v選択次元を再加算しない。
Q零ならD零で既約quinticの対象外である。
これによりE0/E2 EXACT_Bの両**whole α line**はG非優越まで閉じる。
constant-factor補題は `PLANE_DEPENDENT_ORE_SEPARATE_REVIEW.md` §3に独立証明済み。

## 3. TRIPLE_LINE_SQUARE_PLUS_LINEARは全枝を閉じる

receiptのQ-coordinatesをq0,…,q7、central F=q5z²+q6zw+q7w²とする。
全bracketsは

```
Qab=2q0az+(2q1−3q2)aw+F,
Qaz=0, Qaw=q3zw,
Qbz=−q0z²+(q2−q1)zw,
Qbw=−q0zw−q1w²+q4zw+3aw, Qzw=0.
```

strict Jacobiの全idealは `q3(q1−q2), q0q3` である。
complex point/reduced supportについて二分する。
q3零ではQはb非依存の∂b定数因子を持ち、上記global essential三へ入る。
q3非零ではq0零・q1=q2で、Qのz-column全零となる。

元G3 line `C*a²z²w` の条件をここで戻す。
全W十八sourceのうちWz=ρzw+σw²とすると、actual D3のa²w³係数は正確に3σ。
他のQ/W係数からこのmonomialへの寄与は零なのでσ零が必要である。
従ってΠのz-columnは全z倍、omitted-zのPfaffianもEulerのz倍を持つ。
actual Dはzで割れ、既約quinticには矛盾する。
これでq3零・非零の全補集合が閉じる。
この枝でW十八supersetのline制約を省略して結論を出していない。

## 4. TRIPLE_LINE_LINEARはconstant factor又は四fixed Q

全Qはreceiptの座標で

```
Qab=(−2q1+q2)aw−q1bw+2q3az+F,
Qaz=−q1zw, Qaw=−2q1w²+q3zw,
Qbz=−q0z²+(q1−q2)zw,
Qbw=−q0zw+(q1−2q2)w²+q4zw, Qzw=0.
```

strict Jacobiのw²z係数に−2q1²がある。reduced supportではq1零。
残る全式は `q2q3=0, q3(q0−2q3)=0`。
q3零では再びb非依存の∂b定数因子である。
q3非零ならq2零、q0=2q3。
actual b→b+h aは元P0を保ち、q4→q4+h q3なのでq4を零へ直せる。
actual diagonal/tensor scaleでq3を一へ戻せる。
actual b→b+D wはcentral Fのzwを任意に変え、z²/w²は変えない。
残二central coefficientsは独立nonzero diagonal scalesにより各零又は一へ直る。
全inverseはnonzero diagonal inverses、minus shearであり、全零branchも別に保持する。

従ってnonconstant-factor全源はε,η∈{0,1}の四fixed Q

```
Qab=2az+εz²+ηw², Qaw=zw,
Qbz=−2z², Qbw=−2zw, Qaz=Qzw=0
```

と、全W0+W十八のcoverへ入る。W0=(0,−3a²,0,−3aw)。
全固定moduleの原Jacobi/actualDによる結論は次の有限義務であり、ここでは仮定しない。

## 5. E1 EXACT_Bのwhole αは四whole Q line

同じ全sourceのnoncommon Q coefficientをuに取り直すと、全αで分母なしに

```
Qab=−u az+F, Qaz=0, Qaw=−u zw,
Qbz=u z², Qbw=u(α+2/3)zw, Qzw=0
```

となる。旧basisがpoleを持ったα=1/3も、元kernelはこの四次元spaceに一致する。
actual D3=3u ab z²wだから対象m三ではu非零。
actual diagonal a'=sa a,b'=sa b,c'=sc c,z'=sz z,w'=sw wとtensor scale sa/scはP0を保ち、
uをsa²/(sc sz)倍する。従ってu一へ直せる。
b→b+D wでFのzwを零へ戻す。
残f z²+g w²は独立diagonalsで各零/一へ直り、全補集合は四whole α lineで尽くされる。
全αは保持し、特殊値にgeneric rankを外挿しない。

正規化後は上記u一、F=εz²+ηw²、ε,η∈{0,1}、
W0=(a²,−2ab,az,3αaw−3bw)と全W十八である。
全moduleのfactor/上界・全α分母補集合はまだ未証明。

## 6. 残義務

本noteで六残源のうちE0/E2 EXACT_B、TRIPLE_SQUARE_PLUS_LINEARは全G非優越。
TRIPLE_LINEARはconstant-factor枝と四fixed Qへ、E1 EXACT_Bは四whole α Q lineへ還元した。
E2_LOGの全Q strict-support分岐、特殊αのextra basis、全module結論は別に残す。
旧六row-factor criterionのfalse・十四linear粗上界のfalseは不変。
他primitive row degrees、generic triple root、m二、reduced Rは本source外。
