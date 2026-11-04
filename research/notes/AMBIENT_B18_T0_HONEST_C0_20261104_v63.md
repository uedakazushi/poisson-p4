# 残るT0のactual入力：C0 free9/28、Q0全612列、D4へのmap

2026-11-04。T1,T2,T3の消滅後に残るT0を調べるため、degree0の実環D0を
正確な有限module presentationへ戻した。T0のkernel、socle、full Dのdepthは未計算。

## 1. Regular sequence八本のC0

共通physical R0は九変数regular局所環。actual CA F9の最初七つをz1,...,z7、
残りをr,sとする。f63をz=0に制限したg_A(r,s)は非零でorder3。
従って(z7,g_A)はR0のregular sequence八本であり、

```
C0=R0+(z1,...,z7,g_A)*tau_A
```

のfree9 coverは、1,z1*tau_A,...,z7*tau_A,g_A*tau_Aをgeneratorとする。
normalizerのbasis(1,tau_A)がfreeなので、その全kernelはregular sequenceの
28本のKoszul first syzygyで生成される。全28列を保存した。

## 2. 新しいsource計算無しで四つの正則商の原点値を決める

degree0のCA zero graph rowsは64,68,73,74の四つ。
各odd coefficient B_alphaは、actual zero-image algebraの定義から
I0_CA=(z7,f63)に入る。z=0に制限すると

```
Bbar_alpha = g_A*H_alpha,
H_alpha in QQ[[r,s]] algebraic and regular.
```

その正則性は全次数のsource-image membershipに基づく。
既保存のcombined105三次解を使って元source DAGをこの五関数だけ再評価した。
g_Aの三次初期形は既存値と一致し、四つのBbarは三次まで全て零。
従って既知の正則除算恒等式の最低次成分から、**四つのH_alpha(0)は全て0**と厳密に決まる。
jetから全membershipを推測していないし、H_alphaが恒等零とも主張しない。

## 3. 七つのpure CA prefix系と35商

r,sを常に保持し、zの最初k個を残すCA57系をk=0,...,6について構成した。
z七つ全部を残すk=7は元のsourceとpointed unit一意性で一致する。
各系のunit determinantは既受理のCA57 determinantである。
新しいunit Hensel入力を保存しただけで、source IFTのNewton solveや新しいjetは実行しない。

f63と四つのB_alphaに対して各七coordinate quotientを作り、全35本を保存した。
prefixがz_i=0で一致するので、その商は全次数でregularである。
原点値は既存FULL9 source二次証明書の一次係数を読み戻した。

```
f63 = g_A + sum zi*e_i,
B_alpha = g_A*H_alpha + sum zi*h_alpha_i.
```

この式はtelescopingにより全次数で成立する。
四つのH_alphaは上記の正則g商nodeとして保存し、今回決めた原点値0を付した。
旧CB側のr/g関数は、元のactual R0 scalar cofactorのまま保持しており、
それらをcompose・再除算・未知微分0へ置き換える操作は無い。

## 4. 全612列のliftとD0の正確な表示

Q0はdegree0のCB kernel64本とzero graph差4本をC0の九generatorで掛けた全612列。
odd係数はCA zero八関数と上記四graph関数の和・積としてfactorized DAGに現れる。
前のhonest10と同じ、和・負号・既知ideal factorへのscalar乗算の厳密な規則で
全612係数を(z7,g_A)に持ち上げた。二重負号のDAG簡約も同じ正確な符号規則で扱う。

従って

```
D0 = coker(R0^640 -> R0^9),
640列 = 完全28Koszul + 全612actual Q0 lifts。
```

これは全degree0 idealのpresentationであり、positive weightの関係を勝手に省いた近似ではない。
非負gradingなので、degree0部分に寄与するのは正確にこれら68 relationだけである。
全640列の原点matrixのrankは0。従ってmu_R0(D0)=9。
元の二行matrixをidealで減らすTorを捨てた表示ではない。

## 5. 同じD4への明示map

既存honest10はf63*tau_Aの正確なH-basis liftを持つ。
そこからsum e_i*(zi*tau_A)を引くことで、g_A*tau_AのH free10 liftを得た。
1と七つのzi*tau_Aは直接のbasis vectorである。
これへ既存の単位Schur射影9×10を適用し、C0 free9からD4 honest9への
全9×9 mapを構成した。
各列は同じactual trace pairを表し、H coverの全Koszul kernelもD4 presentationに入るため、
mapはD0へ降りる。

その原点rankは7だった。ただしこれはT0のfiber次元、非零性、supportの判定ではない。
正確に

```
T0 = ker(D0 -> D4)
```

を指定する有限入力が得られた段階であり、このkernel自体は計算していない。
T0をnilcurveだけのmoduleと仮定せず、rank2被約成分の寄与も保持する。

## 保存・実測

`code/ambient_B18_T0_honest_C0_20261104_v63.py`、
`reports/AMBIENT_B18_T0_HONEST_C0_20261104_v63/honest-C0-D0-D4-DAG.json`。
全567193 scalar node、七57系、35coordinate商、四正則g商、全28+612列とD4 mapを保存。
証明書11092247bytes。一回1.705280秒、最大RSS430136KiB。
60秒/2GB/32MiB内。GB、新点、新しいsource IFT solve/jet、old nonunit pullbackは無い。
報告metadataの「新しいIFT無し」はsolveが無い意味と明確化し、
literalな新unit入力七系が存在することは別fieldと本節に明記した。
