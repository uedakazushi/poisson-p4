# Actual B18の最上w-torsion T3は零

2026-11-04。これはactual full relative入力の**T3=0**を証明する。
低次jetのmembership試験が陰性だったことだけから零性を推論しない。
literal tailと同R0 N moduleが与える全kernelの44生成上限を使い、43個の独立な一次項と
残る一個の二次項によって、有限jetから全次数のNakayamaへ移す。
T0,T1,T2と元B18全Dのdepth19/20は未判定である。

## 1. 正確なmoduleと実際の二つの全次数入力

R=R0は実際のregular完備局所九変数環、mはその最大idealとする。
既独立受理のhonest9表示を

```
M=D3=R^9/Q3,
N=D4=R^9/J,       J=Q3+(top12),
T3=J/Q3
```

と書く。全36 Koszulを保持した表示から単位10/3を消したものであり、
二行trace ambientの像に置き換えていない。mu(M)=mu(N)=9。

literal tail bindingと独立reviewは、同じR上で

```
0 -> L -> N -> A_N -> 0,
L=(R/P)*F,      R/P=QQ[[Z]],
```

を与える。Pはregular8生成。Nのnilmodule全体がこのcyclic moduleであり、
whole Nを一変数二重環と置き換えたものではない。

同R0 N bindingとliteral coneの独立reviewは、actual formal coordinatesを

```
R=R8[[v]], R8=QQ[[z1,...,z8]],
A_N=R8+m8*tau, tau^2=nu in m8,
v=z1*tau in A_N
```

として指定する。z1=b61のnilcurve方向微分は非零なので、Z=z1をLのparameterに取れる。
従って、R内で

```
P=(z2-psi2(Z),...,z8-psi8(Z),v-psiv(Z)).
```

この全次数結合はAc型の一致やsource rankから仮定したものではない。
同じphysical R0とunit56/8入力が独立検証されている。

## 2. Jの44生成と、43個の一次項

A_Nの八つのminimal module generatorsをn0=1,n1,...,n7とする。
literal coneの全36関係は、次の形である。

* 七つの `z1*n_j-v*z_(j+1)*n0`。
* 残るz2,...,z8の21 Koszul関係。
* 七つの `v*n_j-nu*z1*z_(j+1)*n0`。
* 最後の一つ `r0_N=(v^2-nu*z1^2)*n0`。

最初の35個の一次vectorは独立である。
七つのz1*n_j、七つのv*n_j、21個のz_i*n_j-z_j*n_iは、異なるmonomial/外積成分で
独立なことが、保存した普遍行列から直接分かる。
これら全てのn0成分はm²内であり、r0_Nのn0二次項はv²である。

この八generatorをactual physical functionsとしてN=D4へ持ち上げ、Fを加える。
これら九個はNを生成する。mu(N)=9なので、この九個はminimalであり、元honest9との
generator変換はR上のunit9行列である。

上の36関係を持ち上げた値は、一意にc_alpha(Z)*Fとなる。
従ってfree module R*F plus sum R*n_jからNへのkernelは、全て

```
P_i*F (i=1,...,8),
r_alpha_N-c_alpha(Z)*F (alpha=1,...,36)
```

という44個で生成される。これはextension完全列から直接従う。
任意のkernel vectorをA_Nへ射影して36関係で引き、残りをann(F)=Pで引けばよい。
36個のc_alphaを数値計算したり、extensionの分裂を仮定したりする必要はない。
mu(N)=9により、c_alpha(0)=0。

nil側の八関係の一次項はPの八独立linear formsにFを掛けたもの。
N側の35関係の一次項はN成分が独立なので、合計43一次項は常に独立である。

独立に再構築した実際のJの一次spanのrankは**43**である（§4）。
残るrelationを

```
r0=(v^2-nu*z1^2)*n0-c0(Z)*F
```

とする。もしc0'(0)が非零なら、その一次項は非零なZ*Fである。
Zはnilcurve parameterなので、Pの一次spanに入らない。
N成分を零に保って35個の独立N一次項を使うこともできない。
従ってrankは44になり矛盾する。よってc0 in (Z²)、r0の一次項は零。

以上の43一次独立generatorとr0は全44個を尽くす。
r0はm*Jに入らない。なぜなら全generatorのn0成分はm²内なので、m*Jのn0成分はm³内だが、
r0のn0二次項はv²!=0だからである。
特にJ/mJの一次kernelはr0のclassによる一次元だけである。

## 3. top12をmJへ入れる二つの異なるgate

独立な全九変数計算とorder検証により、元honest9で全top vector tは

```
t in m²*R^9,       its old E0 component in m³.
```

を満たす（§4）。
上の新しいgenerator basisで、n0係数rowの原点値は、元honest9 basisに対して
`(1,0,...,0)`である。
元E0はtargetの単位1、他八個はactual CA odd functionで点で零。
新しいn_jとFも点で零なので、ring residueを取るとこの原点rowが決まる。
ここでmodule fiberとring residueを同一視しているわけではない。
R上のunit generator変換のn0係数だけを、ring residueにより測っている。

従ってt in m²なので、変換行列の正order部からn0成分への寄与はm³に入る。
新しいbasisでも、全topの一次項は零、n0二次項も零である。

任意のtを全44 kernel generatorsで書く。
一次項が零で、最初の43一次項が独立だから、それら43のconstant coefficientは全て零。
従って

```
t = lambda*r0 modulo mJ,       lambda in QQ.
```

n0二次項を比較すると左辺は零、右辺はlambda*v²なのでlambda=0。
全topはmJに入る。よって

```
J=Q3+mJ,
J/Q3=m*(J/Q3).
```

J/Q3は有限生成R-moduleなのでNakayamaから**J=Q3、T3=0**。
finite-jet membershipを逆向きに使ったのではなく、全kernelの44生成構造と
異なる次数の独立な識別子を使った全次数の証明である。

## 4. 今回の独立算術検証

`AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63`の保存データを作者実行なしで照合した。
元pure-unit scalar DAGを独立ADし、CA/CBの各57×57 Jacobianと逆行列を再構築。
保存された各57 source seriesを、作者の全monomial convolutionとは異なる
constant/linear/quadratic blockの積公式で評価した。
両方の57式について全55係数が正しいこと、全pure function2係数が一致することを確認した。
新しいNewton/IFT solveは行わず、保存された二次係数の検証だけを行った。

同じphysical x9を通じたCB→CA adapterとdet1/4096を再構築した。
全792新Hadamard係数、全784旧coordinate係数のfirstjetを独立に戻した。
四つのg商の微分は独立な未知正order函数のまま保持し、全2069関係とtop12から
それらの感度が恒等的に消えることを再確認した。未知微分を零へ置いていない。
全2069×90 firstjet係数を再構築して凍結certificateと照合し、独立疎Gaussでrank43を得た。
全topのfull9 firstjetは零である。

さらに、原点値だけを用いるtyped order propagationでも全top entriesがm²に入り、
E0が少なくともm³に入ることを検証した。
旧honest10の全12 topでunit消去対象component5がliteral node0なので、Schur後のE0は
raw E0と同じpointerである。これは未知Hadamard微分やnil抽出器に依存しない。
高いjetを計算していない。

## 5. 範囲

**actual T3=0、Q3=Jinf**を証明した。
actual w-torsionの他の三つのgraded部分T0,T1,T2の零性はまだ判定していない。
従って元B18全Dのdepth19/20、全nilradical、primary構造の解決は主張しない。
個別成分や他のB18点への無断transportも無い。
36個のactual extension coefficientの計算は、今後他のgraded部分を調べるためには有用だが、
今回のT3零性の必要条件ではない。

Code `code/ambient_B18_T3_zero_structural_20261104_v63.py`。
Report/folder `AMBIENT_B18_T3_ZERO_STRUCTURAL_20261104_v63`。
独立再構築したfirstjet matrix、43 pivot、全top order gate、構造的適用を保存する。
