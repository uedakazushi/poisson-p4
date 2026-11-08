# tangent conic・m二 moving-double：全t源を二つのwhole壁へ

2026-11-08、`/root/quartic_review`。PRIMARY E。
旧QQ(t)全源・全七壁の完全kernelを使い、実shearでhigherを正常化した。
generic全源と五つの壁のactual Dは可約又は零。
**t=−3及び−9/2の二whole族は非因子源として保存し、E帰属／成分性は未解決**。
新しいM1、点、profile又はnormal jetを使わない。

## 1. 既済固定源を再開しない

`PLANE_M2_CONCURRENT_DOUBLE_SEPARATE_REVIEW.md` §2は
`PLANE_M2_DOUBLE_EQUAL_DIAGONAL_NIL_OFF_EXCLUSION.md`を独立再導出し、
tangent h=fの**Z0を含む全source空性**を明記している。
normal k!=0はactual G2=3fku²とnormal位数二から従う。
全raw36/9を保ったratio-Hamiltonianとnilpotent inverseにより、
fraction field C(u,v)で実Fがu²f(base+rational translation)に分解するためである。
したがって旧tangent completion／frontierに残るZ0のB表示は、この既受理適用より古い。
新計算は行わず、この既済の全空性を再利用する。
ZNの全u因子空性、h³のV2 whole essential-three被覆も不変に保つ。

## 2. moving全P0と全source逆

座標(a,b,c,u,v)（旧noteのx,y,h,u,v）で

    Pab=3(ac−b²), Pau=(2t+3)bu−3cv,
    Pav=−3au+(2t+3)bv, Pbu=(t+3)cu,
    Pbv=−3bu+tcv, Puv=(t+3)u²,
    Pc*=0.

actual D0=3t c(bu+cv)²、domainはt!=0。
旧全30/9基底・元全350 J1/J2・全complex例外を不変入力とした。
実shearはa,b,c各々にu,vの線形結合を足す全六つ。
full150 push後のweight一Euler6をraw36上の定数逆で除く。
旧完全K1に対するshear＋残差のmatrix C及びC逆を保存し、
任意の全kernel係数から独立shearと残差へ戻せることを厳密に確認した。
全六shearを常にfreeと仮定せず、t=−6ではrank五として残差二つを保った。

QQ(t)ではrank六、残差零で、det C=−(t+3)²(t+6)⁴。
C/C逆・旧K1・全9P2 forcing逆に現れる全因子は
t,t+3/2,t+3,t+9/2,t+6だけで、全て元七壁＋t0に含まれる。
元壁6,−3/2,−5/2,−3,−9/2,−6,−9を全て元完全kernelで復活し、
新しい複素根補集合やgeneric extrapolationは残していない。

逆shear後も全9P2を未知として元350式を展開し、旧J2最大minorでpivotを解いた。
inhomogeneous particular、homogeneous K2、全残Jacobiを同時に保持する。
特にt=−6ではP2_abのv²係数=−r0²/3が必要で、これを落としていない。

| source | shear rank | residual K1 | free P2 | 元350残条件 | actual Dの結論 |
|---|---:|---:|---:|---|---|
| QQ(t) | 6 | 0 | 0 | なし | `3tc(bu+cv)²` |
| t=6 | 6 | 0 | 0 | なし | `18c(bu+cv)²` |
| t=−3/2 | 6 | 1 | 0 | なし | `−c(−3bu−3cv+r0v²)²/2` |
| t=−5/2 | 6 | 0 | 0 | なし | `−15c(bu+cv)²/2` |
| t=−3 | 6 | 2 | 0 | なし | 下記A、未帰属 |
| t=−9/2 | 6 | 0 | 2 | なし | 下記B、未帰属 |
| t=−6 | 5 | 2 | 1 | `r0 r1=0` | 全c因子 |
| t=−9 | 6 | 1 | 0 | なし | 全c因子 |

proper因子の列は全residual値を含むactual Fの恒等式である。
表示した元350残条件は必要十分な式で、全linear supersourceを無条件Poissonと呼ばない。

## 3. 残A：t=−3の全二parameter Poisson族

r,sを旧residualのr0,r1として、全normalized源は次の十括弧。

    Pab=3(ac−b²), Pac=a(ru+sv),
    Pau=−3(bu+cv), Pav=−3(au+bv),
    Pbc=asu+crv, Pbu=v(ru+sv),
    Pbv=−3(bu+cv), Pcu=Pcv=Puv=0.

全(r,s)で元350式は恒等零、全9P2は零。
actual Dはreportに全展開して保存し、generic(r,s)では単純なproper factorを得ていない。
s=0での因子式は

    D=−[r a u²+3c(bu+cv)] [3(bu+cv)+r v²].

したがって既約quinticを調べるならs!=0だけが必要だが、全s0源を保存する。
この二parameter族へのfull-source逆は上記C逆であり、chosen点の局所モデルではない。

## 4. 残B：t=−9/2の全二parameter Poisson族

p=q3、q=q6とし、全normalized源は

    Pab=3ac−3b²−p uv−9qv²,
    Pac=u(pu+8qv), Pau=−3(2bu+cv), Pav=−3(au+2bv),
    Pbc=q u², Pbu=−3cu/2, Pbv=−3(2bu+3cv)/2,
    Pcu=Pcv=0, Puv=−3u²/2.

全(p,q)で元350式は恒等零。
これは全9P2を先に戻してから、q0=0,q1=−p,q2=−9q,q4=8q,q5=q7=q8=0とした結果である。
actual Dは

    −(3/2)[−aq u⁴+9c(bu+cv)²+bp u⁴
             +8bq u³v+cp u³v+9cq u²v²].

q=0では

    D=−(3/2)(bu+cv)[9c(bu+cv)+p u³].

よって既約quinticにはq!=0が必要。全q0補集合を保持する。
この族のactual normalizerとknown/component認識は次の純構造義務であり、
新しいambient／profile又はjetを本jobでは追加しない。

## 5. 実行と欠陥の正直な保存

V1はsaved RF辞書をSymPy Matrix作成前にdecodeせず、source算術前に停止した。
V2は同じ八sourceの全算術を完了したが、factor exponentのSymPy IntegerをJSONへ書けず停止した。
V1/V2 input/job全文は不変に保持し、各failure receiptを別identityで保存した。
V2の不完全な出力はJSON reportと呼ばず、exact bytesを`.txt`のpartialへ保存した。

**NEW V3**はRF entryの事前decodeとfactor exponentのint化だけを修正し、
同じQQ(t)＋元七壁を再実行した。点、M1、方程式又はprofileは追加していない。
V3 gate PASS、4.152128秒、RSS66,988KiB。
report `CLASSIFICATION_PLANE_M2_TANGENT_MOVING_WHOLE_SOURCE_20261108.json`は
exact実行program、全C/C逆、全残parameter、全9P2解、元350全残差、actual Dと全input hashを保存する。
新算術・全源逆の別作者reviewは本稿作成時には未了。

全接線m二E完了や新成分をまだ宣言しない。
既済hf／h³／bad-columnに、今回のgeneric＋五壁空性を足し、残A/Bを全源で認識する必要がある。
