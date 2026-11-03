# NEW25 の8次特性基底への優越性：独立監査

2026-11-03。対象 `Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63`。
同じ既存 source 点[0]だけで、全source微分と普遍係数変換を独立に検査した。
新点・degree24以上の評価・author codeの再実行は行っていない。

## 1. 全61×18微分を逆方向ADで再構成した

6079-node dense mapの全値をFractionで作り、全逆数nodeが定義されることを確認した。
authorの18-source forward ADとは異なり、各61出力から一つずつreverse sweepを行った。
積の局所微分はprefix/suffixで作り、零因子がある場合にも割算でlog微分する誤りを避けた。
全61値と全1098微分係数が保存証明書に一致する。
value pivot=2051/12、Pfaffian=117727/2592も同じ既存点の値である。

primitive変換は、旧S61の厳密なprimitive-beta係数式から24×16線形写像を別に作り、
元の24×16基底の独立な左逆を用いて16×16変換を再構成した。
全24係数の写像等号、beta16の全値、全16×18 source微分を照合した。
従って、この微分は抽象的なW上の微分ではなく、受理済み actual Poisson sourceへ結び付いている。

## 2. 五つの8次traceの微分とQQ minor

使用するIは4×4のC,Dからの
tr C4,tr C3D,tr C2D2,tr CDCD,tr CD3である。
先に独立受理した全二次共変量係数を直接微分し、各beta方向のdC,dDを得た。
4×4のforward dual matrix productから5×16 beta gradientを作り、
authorのpair-gradient×pair-Jacobianとは別の入口で一致を確認した。
さらに全5×20 pair gradient、5×18 actual source微分も照合した。

10本のSp4 source列は全係数が厳密に0である。これは実装の整合検査であり、
普遍的なSp4不変性そのものは既存の共変量構成とtraceの共役不変性に依存する。
選択列(n0,n1,n2,scale,base1)の5×5 minorは独立Fraction消去で同一の非零値となった。
分子461bit・分母329bitの完全値を保存した。rank5はQQ上の事実であり、modular rankの外挿ではない。

## 3. 全係数の基底変換

Iからpencilの五係数pへの保存5×5定数行列を独立に読み、det=−1/128を再計算した。
既存に全係数で独立受理された二つの表を使い、5式すべてについて

`change * (256 I coefficient dictionary) =256*(p coefficient dictionary)`

を全単項式で確認した。256は共通の分母払いであり、不変式の尺度を混同していない。
従ってIとpは同じ5次元affine座標空間を可逆線形に記述する。

## 4. 優越性とfull characteristic base

actual sourceの非零5次minorから、char0の微分判定により五つのIはsource上で代数的独立である。
sourceの像は受理済みNEW25 first-jet像B_NEWに入るため、B_NEWもA5_Iを優越する。
可逆な係数変換によりA5_pも優越する。B_NEWの滑らかさをこの一点で仮定する必要はない。

F4も含むfull characteristic baseは

`B={F4²=400p0+1600p4}`。

この普遍係数恒等式には既存Spin5独立監査を使う。
p0=(F4²−1600p4)/400として消去すればBはA5と同型で、既約・次元5である。
p-spaceへの射影は、F4のmonic二次式による有限写像である。

B_NEWのB内の像の閉包は既約で、p-spaceを優越するので次元は少なくとも5。
既約5次元Bの真の閉部分集合は次元4以下だから、この閉包はB全体となる。
従って **NEW25 first-jet像はfull characteristic base Bをも優越する**。

よってB上で恒等的零ではない、base係数だけのpolynomial条件の零集合にNEW25全体を入れることはできない。
同じことはBの任意の真の閉部分集合に対して成り立つ。
特性pencilやそのbranch quarticの係数だけから作る真のdivisorは、NEW25方程式にならない。
既存base関係を使って恒等的零となる式は、そもそもdivisorではないと区別する。
この結論は一つの高次数判別式の標本評価によるものではなく、base全体への優越性による。

## 5. multisectionという範囲

既受理のregular-D/gauge chartはNEW25にも交わり、その有限軌道切片では周囲のrational quotientは
6次元、NEW25のrational quotientは5次元である。
前者のBへのdominant mapの一般fiberは1次元、後者のBへの関数体拡大は同じ超越次数なので有限となる。
したがって、このような局所的なrational quotientモデルでNEW25は残る曲線方向の
generically finite multisectionとなる。有限Weyl coverでその引戻しが既約とは仮定しない。

multisectionの次数、曲線の種数、商の有理性、全域での有限性・properness・方程式は未決である。
F4と全Inv8だけの式がNEW25を定義しないことと、全invariant環がその部分環で生成されることは異なる。
今回24次基底や候補式を一点にも評価していない。

入口 `code/q_new_family_octic_base_dominance_review_20261103_v63.py`。
60秒・2GB cap、凍結primary・原稿・台帳は変更しない。
