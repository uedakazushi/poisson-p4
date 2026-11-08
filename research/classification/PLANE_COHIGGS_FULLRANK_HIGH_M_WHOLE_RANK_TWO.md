# Co-Higgs rank M三・高multiplicityの全二源はrank二以下

2026-11-08、PRIMARY E、log_manuscript。producer全source純証明とONE有限readback。
独立受理・global採録はrootの別義務。新成分、全co-Higgs又は全Eの完成は未主張。
1054主稿／frontier、旧runner／旧契約／旧証明書／Gitは不変。

## 1. 全源定理

rank零planeを持つPoisson tensorで、共通co-Higgs imageのsaturationがO(1)、
M:C3→sl2のrankが三、実五次のplaneに沿う位数が三以上なら、
全actual projective tensorはrank二以下である。
従って指定するrank四・被約既約quintic sourceは空である。

全multiplicity local modular入口と全actual Lax normalizerは
[先行全m独立review](PLANE_HIGH_M_RANK_AT_LEAST_TWO_SEPARATE_REVIEW.md) §7による。
全M可逆、全p、全normal frame、全higherと実正逆を戻す。
この旧入口でm≥3はtau零と正確に同値で、
semisimple／nilpotent二つのwhole first symbolを残す。
旧m二の除外tau零を点だけとして使ったのではない。
全式は[今回の有限入口](PLANE_COHIGGS_FULLRANK_TAU_ZERO_TWO_NORMAL_ENTRANCE.md) §1に記した。

## 2. 既全kernelを六shear＋二方向へ完全分解

元の全30／9、全350式を保持する。
保存済み両D1核は八、actual shear rank六、D2はrank九。
いずれも先行別作者が全行・complete kernelと同じnonzero minorを確認済み。
旧Sの同じ六minorのinverseから全30×8の

~~~
K=S A+N1 b1+N2 b2
~~~

をliteralに確認し、全A、b1、b2をreportへ保存した。
N1,N2は旧核の列を既存shear方向から引いただけであり、
新M1、rank、minor又はpointを作っていない。
N1,N2はその同じ六rowで零、二pivotで独立なのでshearとの直和である。

任意P1=K alphaを全六係数のactual finite inverse shear
x_i -> x_i+s_i z+s_{i+3}wで移せば、
P1=xi N1+eta N2となる。
逆は全sの符号反転、det一。全六sを独立変数として正逆、
保存SとのEuler商等式と、全xi/etaの150係数push/inverse roundtripを確認した。
このfinite作用はP2も同時に動かすので、P2の非斉次particularを落としていない。

## 3. Semisimple tau零の全二normal模型

元のweight零tensorは

~~~
P0=∂b∧((bz+cw)∂z+(az−bw)∂w)+zw∂z∧∂w.
~~~

保存順八kernel座標alpha0,…,alpha7から
xi=alpha0+alpha4、eta=−alpha3+alpha7。
完全な残方向は

~~~
N1=z²∂b∧∂w,   N2=w²∂b∧∂z.
~~~

両Nは同じconstant方向∂bを持ち、係数はbに依存しない。
従ってJac(xi N1+eta N2)=0はstrictに全xi,etaで成立する。
元weight二方程式はD2(P2)=0、旧D2単射性から全九P2は零。
全tensorP0+xi N1+eta N2は∂b、∂z、∂wの三constant方向だけに支持され、
skew rankは高々二。全xi/eta、全zero境界を含む。

## 4. Nilpotent tau零の全二normal模型

元のweight零tensorは

~~~
P0=∂a∧((bz+cw)∂z+(az−bw)∂w)−z²∂z∧∂w/2.
~~~

残normal係数は
xi=−2alpha2/3+alpha3/3−4alpha5/3、eta=−2alpha6+alpha7。
完全な二方向は

~~~
N1=∂a∧(bw∂b+cw∂c+3zw∂z),   N2=w²∂a∧∂z.
~~~

共通constant方向∂aとa非依存係数から、ここでも全strict Jac(N)=0。
同じ旧D2単射性で全P2零である。
Eを全radial fieldとすれば、全xi/etaのliteral Euler同一性は

~~~
P0+xi N1+eta N2 + E∧(xi w∂a)
 = ∂a∧((bz+cw+2xi zw+eta w²)∂z
          +(az−bw−xi w²)∂w) − z²∂z∧∂w/2.
~~~

右辺は∂a、∂z、∂wの三方向だけに支持されるため、
全projective rankは高々二。rawの単一normal括弧だけを見たrank判定ではない。
Euler項はlinear fieldとの外積で、projective tensor及び実signed Pfaffianを変えない。

§3–4と全finite shear逆、全first-symbol normalizer逆を戻して冒頭の全源定理を得る。
元の全m可換conormalとm≥3 iff tau零を用いたため、
任意normal multiplicity、zero higher、irregular sourceを失っていない。

## 5. 実行した有限batchと保存範囲

新しい単独wrapperは増やさず、distinct job内にinline caller全文を固定した。
[plane_cohiggs_fullrank_tau_zero_two_normal](jobs/plane_cohiggs_fullrank_tau_zero_two_normal.json)
のgate PASS後、whole二源を一回だけ実行した。
予算180秒／3GiB、実測1.241912秒、最大RSS14,508KiB。
旧独立Fraction sparse primitivesと直前の凍結callerのbasis／encode helperだけを再利用した。
producer／CASをimportせず、新Jacobi operator、rank、minor、anchor、GB又はjetは作っていない。

[結果report](../../reports/CLASSIFICATION_PLANE_COHIGGS_FULLRANK_TAU_ZERO_TWO_NORMAL_20261108.json)は、
両全K=SA+N1b1+N2b2、旧D2 inverseから得るQ20=Q11=Q02=0、
全二normalの元350Jacobi恒等零、actual Euler-Pf全126恒等零、
全finite shearとその全global逆、実行inline sourceと全入力hashを保存する。
Pf零の有限確認だけでなく、§3–4の共通三方向／Euler恒等式が全rankの手証明を与える。

全rank M三のm二にある他の有限壁、rank M二M(p)非零、rank M一以下、
O(0)、scalar、zero-symbolその他のE源はこの命題から閉じない。
1054原稿へ新結果を混ぜていない。
