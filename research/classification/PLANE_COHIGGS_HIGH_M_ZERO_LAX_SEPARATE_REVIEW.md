# Co-Higgs rank三高mとrank二nonkernel零Lax：独立純証明査読

2026-11-08、quartic_review、PRIMARY E。**両全源定理を受理する。**
新operator、M1、rank、minor、point、profile又はCASを実行していない。
旧全kernelと全normalizerの独立認証を固定入力とし、今回の全源への論理を別作者として検査した。

## Rank M三、tau零の全higher

全mの可換conormal／Lax入口と旧全normalizerから、M可逆ではR=tau M(U)、
M(U)は非零でsemisimple／nilpotent二型へ尽くされる。
保存したleading actualDはtauを非零多項式に掛けたものなので、m三以上は正確にtau零。
これは旧m二で除外した一点をgeneric結果だけで復活した推論ではない。

両型の完全D1核八、actual shear六、D2単射は既別作者が全350行から認証した入力である。
今回のK=SA+N1b1+N2b2は全30×8等式であり、
同じ六rowでN1/N2が零、残二pivotで独立なので真の直和となる。
有限shearの全六係数と符号反転逆を使い、任意のK1を全二normalへ戻す。
任意のPoisson solutionを移したとき、同じ実作用でFULL九P2も移すので、
homogeneous P2だけを残す誤りはない。

semisimple残方向はz²∂b∧∂w、w²∂b∧∂z。
nil残方向は∂a∧(bw∂b+cw∂c+3zw∂z)、w²∂a∧∂z。
いずれも共通constant方向を持ち、その方向に係数が依存しないため、
二normalの任意一次結合の**strict** Jacobiは零。
linear項は既D1核条件、weight二はD2P2=0となり、単射性から全九P2が零となる。
従って全二normalの全値がactual Poissonであり、Poisson subsetだけの主張ではない。
有限shearを逆に施せば、元の全八kernelの任意の値がPoisson tensorへ延びることも得る。

semisimpleでは全tensorが∂b、∂z、∂wの三constant方向に支持される。
nilではEuler補正の符号を直接検査した：

```
P0+xi N1+eta N2+E∧(xi w∂a)
=∂a∧[(bz+cw+2xi zw+eta w²)∂z
      +(az−bw−xi w²)∂w]−z²∂z∧∂w/2.
```

右辺も三方向に支持されるので、いずれもprojective rankは二以下。
全xi/eta零、全六shear、全first-symbol inverseを戻す。
rawの一entry又はPf部分列だけのrank判定を使っていない。

## Rank M二、M(p)非零、T=R零の全四source

旧全Lax／normalizer入口で対象はNDSS0、NDN0、BSS0(0)、BN0(0)の四whole source。
次のsymbol familyを手で検査した。

|source|U|M_epsilon|全epsilonで固定のM(U)|
|---|---|---|---|
|NDSS0|∂b+∂c|epsilon aH+bF+cE|F+E|
|NDN0|∂b|epsilon aH+bF+cE|F|
|BSS0(0)|∂b|epsilon aF+bH+cE|H|
|BN0(0)|∂c|epsilon aF+bH+cE|E|

H=diag(1,−1)、E=E12、F=E21。
T=R=0、K=−q_(M(U))/2なので全Lax条件は全epsilonで成立する。
epsilon零で正確に元source、epsilon非零でrank M三、かつM(U)の非零ss/nil型は不変。
一般fiberの全源に前節を適用でき、必要なnormalizerがepsilon零で発散しても、
以下ではそのnormalizerを特殊化していないので支障はない。

D1/D2はこのpolynomial P0 familyからのpolynomial係数写像である。
一般fiberのrankは22／9。零fiber四型のrank22／9と既非零minorは旧独立認証済み。
旧minorが単位となるepsilon零を含む積分的open Bへ縮める。
全23minorが一般fiberで零なのでB全体で恒等零となり、
rankD1は全Bで22、kerD1はrank八vector bundle Kとなる。
D2の選択九row逆もB上正則である。

従ってQ_epsilon(s)=−D2_selected(epsilon)^(-1)Jac(s)_selectedはK上正則二次map。
全元Poisson solutionは、D2単射により必ずP2=Q_0(P1)を満たす。
全九P2を保ち、選択九式以外の条件を捨てていない。
一般fiberでは前節の**全K1**がactual Poissonへ延び、全Pfが零である。
ゆえにK上正則な全350Jacobi／全126actual Pf係数は、稠密なepsilon非零の全K上で零。
Kは積分base上のvector bundleなので積分的であり、全K上恒等零となる。
任意の零fiberの八係数にもこれが適用される。
これは元の任意solutionを正則延長できる全候補族を与えており、
一般Poisson subsetだけから閉性を誤用したものではない。

**結論。** rank M三の高multiplicity全源、およびrank M二M(p)非零でT=R=0の
四whole sourceは、全projective tensorがrank二以下。
指定されたrank四・被約既約quinticは全空。
BSS0/BN0のtau非零、他co-Higgs source又は全Eの完了を本receiptからは主張しない。

## 固定入力と認証範囲

二producer noteとその全manifest、旧high-m独立review、旧全rank三／rank二report、
finite shear補題、全二normal readback reportを本receiptのjob／manifestにhash固定する。
今回の認証は別作者による純証明の読解であり、新しい有限算術runの再現とは区別する。
producer readbackの実測値を今回の新計算時間として表示しない。
