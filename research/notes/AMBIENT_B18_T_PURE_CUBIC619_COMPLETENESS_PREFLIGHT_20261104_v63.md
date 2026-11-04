# Pure r³ の619本判定：全 Schur 補正と完全性

## 結論の範囲

実 T は五生成825関係で完全に表示され、既証の三次初期関係により候補は
最後の target `e4*r³` 一つだけになっている。この dual は actual五生成
frame における literal な source55 の r³係数である。

本 packet は、619本の **補正した実定数 quadratic-kernel 関係**の pure
値が、この候補を判定する完全な test となることを証明する。新しい619値
自体、source3/5、gauge3、native X3、348行列は計算しない。既存cacheにない
正則商の微分を零と置くこともない。

## 固定実 Schur frame の r 軸

全次数の元768比較列を、定数27 column操作で pivot Xp27 と741列へ分け、
さらに定数84操作で実単位 wedge84 と657列へ分けた frame を使う。元の
自由32 quotient の最後の成分 pi5 は、literalな元source55である。

原点の pivot下五行 Bp および wedge-unit下五行 Ystar は零である。
実 source第一係数の直接検査により、r軸ではさらに

```
(Xp55)_1,r=0,       (Ystar55)_1,r=0
```

となる。固定 literal d3 は定数項零で、その一次係数は z0,...,z6だけに
依存し r係数を持たない。この最後の事実は、r二次成分が零という主張では
ない。実 d3にはr²成分があってよい。

620本の零-linear sourceを x_j とし、第一27 Schur係数を a_j とする。
既知の一次pivot係数を a1、第二wedge係数を beta とすると、全次数の
行列積を r³まで展開するだけで

```
R_j,3,pure = rawP55_j,3
             - sum_i Xp55_i,2,rr * a1_i,r
             - sum_s Ystar55_s,2,rr * beta1_s,r.
```

となる。Xp55,0=Xp55,1,r=0 により未知 a2,a3 は消える。Ystar55,0=
Ystar55,1,r=0 と wedge残余定数零により、未知 beta2 と U1も消える。
全て ordinary polynomial coefficientであり、Hessianの2倍規約はない。

## beta1,r と m1/Euler の消去

前の実 quadratic theorem に従い

```
C2_j = X_j2-Xp1*a1_j,
beta1_j,r = [ U0^(-1) kappa_linear(m0*C2_j) ]_r.
```

ここでは単に quotient class の式だけでなく、r一次係数として正確な
等式である。実 a1を引いた sourceは二次から始まり、m1が消える。Euler
contractionの非R線形差の一次 wedge係数は `d3_linear*h0` であるが、
d3_linearにr係数がないため、この差のr係数は厳密に零になる。

この式は元の形式函数で定義し、I二乗や348商のnormal formを微分しない。
同じ実単位frameを保つことが必要であり、UやYstarの一部だけを別の形式
基底へ置き換えてはならない。整合的な変更の下では全Schur恒等式が同じ
値を与え、liftの差は以下のcorrected boundaryに含まれる。

nativeのsource55は全次数でliteral零なのでrawP3にnative X3は現れない。
beta1には実 X2だけが必要で、適用するconstant sourceが零-originなので
既存 effective native `(L0+L1)S` をその証明範囲で使える。ただし必要pairの
supportがB11に入るかは、実adjointで別途確認する。

## Theta：quadratic のliftを一意に補正する pure値

37個の選択した実linear関係 R_i の r²最後成分を h_i と書く。296本の
`variable*R_i` の quadratic 初期形のspanは既知177次元空間全体である。
そのliftのpure cubic値は、variableがrなら h_i、zなら零である。

この296列のkernelは119本で尽くされ、その全pure cubic値は独立監査で
零と確認されている。従って、既知177空間上にwell-definedな線形函数
Thetaがある。本 builderは固定receiptからThetaを独立に再構築し、119値
および全177基底に対する値を確認する。

actual minimal d3列のpure cubicはすべて零である。全次数式
`-Ystar*U^(-1)*d3` において、三次になり得る項の両方が
`Ystar1,r=0` または `d3_1,r=0` を含むからである。Thetaの全168 d3 quadratic
値も零であると固定有理数係数で確認する。従ってd3を使うquadratic liftと
296列だけを使うliftのpure値は一致し、曖昧さはない。

Thetaを180次元quadratic空間全体へ任意に延長し、Theta_extとする。保存
builderは固定177基底の消去による一つの明示延長を作る。前のm1 theorem
の二次proxy

```
Qproxy_j = pi*C2_j-Ystar1*U0^(-1)*kappa_linear(m0*C2_j)
```

と真のq_jの差はactual d3 quadraticのspanに属し、Thetaがこれを零にする。
したがってpure補正にこのproxyを使える。

## 619本の補正値と完全性

620本のq_jの、既知177空間を除いたimageはrank1である。その定数kernel
の次元は619である。lambdaをその固定basisとして

```
F_j = R_j,3,pure-Theta_ext(Qproxy_j),
F_lambda = sum_j lambda_j*F_j
```

を定義する。lambdaの和は既知177空間に入るので、Theta_extの選択を変えた
差も消える。マイナス符号は、元のquadraticを消すliftを引くことによる。

完全性は元の825関係から直接従う。定数column操作でこれらを37本のlinear
関係、620本の零-linear関係、168本のd3関係とする。任意の形式係数結合が
次数三以上で始まるとする。

1. 次数一で37初期形が独立なので、その37係数の定数項は全て零。
2. 次数二を既知177空間で割ると、620係数の定数項は619 kernelに入る。
3. その619関係とd3のquadratic liftを引くと、残る37係数の一次項は
   296列の119 kernelに入る。
4. それ以外の三次初期項は、linear初期関係のquadratic倍、quadratic初期
   関係のlinear倍、119 lift、168 corrected d3 liftの和である。

pure dualは最後の既知項を全て零にする。よって未知のpure三次初期関係は
上記619補正値で尽くされる。単にspecialized kernelを取ったものでも、
一部のT sourceを捨てた testでもない。

**全619の実値が零**と証明できれば、このdualは全三次初期関係を零にし、
残るpure classが生きる。既知上界と合わせactual HF3=1となる。
**一つでも非零の実値**があれば最後の候補も消え、m³T=m⁴Tとなる。
Nakayamaによりm³T=0、既知HF(5,3,2)から全長T=10となる。未知微分を含む
形式式を恣意的に零とした結果にはこの判定を適用しない。

## 実 source-order gate

84個のYstar55,2,rrは既存108 source cacheから全て得られ、非零は13個。
Xp27の後23個はnativeなのでliteral零である。最初四個はそれぞれ

```
-gdiff2[567067], -gdiff2[567096], -gdiff2[567125], -gdiff2[567154]
```

で、元nodeはliteralにoldg-newgである。これらa1_rを使う620 indexは13個
だけである。builderはこの式とnodeIDを検査し、四つの未知値を保持する。
sourceのlabel名だけによる同一視はしない。

既存のcommon105 source3はCA F9零7、rA=sA=rの二変数cacheである。own
CA/CB57 source4もそれぞれown last2座標のcacheで、copied FULLINV15の
thirdは実YB7/8第三係数へcomposeできる。しかし一般のHadamard h2には
mixed zi*r² numerator3、h3にはmixed zi*r³ numerator4が必要で、planeだけ
のcacheからは従わない。cubic denominatorのregular quotient g2には
numerator5、g3にはnumerator6が一般に必要である。

これらは今回新しく求めない。実adjointを通した後に依存が消えるかを
確認し、残れば正確なmissing-source gateとして保存する。新しいphysical
source jet、gauge3、numerator5 jobの許可をこの証明から推論しない。

## 再現

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_T_pure_cubic619_completeness_preflight_20261104_v63.py
```

60秒・2GB・32MiB・150000参照の固定検証。実619第三値、HF3、全長、nil三乗
の新たな数値結論は本packetに含めない。
