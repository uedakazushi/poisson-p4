# B18 actual係数preflight：有限DAG・48+5単元系の独立監査

2026-11-04。**actual nu/I0/J0と両族individual cubic係数の有限構成を受理する。**
相対CB graphをCAへ評価する六membership liftと要求12pairは、このstageでは未構築。
後続packetの成果をこの凍結stageへ遡及して書き込まない。depthやnilradicalの判定も行わない。

## 検算方法とactual source binding

独立checkerは全41,315 scalar nodeを別のFraction interpreterで評価し、全105変数について
原点微分を疎なautomatic differentiationで追う。全inverse nodeの引数は非零origin値を持つ。
本体のc/v/a/m/n/i以外のoperationやfuture referenceは無く、formal parameterの逆元を混ぜない。

元の全150 source polynomialは独自のAST→sparse rational monomial parserで読む。
元sourceのconstant、各k一次、p一次・p二次係数を取り、physical150上のGL Lie作用を
index式から別に構成する。全150×19 orbit列と全576成分のM24を照合する。
さらに一回だけ有限DAG builderをignored workへ再実行し、全source rootsを含む
actual-source-coefficients.jsonがbyte一致することを確認する。
このreplayは新しいformal係数を解くsolverではなく、凍結した有限演算のbindingである。
旧ungauged R9 rank7 job、大きいGB、source jet、新点を再実行していない。

source pushforwardの規約も確認した。座標はG^-1で代入し、bivector基底はGでwedge変換、
scalar群係数を全体へ掛ける。block inverseのcofactor転置はこの規約と一致する。
二次basisの積はp²=sp+tから、normはa²+sab-tb²。
inverse pair `(a+sb,-b)/norm` とtrace変換 `(a+bs/2,b)` はこの関係の直接計算である。

## 共通physical baseと48変数Hensel系

実切断はweight0の群15と共通even base9
`27,108,109,112,113,117,118,119,61`、positive baseは別にrow140である。
γ15とsource k0,k1,k3,...,k8の23未知数をζ_j=a_j+b_j pとし、s,tを加えて48未知数。
source24式は上記24 physical rowsの**全て**から来ており、basis(1,p)への係数化が48式。
全原点値は零である。

独立DAG微分で、row順even24/odd24、column順(a23,t,b23,s)のJacobian全blockを照合する。
上左と下右は、元sourceから再構築した同じM24、上右は零。
下左cross項は零と仮定していない。p一次imageが零、p二次がWであるため、このblock構造は
原点の二次basis関係からも従う。

```
CA detM24=-9216/25, detJac48=84934656/625,
CB detM24=  576/25, detJac48=331776/625。
```

全inverse nodeはorigin-unitなので、必要なら補助inverse未知数を足してpolynomial systemにできる。
形式的IFTはR0上の一意なpointed algebraic power-series solutionを与える。
既受理のactual free quadratic normalizerの(1,p)係数もこの系と同じ原点条件を満たすので、
一意性でactual係数と一致する。抽象別模型のrootを選んでいるのではない。

## Positive5系と全次数w0..3

全sourceの各monomialを読み、weight0部がey∧V(X)、k2部が四座標X上のV2(X)で、
両者がyに依存しないことを確認する。zero群G0はy/Xを混ぜず、positive四方向はE_(i,y)。
従ってG*G0^-1=I+wM、M²=0であり、positive五未知数についての式は正確に線形。

checkerは保存pair matrixの全Gaussian行基本変形を**既存DAG内でsymbolicに再生**する。
各pivot node、row交換、quadratic inverse、最終solution nodeを全て一致させる。
五つのpivotのorigin値は両族とも1,-1,-1,-1,1で、全て単元。
五列のorigin値は独立に構成したphysical Lie/k2列とも照合する。
通常の消去法の恒等式とquadratic inverseのnorm identityにより、これは全次数での解であり、
数値的residualや有限jetだけの検査ではない。p、s、tそのものを反転していない。

実際のshearは
`(ey+w*l)∧V(X-w*l*y)+w*k2hat*V2(X-w*l*y)`。
V,V2は二次なのでw次数は正確に高々3。全150の四係数が有限DAGで保たれている。
高次sourceを切り捨てたjetではない。

## Own trace基底と未完の相対reader

各族で `tau_i=p_i-s_i/2,nu_i=t_i+s_i²/4`。
全出力rootのdependencyを調べ、CAにCB変数が混ざらず、CBにCA変数が混ざらないことも確認する。
nu root、positive odd rows46/95のJ0 root、八I0、六raw cubic rows9/39/54/84/99/144が
全150×4の保存出力と一致する。originでj=0でも係数函数が恒等零とは結論しない。

相対readerの共通1×10 matrixと六右辺nodeは、CB自身の保存出力へ個別にbindする。
zero oddのraw source rootとshear出力rootには`-((-1)*root)`という符号wrapperの違いがある。
checkerはconstant倍とnegationだけを正確に除き、同じrootであることを確認する。
単なるnode ID不一致や原点値の一致でmembershipを判断していない。
その60個のmembership係数hはこのstageで未解決であり、tau_B=tau_Aとはしていない。
individual cubicの二つのpairを単に引いても交叉idealのgraph差にはならない。
membership liftを得た後のtemplateは正しいが、liftの選択差を吸収するため全CB kernelと
lower graph relationsを残す必要がある。

受理したのは有限actual係数入力の構成であり、top T3への所属、四socle test、
全Dのnilradical、depth19/20、後続telescoping graph liftはこのreviewの結論に含めない。

再現（outputとworkは新しい保存先を指定）:

```
python code/ambient_B18_actual_top_coefficient_preflight_review_20261104_v63.py \
  --output /tmp/B18-actual-coefficient-review \
  --work /tmp/B18-actual-coefficient-review-work
```
