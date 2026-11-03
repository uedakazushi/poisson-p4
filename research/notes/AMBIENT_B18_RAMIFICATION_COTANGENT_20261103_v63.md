# generic B18でのnormalizer ramification：旧Ac8/Ac7型は延長しない

2026年11月3日。前稿でnormalizerが有限となることを証明したgeneric B18点P0を使う。
個々のCA,CB,Nの**normalization mapのintrinsic ramification ideal**について、
全QQ計算でそのcotangent rankが7と判明した。
この結果は、QaのAc8/Ac7完成型をrank2境界へそのまま延長することを排除する。
actual Dの全localtypeを決めるものではない。

## 結果

|normalizer|source cone dim|ramification idealのcotangent codim|projective critical tangent dim|
|---|---:|---:|---:|
|CA|29|7|21|
|CB|29|7|21|
|N|28|7|20|

表の最後は**tangent dimension**であり、critical locusの実次元の断定ではない。
従来のQaでのAc8（CA,CB）ならramification idealのcotangent codimは9、
Ac7（N）なら8である。従ってgeneric B18でこれらの完成型は成立しない。

## 1. intrinsicなcritical idealを小matrixへ戻す

source coneは旗9変数と固定旗fiber係数20（Nでは19）のsmooth chartである。
固定fiberからactual targetへの係数matrix Cはinjectiveで、一定の非零minorを持つ。
全normalization mapのsource Jacobianのmaximal minorsによるidealは、
このfiber blockを消した後の9 flag Lie列Qの全9次minorsのidealに等しい。
これはnormalization mapのFitting ramification idealであり、target componentのtangentを代用するものではない。

P0ではQ0のrankは8で、そのkernelは一つである。
CA,Nのkernel旗方向はMA=E_xh−2E_hz、CBではMB=E_xh+(1/2)E_hzである。
Q0の8×8 unit blockを反転して残る列をSchur消去すれば、
critical idealは残る全entryに厳密に生成される。
その一次部分は

    L = R − Q0[:,pivot]*(Q0[pivotrows,pivot])^−1*R[pivotrows,:]

である。Rは固定fiber変化に対するQのkernel列の微分である。
sourceの旗方向はGL同変性でsmooth factorとなるため、固定fiber方向でこの一次rankを計算すれば全sourceでのrankになる。

## 2. CBで必須のHessian correction

CA,Nの固定fiber mapはlinearなので、RはLie_M(deltaP)をfiber Cで割ったものである。
CBはg=−q*fというratio sourceを使うので、fとqのmixed second derivativeを含める必要がある。

P0でLie_M(P0)=C*co、kernel source-fiber velocityは−coである。
CBではcoのbinary f_hx成分が5/2、ratio成分は0である。
従ってratio方向deltaq=1に対する追加termは、正確に

    +(5/2)*raw_CB(g=h*x)

となる。一般式ではLie_M(deltaP)−deltaC*co、または
Lie_M(deltaP)+Hessian(source)(deltaalpha,−co)である。
producerはこの全150係数の追加列を保存している。
単純なtangent fixed-space変化でCBをlinear化してはいない。

## 3. exact rank証明

各Q0の全150×9行列、全微分R、補正済み一次matrix Lを保存した。
Q0はrank8、Lは三族全てでrank7である。
各非零minorだけでなく、全列を選択basisへ戻す有理係数matrixも保存し、
全matrixのfactorizationを照合したため、rank上界もexactである。
有限体・浮動小数による推測ではない。

## 4. 旧pinch型の排除とodd idealへの制約

旧環Acをsmooth factorとともに書くと、normalizerはB=C[[u,s1,...,sc,p]]で、
Aはu,si,p²,p³,si*pで生成される。
従って

    Omega_(B/A) = B/(p,s1,...,sc) * dp

であり、ramification idealのcotangent codimはc+1である。
normalizationとこの相対微分moduleは完成環のisomorphismで保存される。
実計算の7はCA/CBの旧9、Nの旧8と異なるため、これらの型の延長は不可能である。

共通rank2底RでCA/CBのimageをR+J*pと書くと、trace-zero p'を使ったramification idealは

    (p',J) in the regular normalizer.

Rのcotangentのnormalizerへの像はcodim1で、p'はその外の独立方向である。
従ってJの一次部分はbranch discriminatorを法としてrank6であり、
R自体でも独立な一次odd coefficientsは高々7本である。
九つのregular odd generatorsを選んで旧Ac8型を得るshortcutはこの点では成立しない。
高次・非regularなodd idealの構造を実際に調べる必要がある。

## 証明書

* `code/ambient_B18_ramification_cotangent_20261103_v63.py`
* `reports/AMBIENT_B18_RAMIFICATION_COTANGENT_20261103_v63/report.json`

QQ計算1.88秒。全ramification idealや全odd idealを展開したという主張ではない。
actual Dのnilradical・associated primes・depthは未解決のままである。
