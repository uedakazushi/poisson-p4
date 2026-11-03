# 全非零Qaでのactual CA/CB交叉：共通target translationによる完全なtransport

2026年11月3日。Dは既知の**被約成分**CAとCBの実際のscheme intersectionである。
個々の成分の完成環が同型という事実だけからDの同型を推論しない。
以下では同じactual target上の両全source mapsに可換する、明示的な共通automorphismを作る。

## 定理

全ての非零複素aについて、物理的Qaでのactual Dの完成局所環は、
既存のQ1での49変数・218関係式の全表示に、同じ七つのS61群因子と九つの値・perspective因子を加えた環に同型である。
全非零Qaで次が成り立つ。

* reduced quotientはNだけで、minimal primeは一つである。
* nilradicalは非零主ideal(F)、F²=0で、正確な冪零次数は2である。
* S61でnilmoduleはC[[u11]]*F、full projective chartではC[[u20]]*Fである。
* nilmoduleのannihilatorは滑らかな11次元（全射影では20次元）部分のprimeである。
* associated primesはNのminimal primeと、このembedded primeの正確に二つである。
* S61ではdim18、depth11、full projectiveではdim27、depth20である。

環の直和分裂は主張しない。a=0の物理的rank2点への延長や、全Poisson空間の全枝の分類でもない。
a=2/9でも、このactual CA/CB交叉の局所型自体はQ1と同じである。

## 1. universal actual-target Kにおける共通有限変換

既存の全critical-beta source slice

    tau=p,A=0,B=C=1,D=−p/2,F_base=s

を使う。CA,CBの源fiber座標は共通の順序で

    (r,k0,k1,k2,k3,k4,k5,k6,k7,k8)

である。
既存の十個のactual target q座標K=(v0,...,v9)を選ぶ全係数行列Lは、両族とも恒等的にdet L=−72である。

ここで**actual Kのv4をlambdaだけ増やす**ために必要なsource変換は、両族で全く同じで、正確に

    k0 -> k0+lambda,
    k2 -> k2−lambda*p,
    k5 -> k5+lambda*p²/3,
    他の全source座標は固定

である。これは単にk0だけを変える変換とは異なる。
対応するpotentialの差は

    lambda*(x²z−p*xz²+p²*z³/3)

である。

全sourceにこの有限変換を代入し、CAとCBのそれぞれの**全61出力**について厳密に

    u38 -> u38+lambda,
    その他の全60座標は固定

となることを確認した。全次数の多項式恒等式であり、一次微分だけの確認ではない。
inverseはlambdaを−lambdaに置いた同じ変換である。

従って共通actual targetのこのtranslationは両成分の全image idealsを同時に保つ。
全kernelであることが既に証明されたIA,IBについて、source translationとの可換性とinverseから

    translation(IA)=IA, translation(IB)=IB

となり、actual intersection ideal IA+IBも保たれる。
normalization同士のtensor productを使った代用ではない。

## 2. 既存218式とcolonでのさらに直接な検査

既存の49変数環は

    QQ[z0,...,z34,bb,cc,v0,...,v9,ss,aa]

で、aa=p²,bb=p³,cc=spである。
上のtranslationはこの環で単にv4 -> v4+lambdaである。
保存済みのCA109式とCB109式を全て検査したところ、**全218式は文字通りv4を含まない**。
従ってこのfinite presentationはv4について正確なpolynomial productである。

nilpotentはF=z25+v0+2aaで、固定される。
保存されたG=pFのactual target expressionもv4を含まず固定される。
exact colon ID:Fの45本のprime generatorにもv4は一つも現れない。
従ってF²=0、ID+(F)=IN、colon primeとその全module action、正確な二つのassociated primesは、
v4の全値に対して同じままtransportされる。
この検査に新しい大きなGroebner計算は必要なく、旧exact ideal receiptを保持して使う。

## 3. 全非零Qaへの適用

前稿で、a≠0の物理的Qaは正確にS61曲線

    u25=−4/3,u36=3,u38=9a/2, 他は0

へ移ることを確認した。
そのβは固定で、閉sourceではr=0,k0=9a/2、他の全kは0である。
一度だけ固定のSp4 diagonal torusを作用して、このβを本稿のB=C=1の基準βへ移す。
そのtorusはaに依存せず、monomialを定数倍するだけなので、この曲線は依然としてv4だけを動かす直線である。

本稿の共通translationはこの直線の任意の二点を結び、IA,IBの両方を保つ。
したがって既存のQ1のactual intersection完成環が、同じsource/image座標変換により全非零Qaへ移る。
七つのβ-orbit因子と九つの値・perspective因子の形式的productは既存の全slice定理に従う。

embedded layerは共通nonnormal層より一方向大きいことも保たれる。
S61のcolon自由座標はv1,v4,v5,v7で、p=s=0では

    v1=−k1,v4=k0,v5=2k2,v7=3k5.

従ってfree k0でのtranslationは、nilmoduleのsupport内の正則な方向を動かしている。

## 証明書と依存関係

* `code/q_CA_CB_actual_intersection_translation_20261103_v63.py`
* `reports/Q_CA_CB_ACTUAL_INTERSECTION_TRANSLATION_20261103_v63/report.json`

QQ上の新規検査0.36秒。
二族の全61有限変換、共通source速度、det−72、全218式と45 colon生成元のv4非依存性を確認した。
基礎となる全kernel・colon・nilpotenceは、凍結済み
`AMBIENT_INTERSECTION_EXACT_IDEAL_20261103_v63` とその独立reviewをそのまま使う。
原source/reportに変更はない。

物理的a=0の全完成型と、generic Poisson rank≤2の全域交叉境界は本稿の対象外である。
