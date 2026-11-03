# 独立review：二次polarizationによる全affine代数の既知五点分離

2026年11月3日。`S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63` を独立検査し、下記の範囲で **ACCEPTED** とする。

## 受理した大域的な結論

一般betaの係数体K、またはbetaだけの非空開で定義した係数環Rについて、元45q・全134式の**全actual affine代数**Aに、degree_q≤45の明示的なidempotent recipe eがある。

    A ≅ E_known5 × A/(e).

E_known5は先行したactual P29 finite-family quotientのfinite etale rank5代数である。
一般K上で五つのK有理点に分かれると仮定していない。
未知残余A/(e)は全scheme structureを保持する。
cyclic determinant、quintic H、未知qに依存する何かを反転した上だけの分離ではない。
反転するのは既知五点のJacobian normや分離元の判別式など**betaだけ**の非零式である。

元多項式環ではknown idealはI+(1−e)、未知残余idealはI+(e)で、二つはcomaximal、交わりはIそのものである。
未知残余の空性、有限性、reducedness、全五点尽くしは証明していない。

## 新規exact計算の独立監査

fresh producerを再実行した後、そのpolarization関数・determinant関数を使わずに次を確認した。

1. 元134式の全termを読み、16個のbeta係数を変更せず、q-degreeが0,1,2だけであることを検査した。
2. constant、linear、square、mixed quadraticの四つのuniversal templateについてSymPyの厳密QQ恒等式を証明した。
   全termはこの四類のいずれかで、beta-monomialは単なるscalarとして保たれるため、全134式・全betaでpolarizationが成り立つ。
3. frozen original134から別のHessian-array engineを作り、25組全ての134×45 midpoint Jacobianを再構成した。
   全134行についてJ_mid(q−p)=0、全25個の45次determinantがproducerと一致することを確認した。
4. 全five actual q座標から分離hを独立に評価し、Jacobian normとh discriminantを再計算した。

25個のmidpoint determinantsは対角以外全て0で、対角は

    19164,11398,1126,14420,13343  (mod32003)

である。積は16172で非零である。
h=Σ(j+1)q_jの五値は1484,12233,17329,20305,22439で相異なり、判別式は17166である。
これらのp単元性は、一般係数recipeで必要なbeta-only inverseが非空開に定義される証人である。
有限体の五標本から未知支持の空性を推論してはいない。

測定4.32秒はfresh producer replayと全独立計算を含む。
actual rank5 quotientそのものは、先行した有限family/Nakayama/rank下界の受理済み証明を使用し、今回その重い元計算を再実施したとは言わない。

## 全scheme分離の証明の監査

全Fiはq-degree≤2なので、任意のE係数普遍既知点pに対し正確に

    F(q)−F(p)=Jac_q F((q+p)/2)*(q−p)

となる。選択45 original行のmidpoint行列をMp、Dp=det Mp、dp=Dp(p)とする。
unit dpを反転すると、adjugateによって

    ep=(Dp/dp),  ep*(q_j−p_j)=0 in A tensor E

を得る。このidentityから任意のpolynomial fについてep*f(q)=f(p)*epとなり、f=epを代入してep²=epを得る。
全134式のidealへの所属であり、選択45式を全idealと仮定してはいない。

finite etale algebra Eを忠実平坦に分裂させるとepは五つのeiに分かれる。
分離hの値の差はunitであり、ei*ejはh(q)−h_iとh(q)−h_jの両方に殺されるため、i≠jで0である。
ei(p_j)=delta_ijも同じannihilationとdpの正規化から得る。
従ってe=Tr_E/R(ep)=Σeiは大域idempotentで、degree_q≤45である。

評価写像A→Eは先行したactual quotient写像である。
分裂後にf→Σf(p_i)eiで与えるsectionは降下し、**E→eA**はunit eについてunitalで、評価と互いにinverseである。
これをAへ含めた写像は1をeへ送るため、Aのunitについてunitalとは呼ばない。
このunit conventionを明記すれば、分離と二つのideal等式は全Noetherian scheme構造で成り立つ。
Aを有限次元やreducedと仮定する箇所はない。

## 保存したreview

* `code/s61_affine_polarization_splitter_review_20261103_v63.py`
* `reports/S61_AFFINE_POLARIZATION_SPLITTER_20261103_v63_REVIEW/report.json`

元のproducer、全証拠、先行actual quotientへの依存、proof noteのhashを確認した。
未知残余A/(e)が新しいexactな大域研究対象であり、その空性は引き続き未解決である。
