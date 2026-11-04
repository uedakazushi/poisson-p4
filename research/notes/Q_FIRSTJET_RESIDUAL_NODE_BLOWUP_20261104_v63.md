# 同じ係数点のproper (2,2)模型：商曲線の正規化は幾何的種数1

2026-11-04。元の平面四次smoothness gateが失敗した後、二つの普遍節点を解消する
**異なる模型**を一度だけ検査した。同じ係数点(1,2,3,5,7)、同じ素数32003を保つ。
新しいbeta/source点や不変式job、別係数へのretryはない。

## 1. 二つの有限etale節点を全係数で保持する

generic baseをK=QQ(x,y,a,k1,k2)とする。既証のm=0上の節点多項式P0(r)は二次で、
相異なる根alpha,betaを持つ。必要な係数分母、x,yとdiscriminantを可逆にした
Z_(32003)上の係数baseをSとする。monic P0の根alphaを付加したS'はfinite etale、
betaはtrace−alpha、delta=alpha−betaは単元となる。

同じ特殊係数では、その二次etale algebraは

```
Fp[alpha]/(2alpha²−287alpha+10914),   Fp=GF(32003),
beta=287/2−alpha,       delta=alpha−beta.
```

discriminant−4943は非零。rootがFpにsplitすることも、特定の平方根を選ぶことも
仮定しない。以下の算術は常にbasis(1,alpha)のpairで実行した。

## 2. Properな(2,2)模型と元の四次との関係

元のhomogeneous plane quarticをF(r,m,z0)とし、

```
U=r−alpha*z0, V=r−beta*z0, M=m,
r=(alpha*V−beta*U)/delta, z0=(V−U)/delta.
```

を使う。Fは両節点で二重に消えるため、

```
delta^4 F((alpha V−beta U)/delta, M, (V−U)/delta)
```

のU,V次数はそれぞれ高々2。M=1としたH(U,V)をP1×P1へ(2,2)で閉じる。
今回、保存した元の15係数を一つずつ変換し、alpha²の二次関係で全項をreductionした。
UまたはV次数が2を超える全係数が零であることを確認し、残る全係数pairを保存した。

F(r,0,z0)=lambda(r−alpha*z0)²(r−beta*z0)²で、
lambda=16x⁴y⁸は単元。従ってHのU²V²係数lambda*delta⁴も単元である。
このbihomogeneous (2,2)式はS'上のproper flat familyを与える。具体的には
unit leading monomial U1²V1²によるdivisionで、bigraded coordinate ringの各graded pieceが
S'-freeになる。fiber内で式が恒等零になることもない。

幾何的にはP2の二節点をblowupし、両点を結ぶ直線をcontractしたP1×P1模型である。
strict transformのclassは4H−2E_alpha−2E_beta、結ぶ直線のclassはH−E_alpha−E_beta。
交点数は4−2−2=0。結ぶ直線は元の四次の成分ではないので、strict transformと交わらず、
contractは曲線の近傍で同型である。従ってこのproper (2,2)曲線が滑らかなら、
それは元の平面四次の正規化を与える。

## 3. 同じ特殊係数で四つのpatchが全てPASS

Macaulay2 1.24.11のdefault GBを用い、各P1×P1 affine patchについて

```
( 2alpha²−287alpha+10914, H_patch,
                      partial_U H_patch, partial_V H_patch ) = (1)
```

をFp[alpha,U,V]で検査した。alphaの微分を加えない。alphaはetale定数baseであり、
curveのrelative smoothnessをU,Vの二偏微分で調べるからである。
四つのpatch全てでunit ideal。これにより両root選択を含む**全幾何的特殊fiber**が滑らか。
実測1.909995秒。scriptとlogはそのまま保存する。

S'上でnon-smooth locusは閉、familyはproperなので、そのbaseへの像も閉である。
この特殊fiberが全てsmoothだから、その像を除いたbase openは非空で、generic fiberを含む。
S'が二componentにsplitしていても、全root選択を検査したため各generic componentに適用できる。
従ってgeneric (2,2)曲線はgeometrically smooth。

さらにsmooth (2,2) divisorは幾何的にconnectedである。互いに交わらない異なるeffective
成分のclasses(a,b),(c,d)ならad+bc=0が必要だが、それらを合わせて(2,2)にすることはできない。
smoothでconnectedな曲線はgeometrically irreducible。adjunctionまたはP1×P1の式から
算術種数は(2−1)(2−1)=1なので、**generic商曲線の正規化種数は正確に1**。

元の平面四次もgeometrically irreducibleで、既知二節点のdelta計2が算術種数3から
種数1への全差を占める。従ってgeneric平面四次の特異点はその二節点だけである。

## 4. 残る二次coverと範囲

full ordered residual curveはh²=qR=r²+4fuvという二次coverである。
今回閉じたのはその商のproper normalizationであり、full coverの正確な種数・
選択されたactual componentのm-degree8には別の分岐判定を要する。
qRがK(m)上のfinite algebraの単元だったことを、projective curve全体で零やpoleを
持たないことと取り違えない。

保存先は `reports/Q_FIRSTJET_RESIDUAL_NODE_BLOWUP_20261104_v63/`、専用sourceは
`code/q_firstjet_residual_node_blowup_20261104_v63.py`。旧false smoothness report、
reciprocal packet、genus上界packet、二節点packetのbytesは変更しない。
再現する場合は既存環境をactivateし、新しいoutput先を指定する。
