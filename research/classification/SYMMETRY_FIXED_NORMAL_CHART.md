# 一つの固定有理法線形式による一般支持判定

2026-11-06 / v64-classification-first。Gの構造的還元。全S61の一般既知五点尽くしや単位所属は未証明である。本稿の固定形式は**全W16の関数体上**で定義する。既存疎K5切片の同名法線や同切片だけでの所属判定へ自動的に読み替えない。

## 1. 環と固定frame

W16をprimitive一次部分の空間、K=C(W16)、R=K[q0,...,q44]、Iを元S61の全134式のidealとする。既知五点が張る一つのaffine3planeをPi_betaとする。全体の被約台がこのplaneへ入ることは仮定しない。

Pi_betaの有理点として既知T39のglobal graphのc(beta)を選べる。T01三次被覆のrootを選ぶ必要はない。U_betaをPi_betaの3次元方向空間、W_betaをその45×3有理basisとする。J=[0,1,2]行のminorが非零な底開で

    E_beta=W_beta*(W_beta[J,:])^(-1),   E_beta[J,:]=I3,
    n_j=q_(j+3)-c_(j+3)-E_beta[j+3,:]*(q_J-c_J), j=0,...,41

と定義する。全42個でPi_betaのideal Lを生成する。cをplane内の別の点に取り替えてもn_jは変わらない。

ここで固定形式を

    eta=n0+n34
       =q3+q37-c3-c37
        -(W_beta[3,:]+W_beta[37,:])*(W_beta[J,:])^(-1)*(q_J-c_J)

とする。係数はK、q次数は1である。

**frameの区別。** このJ=[0,1,2]は保存された有理T01中心でのframeである。実疎K5模型の`full134-direct-circuit`はpivot_q_rows=[3,6,13]、別のnormal_q_rows42を使う。従ってそのtyped模型の`n0+n34`とは異なる。以降の定理はこの全K上のetaについての定理である。

## 2. planeが群作用と整合する理由

既知五点はT01、T39、T40の三つのPGL5不変閉族の一般fiberから来る。これらの族がPGL5不変であることは既存の族の構成・成分対応への依存である。連結群が有限個の全空間既約成分を個別に保つことからも確認できる。

既存のuniversal3planeは五点を含む。五点のaffine spanは一般に3次元である。保存中心のplane座標(tau,rho,sigma)では、A=24010/13≠0として

    (0,0,0), (0,0,1), (0,1,-1/2), (sqrt(A),1,0), (-sqrt(A),1,0)

が既存五二次式の全五点である。最初の点から他の点への差は3次元を張る。相異なる五点とaffine-spanの非零minorは非空底開で保たれる。従って一般Pi_betaは既知五点の**affine spanそのもの**であり、任意に選んだ余分なplaneではない。群作用がfiber座標上affine-linearであることから、Pi_betaは群同変である。

この議論は全S61の未知点を五点へ含める証明ではなく、既知部分から定めたplaneの同変性の確認だけである。

## 3. 第二jet差への線形作用

Sp4の線形座標変換の後、保存S61 gaugeへ戻すperspective係数は一次jetだけに依存する。有限perspective変換の二次部分は

    Q'=Q-f L+y∧(L a)+f²J-2f*y∧(J a),  f=a^t y.

したがって同じ一次jetを持つ二つの入力の差deltaQは、この再正規化では変化しない。三次尾部は二次へ降りない。ゆえに固定gの全fiber作用は

    q -> A_g*q+b_g(beta),

で、45×45行列A_gはbetaに依存しない。

A_gは元150係数のweight1=0でのweight2変分から独立に構成した。affine Euler変換を通したdeltaQにbeta2=-J deltaQ Jを適用すると、60係数の閉quadratic2formを得る。保存したM60×45はこの全45次元閉空間のbasisで、deRham行列D16×60についてrank D=15、DM=0、左逆TについてTM=I45をQQで確認した。

採用した有限作用の規約は

    beta2(y) -> g^t beta2(g y) g.

これは物理的なbivector pushforwardではg^(-1)に対応するが、Sp4の元の集合は同じである。60係数上の作用はLambda²(g^t)とSym²(g)の有限tensor積で、A_gはその作用をMへ制限したものとなる。

## 4. 固定形式の軌道が全42法線を張る証拠

Wを保存中心での45×3方向basisとする。J行がI3なので、nの線形部分は(q_[3..44]-W_[3..44]*q_J)である。v=e3+e37とおけば、gによるetaの引戻しのq線形部分は

    a_g = v^t A_g
          -v^t A_g W * ((A_g W)[J,:])^(-1) * A_g[J,:].

必ずa_g W=0である。frameを再正規化した式そのものであり、qの定数項はplane内の中心差によって消える。従ってg etaはLのK線形結合となる。

`classification_symmetry_chart.py`は64個にあらかじめ固定した四transvection積だけを使う。全gは整数行列でg^t Jg=J、det g=1をQQで検査する。元のn0単独では、この有限集合でのrank下界37しか得られなかった。これを上界37や全軌道spanの決定とはしない。

続く一回だけの固定形式eta=n0+n34では、**同じgに対する二行の和**からa_gを組み立てた。二つのorbitの行の合併ではない。固定64個のうち45個までで42本の独立行が得られ、末尾42列の42角minorは

    det = 27110 != 0 (mod32003)

となった。全45行のa_g W=0、60係数でのM A_g=作用(M)、42角minorの逆行列恒等式を確認した。新しいbeta点・素数・group集合は試していない。実行時間は1.21秒である。

WのQQ値とその元betaへの束縛は、保存frame・元150係数・既存universal3planeのexact証明を用いる。各評価済み係数の分母は32003単元であり、このminorの非零像はQQでの同じminor値の非零を証明する。従って全betaの有理関数としてもminorは非零であり、K上でrank42となる。これは有限体の**非零minor**による証拠であって、有限体のunitや空性を標数零へ移したものではない。

## 5. 一つの固定chartによる定理

**定理。** 上述の全S61模型、実際の群作用、既知五点のcanonical affine3plane、正当な元座標との橋を依存とする。上の有限42minor証明書を用いると、全K上で

    L subset sqrt(I)
      iff eta belongs sqrt(I)
      iff 1 belongs (F0,...,F133,1-u*eta) in K[q0,...,q44,u].

**証明。** eta∈Lなので第一条件から第二条件は直ちに従う。逆にeta∈sqrt(I)とする。固定g∈Sp4はW16への可逆な線形作用を通してKの自己同型を定め、Rと元Jacobi ideal IにはK-semilinearな自己同型を与える。よってsqrt(I)も保たれ、各g etaはsqrt(I)に入る。

plane同変性からこれらは上のa_gによるnのK線形結合である。42minorがK上非零なので、全n_jは有限個のg etaのK線形結合に戻る。従ってL⊂sqrt(I)。第二条件と最後のunit criterionの同値はRabinowitschの局所化同型`R[u]/(I,1-u eta)=(R/I)[1/eta]`による。□

この定理にはBの精密構造や元H0・T2の計算は不要である。Sp4がgeneric幾何fiberの各成分を個別に保つという余分な仮定も要しない。必要なのは**全一般代数へのsemilinear自己同型**である。既存plane交叉E5の定理を合わせれば、一つの上記unit証明書が一般五点の被約尽くしに十分となる。ただしunit証明書はまだ得られていない。

## 6. 疎K5へは自動的に降りない

Kへの群×切片の忠実な体拡大では、この定理とetaの所属問題をそのまま引き戻せる。しかしetaは非不変な有理係数形式であり、引戻したetaには一般に群変数も残る。

したがって「群変数を恒等元にして、固定K5切片上でetaがidealに属することを証明した」という条件だけでは、全K上のeta所属を証明したことにならない。全S61の尽くしという結論が不変であることと、この一つの非不変chartの所属条件を群変数なしで検査できることは別である。

次のどちらかが必要である。

* 群と切片の関数体の全変数を保持して引戻したetaのunitを証明する。
* 適切な固定群移動・frame変更や、切片上の別の十分条件を実際の恒等式で証明してから群変数を除く。

その輸送の橋は未構築である。現在のK5上の超越一次形式chartは別の正確な入力還元として保持する。本稿の全K16固定形式を、直ちにより軽いK5計算へ置き換えられるとは主張しない。

## 7. 証明書と補集合

* `code/classification_symmetry_chart.py`：元quadratic閉空間と有限64g、n0の有限span。
* `reports/CLASSIFICATION_SYMMETRY_NORMAL_CHART_20261106.json`：M、D、左逆、QQ frame、全g、全source hash。n0の下界37を保持する。
* `reports/CLASSIFICATION_SYMMETRY_FIXED_SUM_20261106.json`：固定eta、選択42行、非零minor27110、同じ有限群、literal再現コードとそのhash。元producer bytesは変更していない。

底の除外元はcanonical planeの定義分母、J-minor、選んだgに対する(A_g W)[J,:]のminor、最後の42minorを含む。すべてbetaだけの有理式であり、分母払いして底の開へ移すときに補集合をEの義務として残す。未知qを反転するのは最後のRabinowitsch chartであり、その十分性は5節の定理で説明される。

今回、unit・GB・新しいjet・nil指数・全成分上界は計算していない。Gと分類主定理は引き続き未解決である。
