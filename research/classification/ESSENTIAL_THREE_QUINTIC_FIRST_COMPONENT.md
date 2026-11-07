# 最小free quinticから得る本質三変数の23次元Poisson成分

2026-11-07、PRIMARY E。具体的な平面五次

    F=x⁵+x²y³+y⁴z

を扱う。全二次gradient核、全GL2 frame・六磁気源、全150/350周囲空間、
真正一法線切片の非零二次、全prior103との相違を証明する。
独立レビューは別のreceiptで区別し、主台帳の採用・成分数更新は統合担当だけが行う。
他のfree quinticの具体的尽くしや、全Poisson分類の完成は主張しない。

## 1. 一次資料と絶対既約性

このFはA. Dimca–G. Sticlaru,
*Free divisors and rational cuspidal plane curves*, arXiv:1504.01242v4,
Introduction (1.1), p.1 の d=5,a=1,b=c=0 に一致する。
当該primary PDFの式を直接読んだ。
同論文Theorem2.8(i), pp.6–7は既約free quinticが有理cuspidal・全weighted homogeneousで、
exponents(2,2)、τ=12であることも明記する。
本稿ではfree条件を文献の形だけから推測せず、全十八係数核で新たに確認する。

Fはzに関して一次、leading coefficientはy⁴、constant termのy=0への制限はx⁵である。
従ってC[x,y]上primitiveであり、C(x,y)[z]上一次の既約式である。
Gaussの補題からC上で絶対既約である。
三偏微分は独立で、Fの本質変数はちょうど三個である。

## 2. 全十八係数核と全Poisson源

baseを(x,y,z)、vertexを(u,v)とする。
三つの場成分それぞれに全六quadratic monomial

    x²,xy,xz,y²,yz,z²

を置き、全degree6 monomial二十八行に対するgradient syzygy条件を構成した。
全28×18行列の階数は16、核は次の二本で自由に張られる。

    U=(0,−y²/4,3x²/4+yz),
    V=(−9y²/80,3x²/16−yz/4,9xy/40+z²).

直接、U(F)=V(F)=0で

    det(E_base,U,V)=−9F/64.

従ってEulerとU,VはSaito基底で、Fはfree exponents(2,2)である。
全28行・二核vector・determinantをreportへ保存した。

α=divV=7z/4、β=−divU=−y/2 とおくと、任意の全GL2 frame(U',V')と
任意の六係数quadratic Rに対する全源は

    Π=U'∧∂u+V'∧∂v+
      (−div(V')u+div(U')v+R)∂u∧∂v.

全十strict Jacobi、全五divergence、全五Casimir contractionを、
GL2の四係数と六磁気係数を全て残した多項式として零と確認した。
全五PfaffianとEulerによる実五次を保持し、GL2 frameをGと書くと

    D=9 det(G)F/64.

となる。detG≠0だけがこの固定非零Dのrank四源に必要な局所化である。
磁気零・非零はいずれも保持する。

`ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md` のcanonical全源定理から、
この表示は当該Dを持つ全Poisson構造を覆う。
`ESSENTIAL_THREE_QUINTIC_MAGNETIC_REMOVAL.md` とその独立reviewにより、
全六Rを一意の実unipotent shearで零にでき、全frameはvertex GL2で移る。
従って固定Fの全源とPGL5飽和は一つのactual23次元軌道とその閉包に入る。
ここで磁気零をgeneric magnetic点と同一視した理由は、実可逆shearの完全な証明である。

## 3. 整数canonical点と全周囲空間

standard frame、R=0で全tensorを80倍した整数点を使う。
非零括弧は

    {x,v}=−9y²,
    {y,u}=−20y²,
    {y,v}=15x²−20yz,
    {z,u}=60x²+80yz,
    {z,v}=18xy+80z²,
    {u,v}=−140uz+40vy.

この点の実Dは900Fである。
全二次五変数二ベクトル150係数と、元projective Jacobiの全350係数を使い、
全微分D、全Euler25列・GL5二十五列・tensor scale一列のMを計算した。

    rank D=101,   dim ker D=49,
    rank M=48,    D M=0.

全D sparse行、全51作用列、階数の非零小行列を保存した。
全源は既にactual orbit23なので、raw source下界はvertical24＋projective scale1を戻して48である。

C*作用の座標重みは(0,0,0,1,1)で、tensor characterは+1。
g_t pushforwardをtで割る、実際のprojective stabilizer作用を用いる。
二ベクトルのrelative input weightsは−3,…,1であり、
全350×150微分がそれぞれのoutput weightへ交錯することを全非零成分で確認した。

| relative weight | ker D | actual action | true normal |
|---:|---:|---:|---:|
|−3|0|0|0|
|−2|6|6|0|
|−1|19|18|1|
|0|18|18|0|
|1|6|6|0|

これは固有値だけを見た仮のnormalではない。
weight-compatibleな実作用48方向に非零な48座標をlinear cutとし、
切片上の唯一normalをR0、Dが独立な101標準方向をU101と取った。

    C150 = actual48 ⊕ C R0 ⊕ U101

の全150直和determinantを非零と確認した。
R0のrelative characterは−1である。
選択した101Jacobi式と101標準方向のblockおよび全exact inverseを保存した。

## 4. 真正一法線の非零二次で成分性が閉じる

上の48linear cut上で、選択101式の形式陰関数定理により、
101implicit座標は一つのnormal parameter nの形式級数として一意に決まる。
R0の全二次JacobiをQとする。
選択101行を消す一意の二次補正T2を求め、全350残差D T2+Qを保存した。

全150列でλD=0となる十一項の疎functional λがあり、

    λQ=−57897/220448 ≠0.

である。λと全Q・T2・四つの非零残差を公開reportに保存した。
従って真の一変数IFT残式には非零二次初項があり、被約切片の局所支持は原点だけである。
全implicit graph、三次、高次や正確な非被約局所環を求める必要はない。
実行も二次で停止した。

48linear cutのKrull次元不等式から全局所raw次元は48以下であり、
actual source48の下界と一致する。
よって全23次元軌道閉包は全Poisson空間の既約成分である。
PGL5は連結なので各全空間成分を保ち、この点を通る任意の成分はその23軌道を含む。
局所次元上界と成分の極大性から、その被約成分はこの軌道閉包と一致する。
非被約schemeの滑らかさは主張しない。

normal weight−1とinput最小weight−3から、IFT graphの事前次数上界三もあるが、
この有限modelを計算する前に必要十分な二次支持結論を得た。
事前上界を理由に不要な高次計算を実行してはいない。

## 5. immutable prior103全成分との相違

比較入力は `inputs/essential_three_quintic_prior103_20261107.json`。
commit `228039f9f7e13ec850d0a49b0d21636fddbab814`、COMPONENTS blob
`e04cec24c8aed04e88053a57bb4a17cb8c2c1268` の103行から作った最小不変projectionである。
可変な現行台帳のhashに依存せず、全103のrank・次元・受理済みgeneric因子型を保持した。

異rank又は異次元の成分は自動的に異なる。
同じ23次元のMPS二十九成分について、受理済み完成sourceの全periodを残し、全edgeを零にした
実際のSNC点の五次を再構成した。全てκ x0x1x2x3x4、κ≠0であり、
五偏微分の非零5×5minor κ⁵と、各accepted sourceへの所属を保存した。

新軌道閉包全体では実Dの偏微分のspanは高々三次元である。
これは五次の係数に関する閉determinantal条件で、全PGL5作用・全退化で保たれる。
各既知MPS成分は上の本質rank五の実点を含むから、同じ成分にはなれない。
一つの新しい代表点がMPSに一般的でないというだけの議論ではない。

既存の二つの同次元既約五次成分

    R4_GA9_DISTINCT_IRREDUCIBLE_QUINTIC_23,
    R4_NIL22_LONG_IRREDUCIBLE_QUINTIC_23

についても、受理済みcanonical実tensorからEuler/Pf五次を新たに再構成した。
双方で偏微分rank五のexact minorを保存し、同じ閉不変量で区別した。
その実tensorの成分所属は既存の全source認識定理に依存する。

残る同次元非MPS成分は受理済み一般実Dが可約又は非被約であり、
各所定因子次数の積写像のproper closed像に全成分が含まれる。
今回の絶対既約Fはその像に含まれず、同じ成分にならない。
全103行の理由を一つずつreportへ明示した。
異なる成分の境界同士が交わらないとは主張しない。

## 再現と保存

再利用focused runnerは `code/classification_plane_quintic_saito.py`。
四つの有限契約・reportは `essential_three_quintic_first_` と
`CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_FIRST_` のprefixを持ち、順に
SAITO、AMBIENT、COMPONENT_SUPPORT、PRIOR103である。
各reportは実行したproducer bytesとSHA、全exact入力を保持する。
mode追加のたびに既存functionはmain以外AST不変である。
以前のshared quartic runnerと旧不変reportは変更していない。

算術実行時間は順に約0.184秒、4.313秒、1.125秒、0.431秒。
各実行の前にjob gateを通した。

これらは一つの固定Fの全源に関する証明であり、他のfree quintic代表や有限リストの尽くしは残す。
文献の全rational cuspidal quintic分類とfree条件の全parameter補集合を確認することが次のE義務である。
