# Ga9の二つの全源を含む平面H+quarticの24次元成分

2026-11-07、PRIMARY E。producerの全計算PASS、独立再現と統合担当の受理は確認中。
下界・台帳・mainはこのノートで更新しない。
対象はGa9 `c=0` のpower1/power3全源で、linear-w源は未解決として保持する。
型6の全plane quarticの尽くしも主張しない。

## 全六磁気項を持つ二つの実Poisson族

baseを `(h,y,w)`、fiberを `(x,z)`、aを族のparameterとする。

    F1(a)=(3hw+2y²)²−16y(h+ay)³,
    F3(a)=(hw+6y²)²−144y³(h+ay).

平面五次は `D_i=h F_i` である。
全18係数を自由にしたquadratic gradient syzygyの線形系から、
kernelが正確に二次元であることを計算した。
その全basisからEuler部分を引いて `U h=V h=0` となるquadratic derivation二本を得る。
a=0のbase anchorが既存Ga9のcanonical anchorと一致するよう、定数GL2 frameを選んだ。
det(U,V)/Fは順に `12/(2a³−3)`、`252/(18a−7)` で非零。
従って `det(E,U,V)=unit*hF` であり、Saito criterionによりDはexponents(2,2)のfree quinticである。
τ=12だけからfreenessを推測していない。

`[U,V]=αU+βV` の係数α,βもbaseのlinear formとして全parameterで求めた。
任意のquadratic magnetic form

    R=r0 h²+r1 hy+r2 hw+r3 y²+r4 yw+r5 w²

に対し

    P_i(a,R)=∂x∧U+∂z∧V+(αx+βz+R)∂x∧∂z

はstrict Jacobiの全10式をQQ(a)[r0,…,r5,h,x,y,z,w]上で恒等的に満たす。
全六rを残しており、磁気切片の一点だけを変形したものではない。
両族はa=0で正則で、公開Ga9 normalformのpower1/power3 canonical二ベクトルと
**全150係数で**一致するrの値を保存した。
その固定coreではR自体が任意のままなので、磁気零も非零も族に含まれる。

公開canonical点は既済の全十一parameter源の密な23軌道を代表する。
各全源の飽和像上界23とその実軌道23から、全源飽和閉包はその軌道閉包である。
従ってcanonical点をPGL安定閉成分へ入れれば、逆正規化を戻した全parameter源も入る。
磁気零は非零の軌道と同一とは主張せず、その閉包への包含を使う。

## quarticの絶対既約性と一般配置

wに関するdiscriminantは

    disc_w F1=576h²y(h+ay)³,
    disc_w F3=576h²y³(h+ay).

`1−4a³`、`1−4a` がそれぞれ非零ならwの係数はprimitiveである。
a≠0ではyとh+ayの奇数冪により、a=0でもそれぞれh5y、h3y3により、
C(h,y)上のdiscriminantは平方でない。Gaussの補題からFはC上で既約である。
従って一般Dは単純Hと既約quarticからなる被約因子である。

F3には次の有理正規化がある。`Q=s0²−a s1²` として

    [h:y:w]=[Q²:Q s1²:12s0 s1³−6s1⁴].

a≠0、1−4a≠0では、s0/s1=±√aが同じ点 `[0:0:1]` へ写る。
二枝は滑らかで、hをyで展開した二次係数が `±12√a−6` と異なるのでA3 tacnodeとなる。
s1=0ではvanishing order2,3のA2 cuspを持つ。
この二点のdeltaは2+1=3で、既約quarticの算術種数3を使い切るから他の特異点はない。
H=hはtacnodeの両枝に接し、交点はこの一点だけ、全交点重複度4である。
従って一般元の配置は「A3+A2有理quarticと、tacnodeの共通接線H」で記述できる。

## 二族の全磁気parameterの実PGL移送

F1のparameter a≠0からF3のparameter b=a³へ移すbase座標変更は

    H=h,
    Y=−(h+ay)/a³,
    W=9w/a4−12y/a5−6h/a6.

行列式は−9/a6で、厳密に

    F3(a³)(H,Y,W)=9 F1(a)(h,y,w)/a8

となる。このbase変更後の二本のquadratic anchorはF3の全Saito basisの定数線形結合なので、
その2×2係数行列Tを全quadratic係数から一意に求めた。
fiberを `T^{-T}` で変更すると全anchorが一致し、実PGL pushforwardした全150係数が
F3のPoisson族に一致することを直接確認した。
六磁気係数の線形変換も明示し、その行列式は

    117649 a12(2a³−3)6 / [9(18a³−7)6]

で非零である。Pfaffianの一致だけからPoisson帰属を推測していない。

Saito chartの分母 `2a³−3`、`18a−7`、quarticのprimitive unit、
transferに必要なaと `18a³−7` を全て記録した。
transferをa=0で割ってはいない。両実族は別々にa=0まで正則で、generic transferでの包含を
閉包に延長することにより、両Ga9全源を同じF3源閉包へ入れる。
他の除外parameterはこのchart又は既約性の補集合であり、型6全支持の排除には使わない。

## 全24源閉包の成分性：真正49切片の非零二次で十分

F3の固定source点 `a=1, (r0,…,r5)=(1,2,3,5,7,11)` 一点で、
全150変数・全350式のD、Euler鉛直・GL5・scalarのM、
さらに**全七source parameter**の微分を加えたAを計算した。

    rank D=100,  rank M=48,  rank A=49,  D A=0.

非零小行列を全て保存した。このsourceの射影像は24次元以上である。
同じ点のh-Casimir cofactorは冪階数2,0で、hはdual cofactor像の外にある。
これらは一般parameterで成り立つ開条件なので、一般点は型6 nil22-image・base括弧零となる。
既済 `HYPERPLANE_QUARTIC_NIL22_COTANGENT_G_BOUND.md` の全支持上界24が適用できる。
従って既約な全source/PGL飽和閉包Sの次元は正確に24である。

Aの独立な49方向に非零となる49座標をlinear cutとする。
全D核は50次元で、このcut上の核は正確に一つのRのみ。
`A49⊕R⊕U100` が全150空間の真正直和になることを確認した。
この101次元切片でDの独立100式を形式陰関数で消去すると、残る変数は一つ。
Rの全二次Jacobi余核は非零で、`[D|Jac(R)]` の階数は101、非零101小行列を保存した。
従って残式は非零二次係数を持ち、被約局所切片は原点だけである。
高次はunitであり、完全IFTや追加三次・nil指数を求める必要はない。

49本のlinear cutによるKrull次元不等式から、全局所raw次元は49以下、射影次元は24以下。
Sを含む全Poisson成分はこの点を通り、次元は24に等しくなければならない。
従ってS自体が全Poisson空間の24次元既約成分である。
49方向には群作用以外のsource方向も入るため、局所唯一性やfullschemeの滑らかさは主張しない。

## 凍結したprior102全成分との相違

比較台帳はcommit `3fdd4295f818338f23079fa6201bfae6f291e390`、
COMPONENTS blob `f588ae1d60364ca771f746a46cc23fa5e11e6cf9` を固定し、全102行をreportにも保存した。
rank2又は次元24以外の成分とは閉rank又は全成分の次元で区別する。
残る同rank・同次元はMPS8件と既存三件だけである。

本sourceの実五次は全parameterで(h,y,w)だけに依存する。
方向微分写像 `C5→Sym4(C5*)` の階数≤3は閉determinantal条件で、
実Pfaffianが零となる極限も含め、S全体がこの条件に入る。
一方T02,T03,T16,T17,T21,T36,T37,T38のaccepted全sourceで、
全edgeを零にした**実coordinate-SNC**因子 `κ(p)x0x1x2x3x4` を全periodのまま直接計算した。
八件全てでκは非零有理函数で、方向微分の非零5×5小行列 `κ5` を保存した。
各named MPS全成分がこの実source点を含むという既済のsource認識を用いる。
SとそのMPS成分が同じならこのrank5点も閉rank≤3条件に入ってしまうので矛盾する。
MPS23の因子監査の流用や、未証明のgeneric quintic因数分解を用いていない。

残る既存三24成分は、閉 `1+2+2`、閉非被約五次、閉 `1+1+1+2` locusにそれぞれ全体が入る。
これらは本源の一般の単純H×既約quarticを含まない。
従ってprior102のいずれとも異なる。

## 再現と射程

runnerは `code/classification_ga9_cuspidal_quartic.py`、
有限契約は `jobs/ga9_cuspidal_quartic_component.json`、
公開結果は `CLASSIFICATION_GA9_CUSPIDAL_QUARTIC_COMPONENT_20261107.json`。
private work/jet/IFT入力に依存せず、Ga9公開normalform、既済次元上界、
MPS公開source八件と凍結Git台帳だけから再計算する。
全producerは13.084秒でPASSし、既存共有quartic runnerも変更していない。
独立reviewと統合担当の受理後にのみ、新成分の台帳登録を行う。
