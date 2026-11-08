# Triangle・m二double・H21：十whole壁は全て既約quintic空

2026-11-08。PRIMARY E。既受理の全有限壁入口と完全K1を使う。
新しいM1、点、壁、ambient profile、GB又はnormal jetを追加していない。
一回の全源実shear逆／FULL9P2／元350 batchと、その後の最弱の純係数帰結である。
新算術と全source逆の別作者レビューは別義務として残す。

## 1. 全source入口と実逆

`TRIANGLE_M2_DOUBLE_SEPARATE_REVIEW.md` は、HF全Casimir源、二つのH21 fixed源、
oneparameterの壁外を既に全空とした。残る十壁は

    r=−7,−1,−3/4,−4/7,−1/2,−1/3,1/7,1/2,1,3。

元の完全30basis、完全K1、全9P2、全350式を不変入力とした。
P0は座標(a,b,c,u,v)で

    Pab=ab, Pac=−ac, Pbc=bc,
    Pau=−2(r+1)au, Pav=−(2r+1)av,
    Pbu=(r+3)bu, Pbv=(r+1)bv,
    Pcu=(r−1)cu, Pcv=au+rcv, Puv=(3r+1)uv。

全六actual shearsはa,b,c各々へu,vの定数線形結合を足すglobal PGL mapで、
逆は六係数の符号を変えるだけである。
raw36をEuler6と忠実商30へ分ける固定36逆を保持した。
全十源で実shear像はrank六であり、完全K1を六shear＋残差へ分ける定数行列Cと
左右逆を保存した。全K1係数alphaとraw P2九を**同時に**inverse-pushして、
全weight一が保存残差だけになることを確認した。
残P2写像は各九変数について `qnew=qraw+shift(alpha)` のidentity-linear triangular mapである。
したがってglobal polynomial逆があり、P1零のchosen slice又は六接方向だけの議論ではない。

全源のnormalizerはsource係数を一つも反転しない。
元各壁のnonzero rational constants以外を分母に持つ逆も無い。

## 2. 正常化後も全係数を保持する

r非−1の九源でK1の残差は一つ、bad coefficientをmと書く。
その全P1は正確に

    Pau += m v²,
    Pac += 4m a v/(3r),
    Pbc += m b v/(3r)。

r=−1では残差二つであり、同じbad mに加えてnを保持する。
追加項は `Pab += −n av, Pbc += 3n cv, Pbu += n uv`。
これは元完全kernel column三であり、Pbc=a uという別の元kernel方向も
六shearとこの残差の定数全逆を通して回収される。元K1のどの列も捨てていない。

全P2は最初から全九つ、すなわち三base-pairにu²、uv、v²を許す。
特にr=−1/2のhomogeneous `Pbc=q8 v²` は全保持し、
他の全inhomogeneous particularも同じ既存J2最大minorで戻した。
元350の全残条件と全signed Pfaffianはreportへ保存した。

## 3. 最弱certificate：元quad(a,c,u,v)のv⁴係数

projective Jacobiは標準の `E wedge strictJac` conventionを用いる。
全350のrow279はquad(a,c,u,v)、monomial v⁴に対応する。
全十源について、その係数は

    J279 = 4m²/(3r)。                                      (1)

rは上の十値で非零なので、被約Poisson支持でm零が強制される。

この係数式は、実はP2 forcingを解く前の**全raw九P2**でも直接成り立つ。
base a=b=c=0及びu=0でstrictJac(a,c,u)のv³項を見ると、
唯一の寄与は

    Pua * partial_a Pac
      = (−m v²) * (4m v/(3r))。

E_vの外積符号が反転して(1)になる。
P2はbase-pair normal quadraticなのでa微分が零であり、u-rowも変更しない。
r=−1の追加nはこの三つ組の該当項を変更しない。
従ってunknown particular、q8又はnの零／非零を仮定していない。

保存された元350のrow54（quad(a,b,c,u)、monomial b v³）も同じsquareを与える。
一回の全源計算で再現した値は次のとおりである。

|r|mの保存名|J279=J54|
|---|---|---|
|−7|r0|−4r0²/21|
|−1|r1|−4r1²/3|
|−3/4|r0|−16r0²/9|
|−4/7|r0|−7r0²/3|
|−1/2|r0|−8r0²/3|
|−1/3|r0|−4r0²|
|1/7|r0|28r0²/3|
|1/2|r0|8r0²/3|
|1|r0|4r0²/3|
|3|r0|4r0²/9|

## 4. 全actual Dのu因子

m零の正常化後は、全P0、残るn、**全raw九P2**について、u-rowの各成分がuで割れる。
従ってsigned Euler-Pfaffianはu因子を持つ。
uをomitしたPfaffianにも外側にuが掛かり、それ以外の各Pfaffianはu-rowを含むためである。
u-rowを含むPfaffianだけを調べて最後の一つを落としていない。

actual shearはuを固定し、Euler gaugeはactual Dを変えない。
全constant C逆を戻しても、元完全源のactual Dはproper factorを持つか零である。
ゆえに十壁の全rank四・被約既約quintic源は空。
既受理のHF／二fixed／壁外空性と合わせ、
**reduced nonconcurrent triangle、normal位数二、generic double normal quadraticの全源は空**となる。

## 保存と誤読防止

ONE batchはgate PASS後180秒／3GiB、10.485523秒、RSS63668KiBで全算術PASS。
`CLASSIFICATION_PLANE_M2_TRIANGLE_H21_WHOLE_SOURCE_20261108.json` に
全C／逆、全alpha＋raw九の同時実逆、P2 shifts／forcing／全残条件、全150括弧、全Dを保存した。
最初の成功reportは、因子を無条件に得なかったr=−1,−1/2を当初の予定通り
`FULL_FAMILY_SAVED_STRUCTURE_PENDING` と記録している。そのfrozen bytesを直さない。
本稿は**そのreportに既にある元squareを適用する別の数学的帰結**であり、
必要線形supersource上で先に因子があると偽っていない。
非被約Poisson schemeにnilpotent mが残る可能性と、C上のpointwise空性も区別する。
他のR、squarefree normal quadratic、concurrent、nonreduced、normal位数三以上は別義務である。
