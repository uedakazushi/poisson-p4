# tangent conic・m二double root四全源：別作者算術と正しいfiltered fiber

2026-11-07、`/root/quartic_review`。**全四sourceとwhole moving parameterの独立算術はPASS。**
このreceiptだけの結論はhf/Z零の全像≤22、hf/ZNの既約源空、moving全t非零の像≤23。
**h=z³の粗上界27は未閉のまま保持する。** 後続のOre別証明を本receiptへ混ぜない。

## 入力と全sourceの独立再構成

public report `CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_20261107.json`（ebcafe36…）に固定された
producer6ff941e7…全文をinput identityに使う。mutable current runner/jobは参照しない。
既純証明 `PLANE_M2_DOUBLE_ROOT_SEPARATE_REVIEW.md`（cd47db9b…）は再実行・改変しない。
そのscope外だったconstant-model reductionは次節で新しく査読した。

新checker `classification_plane_m2_double_tangent_review.py` はproducer/CASをimportせず、
既独立Fraction・QQ(t)・五変数sparse演算から全式を再構成した。
finite job `jobs/plane_m2_double_tangent_separate_review.json` gate PASS後、
**4.862秒、28,660KiB**で次を完了した。

* 全150 P0、全350 Jacobi零とactualG2。三fixedとmoving全QQ(t)の表示式。
* 既忠実Euler商の全P1三十とP2九。全350×30/J1、全350×9/J2。
* Euler全25＋Levi全13＋scalar一の全150×39 matrix、全kernel、非零minor。
  Euler rank24を確認し、scalar GL5中心を二重に数えない。
* 全complete K1/K2及びnormalizer kernelの元全行消滅と独立性、全rank minor。
* 全39変数の全quadratic polarizations。780可能pair中276の非零Jacobi/actualD係数を全て確認し、
  保存表に無いpairも零を確認。全39actualD一次項も各sourceで再構成。
* movingの全三minorと全kernel分子・分母をQQ[t]のexact divisionでfactor確認。
  非零壁七つ以外のcomplex rootが無いこと、及び各壁の元全行からの再構成。

新report `CLASSIFICATION_PLANE_M2_DOUBLE_TANGENT_SEPARATE_REVIEW_20261107.json`
は実行reviewer全文、input hash、各whole modelのrankと判定を保持する。
高次compatibilityを勝手に捨ててPoisson solution sourceと呼んでいない。
この上界に不要な全quadratic式も保存している。

## 新constant-root正常形の必要完全性

producer `PLANE_M2_DOUBLE_ROOT_CONSTANT_NORMAL_FORMS.md`（fb33d82e…）§§2–3だけを今回適用する。
concurrent §§4–7の新正常形はこのreceiptのscope外である。

f=z(xz−y²)、N=2y∂x+z∂y。
既bad-row必要式はh=fでX=Y=−2N、h=z³でX=−8N/Y=−5Nを強制する。
h=fなら全off-formはclosed logなのでZ=cN。
今回はdet(E,X,Z)もZ(f)も零であり、actualG2=3f k20 u²がk20非零を正しく強制する。
相対normal scaleでk20一に戻した後、weights(x,y,z,u,v)=(2,1,0,0,1)のactual diagonalで
main P0のcharacterは−1、off cNのcharacterだけ−2。
projective scaleを合わせればc零とc非零→一の**二whole strata**となる。
triangleの場合の誤ったk20 unitをここへ類推していない。

h=z³の全closed-form五parameterから
Z=(Ax+By+Cz)∂x+(Dx−Ay+Ez)∂y。
独立checkerは全係数で

`det(E,X,Z)+Z(f)=−18Dxyz+18Ay²z+9Axz²+(9B−18E)yz²+9Cz³`

を展開した。
actualG2/u²=3f k20−この式がnonzero gamma z³に等しいことから
D=A=k20零、B=2E、C非零。
actual normal shear v→v+c uのDelta=−3NでEを消し、相対normal scaleでC一へ戻せる。
従ってh=z³も一fixed P0で尽くされる。全inverseは指定nonzero leading amplitude上で正当。
全higher sourceはこのactual変換で同時に動かし、その支持を固定点だけで代替しない。

以上の三constant模型と、既純証明のwhole moving t非零は全tangent normal-double leadingを覆う。

## actual filtered groupによる上界の橋

Sを保つPGL5のLeviは(GL3×GL2)/C*、dimension12、unipotent部は全六base-from-normal shears。
raw matrixGのrankをrとする。
Euler rank24、tensor scale一を除いたprojective Levi orbit dimensionはr−25なので、
projective P0 stabilizer dimensionは **s=12−(r−25)=37−r**。
単にLevi Lie kernel dimensionをPGL stabilizerと同一視していない。

全六filtered shearsはP0をraw tensorとして固定する必要はない。
associated-graded P0を固定し、変換後の全higher係数は同じ忠実30/9 Euler商へ正当に戻る。
Levi projective stabilizerも必要なtensor scaleを合わせて同じP0 sourceへ作用する。
G×fullhigher-sourceからprojective Poisson空間への写像では、
このstabilizerと六shearsの**group graph**が各fiberへ入る。
第一PGL因子で区別されるので、higher source上の作用に固定点があってもfiber dimension s+6を失わない。
Leviとunipotentの独立性により重複して引いていない。

P1は全kerJ1、固定P1上のP2は空又はkerJ2のaffine torsorの部分集合。
従って全nonlinear solution dimensionはk1+k2以下であり、model parameter dimension pを戻して

`image dimension <= 24−(s+6)+p+k1+k2`

となる。これはactual fullsource mapのfiber上界であり、chosen affine sliceの粗和ではない。

## 全結果と全parameter壁

| 全模型 | k1 | k2 | s+6 | p | 像上界 | 判定 |
|---|---:|---:|---:|---:|---:|---|
| constant hf/Z零 |8|0|10|0|22|G非優越|
| constant hf/ZN |7|0|9|0|22|全既約源空|
| constant h=z³ |9|3|9|0|27|本receiptで未閉|
| moving generic全QQ(t) |6|0|9|1|22|G非優越|

hf/ZNの全complete K1はbad u-row v²Cを持たない。
全P0 u-rowはu倍、全P2にnormal columnは無いので、全solutionのu-rowがu倍。
従ってactualDにu因子があり、指定既約源は全空。
一点でのfactorization又はkernelの一部への制限ではない。

movingの全非零壁はexact factorizationで次の七つ。

| t | k1 | k2 | 像上界 |
|---|---:|---:|---:|
|6|6|0|21|
|−3/2|7|0|22|
|−5/2|6|0|21|
|−3|8|0|23|
|−9/2|6|2|23|
|−6|7|1|23|
|−9|7|0|22|

全壁でs+6九、壁のpは零。
全minor numerators/denominatorsとkernel分母はt零及びこの七rootで完全に尽きる。
t零はactualG2零でm三以上の補集合であり、normal order二の対象に入れない。
したがってwhole moving sourceの像は二十三以下。
seedT18帰属から未知higher帰属へ飛躍する推論は用いない。

## 残余

本receiptのscope内ではh=z³だけを未閉に保持する。
その27は実際の27成分の存在を意味しない。
whole tangent/fullnormaldouble/fullM2 G又はEの完了はここから主張しない。
後続Oreの誤ったstrict-radial案7ea7…は本証明の入力でなく、別V2の正しい全projective式を要する。
triangle/concurrent/nonreducedR/m三/他surfaceも別義務。
