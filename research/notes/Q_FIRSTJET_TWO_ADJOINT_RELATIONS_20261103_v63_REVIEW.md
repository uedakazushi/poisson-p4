# 二つの adjoint 共変量の関係と一般復元：独立監査

2026-11-03。対象 `Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63`。
受理済み Spin5 構造と、既存一点の独立監査を依存として明示する。
新しい first-jet 点、GB、22次以上の不変式評価は行わない。

## 1. 元の物理係数から共変量を再構成した

author module は import していない。元の24×16 primitive 基底から5×4の線形 Psi を作り、
C=Psi^T Q Psi J^-1 と V=Psi J^-1 Psi^T Q の全二次係数を Fraction 演算で再構成した。
五つの primitive 定数2形式から全10個の5次元作用を別に作り、V をこの作用へ戻して D を得た。
全25成分の戻し恒等式、CとDの全16×136係数、全 gamma trace を照合した。
従って V(D) は抽象的な記号ではなく、元の同じ Psi から作った V そのものである。

8次式の検査では author の9進指数符号を使わず、一変数4bitの16進指数符号を使った。
全次数8以下なので桁上がりはない。全二次係数の分母を4で払った整数多項式を掛け、
六つの cyclic trace words と六つの quartic trace products、合計12式の全係数を作った。
全式は共通に256倍されており、異なる尺度の式を比較していない。

七つの関係を全単項式で零と確認した。関係matrixの7角 minor は6で非零である。
12式の係数からの5角 QQ minor も保存値8388608と厳密に一致した。
従って全12式の rank はちょうど5、線形 relation 空間はちょうど7次元である。
五つの積関係を除く primitive relation は f3,f4 の二つで尽くされる。
ここで tr(C²D²) と tr(CDCD) を混同していない。
既受理の dim Inv8=5 により、12式の span は全8次 invariant 空間である。
これは任意次数の非線形 relation ideal の尽くしを意味しない。

## 2. 普遍的な復元恒等式と同じ一点の minors

定義から Gram=Psi^T Q Psi=CJ なので、行列積の結合則だけで

`V(D)Psi = Psi J^-1(Psi^T Q Psi) = Psi(J^-1 C J)`

が全 beta に対して成り立つ。gamma traceも普遍的に零である。
この二条件を未知20係数の Psi に関する24×20線形系 L(C,D) にした。
今回の code は各 matrix unit に作用させて全480成分を作り、author の index 公式と別に照合する。

numeric 検査に使うのは、先の独立監査で DAG→S61→beta まで確認した既存 sample[0]だけである。
19角 QQ minor は非零で、実 Psi が非零の核に入るので、階数は正確に19である。
free coordinate を1にした全20核係数も一致する。

元の二次共変量の全係数を、同じ10個の A_i^T 基底へ独立な左逆で移し、
全20×136 pair coefficient と全20×16微分を再構成した。16角 QQ minor は非零である。
四関係の全4×20 gradientは、非可換な dual-number matrix multiplicationで作った。
4角 QQ minorの非零と gradient*dphi=0 の全64成分を確認した。
既存の一点を別の invariant probe に置き換えていない。

## 3. 一般二重写像と関数体

Wは16次元 affine space、Zは phi=(C,D) の reduced image closureである。
rank16微分により dim Z=16、したがってphiはdominantかつgenerically finiteである。
Lの19角 minor が非零で、CまたはDの指定した成分が非零となる、Zの非空openを取る。
この open の一般点では、実 Psi の核線が L の19次 minors または19角の逆から有理的に復元できる。

K=K(Z)上でこの kernel line の一つの座標を1に正規化し、Psi=lambda*Psi0 と書く。
共変量は二次なので、非零の一成分からlambda²はKの非零元として回復できる。
Psiからbetaへの線形埋込みは単射であるから、K(W)=K(lambda)である。
一方 beta->−beta はWの非自明な自己同型でphiを固定し、lambdaを−lambdaへ移す。
従って拡大次数は正確に2であり、K(Z)=K(W)^{±1}となる。

同じ対を持つ任意の一般beta'もkernel lineに入るのでbeta'=c betaとなり、二次性からc=±1である。
char0でこのfunction-field拡大はseparableであり、さらに非零Jacobian minorのopenを取れば、
一般fiberの二点はreducedである。単にset-theoreticな二点数を一般有限次数と同一視していない。

Sp4の中心−Iは三つのcovariant indexによりW上で−1として作用する。
phiは受理済みの同時共役作用にequivariantなので、

`K(Z)^{Sp4}=K(W)^{Sp4}`。

また一般の二つの対が同時共役なら、共役を施して対を一致させ、残る±も中心作用で消せる。
よって一般orbitの分類を同時共役へ帰着できる。
regular invariant環の同型、全域有限性、退化orbitの分類はここから導かない。

## 4. 四式の一つの局所成分

Y=V(f1,f2,f3,f4)を20次元 pair affine space内に取る。
全係数恒等式によりZ⊂Yである。同じ点のgradient rank4からYの局所環は正則、次元16である。
その点で16次元既約閉集合Zが含まれるため、Yの正則局所整域におけるZのidealは高さ0で零となる。
従ってZとYの局所環は一致する。
sourceとtargetのこの正則16次元局所環の間でdphiは同型であり、完成写像も形式的同型である。

結論は、Zが四式の零点の **一つのgenerically smooth codimension4 component** であること。
Y全体の他の成分、embedded構造、Zのglobal ideal・正規性・閉像の全域記述は主張しない。

## 5. 22・24次はcharacterだけ

元の16個のCartan weightを一つずつ処理する、別の整数対称冪DPを24次まで実行した。
各次数の全weight multiplicityの和を binomial(15+d,d) と比較し、C2 Weyl群の全8項を照合した。
結果は Inv22=21、Inv24=67で、保存された全character recordと一致する。
この段階では22/24次の不変式を一点にも評価せず、基底・生成元数・NEW25上のrankも計算していない。

入口 `code/q_firstjet_two_adjoint_relations_review_20261103_v63.py`、
出力 `reports/Q_FIRSTJET_TWO_ADJOINT_RELATIONS_20261103_v63_REVIEW/`。
60秒・2GB cap、author codeと凍結primaryは変更しない。
