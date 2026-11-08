# Twisted-cubic cone・rank二：全源はT38／T39又は既約Dで空

2026-11-08、PRIMARY E、producer完成報告。独立reviewは別receiptへ固定する。
新componentの増分は零。全global台帳・主稿・Gitは本packetで編集しない。

**定理。** P4上のrank四Poisson構造の実退化因子Dが被約既約quinticで、
非正規支持にtwisted-cubic cone Sがあり、そのrestriction rankが一般に二なら、
この全sourceは既知閉成分T38又はT39の合併に入る。
Rの既約／被約、clean、ordinary normal node又はvertexのsmoothnessを追加しない。

## 1. 完全な元源と全補集合

既独立受理 `CUBIC_CONE_RANK_TWO_SEPARATE_REVIEW.md` により、
元975×150、Euler24商、完全canonical27係数、全350Jacobi／全126actualPfを保持する。
base九＋vertex十八へ分けた被約支持はactual SL2／nonzero scalarと明示逆により
A_ss(λ),A_nil(λ),B_nil(λ),W_ss,W_nilの五whole familiesに尽くされる。
全λ、全十九minor walls、無限遠W-only、base零を戻している。

`CUBIC_CONE_RANK_TWO_VERTEX_SHEAR_REDUCTION.md` と
`CUBIC_CONE_RANK_TWO_CRITICAL_ORBITS.md` の全pure還元を適用すると次の表になる。

| 全source | 指定するrank四・既約D枝の結論 |
|---|---|
| A_ss/A_nil、λ≠2/5,3/5 | 全四vertex shear除去とconstant-rank kernel bundleからactualD properfactor又は零。λ=1/2と全minor wallsも保持 |
| A_ss/A_nil、λ=2/5 | 全六vertex sourceは旧rank零cone sourceとactual canonical liftで等しく、全T38 |
| A_ss/A_nil、λ=3/5 | 下の陽性Fisher認識と全実逆／閉包から全T39 |
| B_nil、全λ | canonical CasimirからN(F)=Z(F)=0。次数五の全六highest vectorsはいずれもproper productなので全空 |
| W_ss/W_nil、base零 | 全ΠがEulerを除いてcommon ∂c wedgeを持ち、rank二以下 |

固定λ=2/5,3/5のsemisimple／nilpotent両base orbitで全kernel rankは六。
nonzero binary quadratic U全体上のkernel vector bundleがirreducibleで、semisimple逆像が稠密。
従って全nil fiberは対応するsemisimple全源の閉像へ戻る。
generic denominator代入だけでnil全源を覆ったものではない。

## 2. λ=3/5のwhole sourceと唯一の代表

保存全六kernelのactualFは

```
F=(c+ell)*q0*q2+a*q2*Cminus+b*q0*Cplus,
ell=(5k4/2)x0−5k1*x1−(5k2/4)x2−(5k3/2)x3,
a=k0, b=k5/2.
```

全D写像のinjectivityにより四actual c-shearでellを消し、tensor自体が一致する。
a0又はb0はproperfactor、ab≠0はprimitive c-linear Gaussにより絶対既約。
actual Cartanとvertex scaleは(a,b)を(gamma*t^5*a,gamma*t^-5*b)へ動かし、
全ab≠0をa=b=1へ戻す。十乗根は有限選択、反転したunitは明示したabだけ。
全actual逆は凍結済みpure noteにあり、sourceを標本一点に置換していない。

したがって残る義務はこの一つの代表の陽性帰属だけである。
canonical integer代表Pは分母消去scale5、k0=1,k5=2、他零。
全sourceに対するactual軌道次元24は既証明だが、そのcomponenthoodは帰属に不要。

## 3. explicit nodal Fisherモデルの全係数陽性一致

`CUBIC_CONE_RANK_TWO_NODAL_FISHER_MODEL.md` のPhi(e)を使う。
五signed PfaffianはWeierstrass cubic y²=x³+x²+eのdegree五モデルの
五quadricsの正確に−一倍。generic e(4+27e)≠0でsmooth elliptic normal quinticである。
既受理Fisher／Polishchukのsmoothモデル定理によりOmega(Phi(e))はT39構造。
Omegaは全model係数の多項式であり、e又はdiscriminantを除していない。

constant actual map z=Mxは前の二noteのものを用い、detM=16384、正逆両積はidentity。
有限receiptは正確に

`M^-1 * Omega(Phi(0))(M*x) * M^-T = 2*P(x)`

を**全150係数**で確認した。Euler補正すら不要なliteral equalityである。
右辺は非零。したがってpolynomial族のe=0でのclassもclosed T39に入る。
singular nodal曲線へsmooth-only theoremを無断適用したのではなく、
実際のsmooth generic fiberからのregular polynomial specializationである。

同finite jobは全eを保つexact global族についてさらに以下を確認した。

* 五Pfaffian、Phi*Pf=0、全modelと全normalization pullback。
* 全十rawJacobi、全350projectiveJacobi、canonical divergenceは恒等零。
* 全126Euler-Pfaffianは五quadricsのJacobian determinantの−25/32768倍。
* e=0の五quadric basis changeはdet1073741824、元q0,q1,q2,r1,r2のJacobianは−2F。

得られたtarget tensor族は実際にeに**一次**である。
これはfinite jetでなく、全eでJacobiを満たすglobal actual curve。
active gate後のONE実行は2.66223秒、RSS59,076KiB。
報告は `CLASSIFICATION_CUBIC_CONE_RANK_TWO_NODAL_FISHER_ACTUAL_20261108.json`。
Casimir full150 operatorは実行していない。

## 4. 全zero枝／nil fiberへ帰属を戻す

ab≠0全sourceは前節代表のactual orbitに入り、T39はPGL-stableかつclosed。
semisimple全六fiberはirreducibleで、ab≠0が稠密なので、ab=0も同じT39に入る。
§1のconstant-rank bundle橋によりλ=3/5のnil全六fiberもT39へ戻る。
この全zero枝包含は既約D外でも成立するが、他のsourceの別factor型分類を完了した意味ではない。
λ=2/5の全T38橋と他の全空枝を合わせ、冒頭の全pointwise E被覆を得る。

## 5. 先行の非閉鎖結果を保ったままの結論

D100/M49とquadratic余核零は不変保存する。高次obstructionを実行していない。
T39をdim28と置く仮のglobal矛盾は誤りであり、実dim25のため使わない。
条件付きglobal≤24監査も、この陽性包含の前提ではない。
今回の成功はexplicit genus-oneモデル／全actual tensor familyのT39認識である。
新component、smoothness、nil指数又は全Poisson分類の完成は主張しない。
