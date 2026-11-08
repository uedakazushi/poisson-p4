# Tangent conic＋line・m二doubleの全点被覆：独立全源算術と正の親の監査

2026-11-08、PRIMARY E、別作者 `exceptional_coverage`。**全source／全点被覆 PASS**。
`PLANE_M2_DOUBLE_TANGENT_POINTWISE_COVERAGE.md` の全依存を読解し、
新しいQQ(t)＋全七壁のsource逆を一つの別Fraction/RF batchで再構成した。
新M1、operator、point、rank、kernel、profile、GB又はjetは実行しない。
新component数は零。

## 1. 全八sourceの独立算術

job `jobs/plane_m2_tangent_whole_source_separate_review.json` はactive gate PASS。
標準ライブラリFraction/RFと旧独立checkerの基本疎多項式だけをimportし、producer/CASはimportしない。
1.82164秒、最大RSS27,548 KiB、3171個のexact identityが全PASS。
実行sourceはjobと結果JSONに埋め込んだ。

元full150で保存された30個のweight一basis、全九weight二basisと全八caseの完全K1/J2を入力とする。
新しいM1又はkernelは構成せず、以前の独立非零minor・完全kernel証明を再使用した。
別にfull Eulerの六gaugeを作り、raw36の定数inverseを確認した。
全六actual unit shear `a_i'=a_i+s_i*u+t_i*v` をfull tensorへpushして、
weight一をこの正しいEuler商へ戻した全30×6 matrixが保存Sと一致した。
weight一partはgradingによりshearに線形なので、これは全六parameterのshear mapである。

各caseでK1*Cが選択shear columns＋全残kernel columnsに等しく、CとC逆の両積がidentity。
したがって完全K1から全残sourceへの実逆があり、単なるdimension比較ではない。
t=-6だけshear五、残K1二つである点もそのまま保持した。

残P1と全九P2を戻し、元J2の証明済みpivot rowsの**実full affine方程式**を再構成した。
保存P2 particularと全free variablesを代入して全selected equationが零となることを確認し、
全normalized150、全350残Jacobi係数、全126actual Euler-Pfaffianを別に計算して比較した。
元J2の非零minorを再使用するため、この解はselected方程式の全解でもある。
homogeneous kernelだけの確認ではない。

| source | 実shear | 残K1 | free P2 | 結果 |
|---|---:|---:|---:|---|
| QQ(t) | 6 | 0 | 0 | c因子 |
| 6 | 6 | 0 | 0 | c因子 |
| -3/2 | 6 | 1 | 0 | c因子 |
| -5/2 | 6 | 0 | 0 | c因子 |
| -3 | 6 | 2 | 0 | 全二parameter Poisson、下のT40認識 |
| -9/2 | 6 | 0 | 2 | 全二parameter Poisson、下のT37認識 |
| -6 | 5 | 2 | 1 | 互換r0*r1=0、全許容値でc因子 |
| -9 | 6 | 1 | 0 | c因子 |

t=-6のP2_abには `-r0^2*v^2/3` が必要であり、これも保持して再構成した。
その全非零残Jacobiはr0*r1の倍で、同条件の下でactualDの全項がcで割れる。
QQ(t)のK1/C/C逆と選択forcing minorの全有理係数について、
分子・分母の非定数因子はt零と元七壁の線形因子だけで尽くされることをFraction多項式除算で確認した。
全七壁は元の完全kernelから別に戻してあり、未処理complex rootを残さない。
t零はnormal order三以上の別入口であり、今回のm二domainへ混入しない。

## 2. -9/2全族のactual orbitとT37への正の適用

全P(p,q)に対し、
`a'=a+2lambda*b+lambda^2*c,b'=b+lambda*c,c'=c,u'=u,v'=v-lambda*u`
はdet一、逆はlambdaの符号反転。
P0は全tで不変、残P2ではPacの追加が10lambda*q*u²となり、
Pabのu²項は-9+8+1で相殺される。
全entryで正確に `P(p,q)->P(p+10lambda*q,q)`。
q非零でpを消し、共通normal scale eta²=qでq一へ戻す。
q零は全A2 polynomial Poisson族の稠密openの閉包として同じ閉親に含まれる。

canonical点の第二planeはS'=(c=u=0)。そのidealを全raw括弧が保ち、restrictionは
`Pab=-3b²-9v²,Pav=-6bv,Pbv=0`。
zero cubicは `3v(b²-3v²)`、三つの相異なるconcurrent lines。
actualDのS'に沿うorderは三、leading cubicのdiscriminantは
`-2916*v^4*(b²-3v²)^2` で非零。
従って既受理の**全point concurrent m三定理**を適用できる。
leading divisor形又はgeneric cleanだけを使った適用ではない。

この閉和はT37、K22の既知23軌道閉包、essential-three八成分。
canonicalDはessential五なので最後の閉determinantal条件を満たさない。
被約特異支持は正確に二planeで、K22の全unit源が作る単一actual軌道の三planeとは異なる。
既保存ambient actual M48が与えるorbit23の下界を再使用すれば、K22軌道閉包の次元22以下の真境界にも入らない。
ここで既存ambient rankを新計算したとは呼ばない。
したがってK4aの全源へ戻り、既actual PGL同定によりT37に入る。
元全unit・zero境界・実shear逆を戻して、-9/2の全源がT37。

## 3. -3全族のactual Hamiltonian逆

producerの全P(r,s)について、canonical divergence-zero liftのratio rho=v/uに対するHamiltonianは

```
Ma=-9a/5+3c*rho²,
Mb=-9b/5-3c*rho-u*(r*rho²+s*rho³),
Mc=6c/5, Mu=6u/5, Mv=6v/5.
```

row divergence Ru=0,Rv=6uを用いて再導出した。
表示されたglobal Borel＋六shiftのcoefficientsを、MaとMbのrho各次数で別に確認した。
特にMaではrho³からf=-2lambda*s/3、rho²からe=-(2lambda*r+4lambda²*s)/3、
rho/定数からB=-2lambda*e+lambda²*f、A=lambda²*e。
Mbのrho/定数はproducerのD,Cに一致する。
よって全actual mapの下でrがr+5lambda*sへ移るHamiltonian恒等式は全parameterで正しい。

これだけからtensor同値を推測しない。
旧全K1は六shear＋r,sの八方向で尽くされ、Maのrho³/rho²/定数/rhoでf,e,A,B、
続くMbのrho³/rho²/定数/rhoでs,r,C,Dを一意に読める。
したがってHamiltonian mapはこの全K1上injective。
残P2は元J2 rank九により存在する点で一意であり、全(r,s)でP2=0の実解が上のbatchで確認されている。
これでactual projective tensorの等式 `g_lambda*P(r,s)=P(r+5lambda*s,s)` が従う。
u非零は等式を証明するchartだけで、g_lambdaはglobal polynomial linear map・det一の三角逆を持つ。
正則tensorの等式はu零にも延びる。

s非零でlambda=-r/(5s)、続いてnormal common scale sでP(0,1)へ戻す。
s零も全A2の同じ軌道閉包に入り、全zero sourceを保持する。

## 4. -3 canonical点の全X122と、実T40 theorem

canonicalDはC(c,u,v)上で、A=a-c(v/u)²、B=b+c(v/u)+u(v/u)³/3により
`3u²*(u*A²-v*A*B-3c*B²)`。
discriminant `v²+12cu` はrank三quadricでsquareでなく、元polynomialのcontentは一。
従ってC絶対既約でactual rank四。

全normalization mapは

`[a:b:c:u:v]=[s*x:t*x-s²*z:s²*y-st*z:3t²*y:3t²*z]`

であり、sourceはX122。
t非零でu/vの二columns、t零でa/b/cの三columnsが独立なので全baseでrank三。
tautological Hはample、H³五、imageのdegreeも五なのでfinite birational、全normalizationである。
u非零又はB非零で導いた有理式だけでboundaryを扱ったのではない。
plane導手CPはt²で、node stに置き換えない。

u一chartの**projective**bracketを用いると
`{A,B}=0,{A,c}=Av+3cB,{A,v}=-3A,{B,c}=A,{B,v}=-3B,{c,v}=0`。
D上でeta=A/BはCasimirで、pi_D=-3B∂B wedge ∂v、従ってX122のbaseにvertical。
dualizing volumeを引き戻してcontractすると正確にdeta/3、CP=t²のRosenlicht generatorである。
rankzero S5もA=B零の射影閉包で、cubic systemの唯一のdouble basepointを解くとH²=五。
従って既受理genus零rank零・degree五theoremの仮定を満たす。

有効な正のtheoremはauthor cc6b209c…と新独立review02acb785…である。
全五sectionのQuot変形、MCM partial-conductor base change、中央一線との一致とambient Poisson liftにより
この任意cusp sourceがT40に入る。
旧失われた80e28… identity、未確認cyclic normal form又は単なるprojection閉包を使わない。
全unit軌道とzero境界・元shear逆を戻して、-3の全源がT40。

## 5. Constant源と最終全点範囲

全leading入口のprimitive root度零／一と、そのzero境界は既独立normal form theoremを再使用する。
constant h=fのZ0/ZNは、全raw36＋full9P2を保ったratio-Hamiltonian／nilpotent inverseにより
actualDがfraction field上でも二つの正次数base因子に分かれ、既約quinticでは空。
これは `PLANE_M2_CONCURRENT_DOUBLE_SEPARATE_REVIEW.md` §2でtangent Z0を含め受理済みで、再runしない。
constant h=L³は正しいOre V2と別reviewだけを用い、essential-three八成分の全点定理へ戻す。
旧V1の誤ったstrict-Jacobi零を使わない。

以上によりreduced tangent conic＋line、normal order二、generic double normal rootの全sourceは、
既essential-three八成分／T37／T40でpointwiseに被覆される。
全generic逆・全七壁、inhomogeneous P2、rank-jump、source zero、cusp、frameを保持する。
normal order三以上、t零の別leading profile、他R又は全Eの完了はこのreceiptから主張しない。
producer、global、main、Gitは変更しない。
