# 三重直線・m二 SF00 全 first symbols の独立確認

2026-11-08、PRIMARY G、別作者 `/root/quartic_review`。
root の凍結純証明 `PLANE_TRIPLE_LINE_M2_SQUAREFREE_CONSTANT_ROOT_SYMBOLS.md`
(`cfd6f8ac453f377c65ca92006df5879cfbdf4d76aa88d47b78245ce4e585c019`) と全六源 JSON を独立読解した。
**四半単純固定源、一 translation 固定源、全 (t,d)∈C×C* nil source への必要被覆を受理する。**
G 上界、higher compatibility、E membership は本 receipt の結論ではない。
CAS、source/operator/minor 計算、jet 又は GB は実行していない。

## 1. 全 irregular pool と strict Jacobi

対象を f=a³、normal 位数二 SF、両 primitive roots が constant の全源に限定する。
accepted discriminant identity ΔG2=constant*f² と実 constant normal GL2 により G2=γfuv、γ!=0。
二根を保存する全 ordinary conormal connection は diagonal であり、log 性は仮定しない。
tracefree mixed field との式 θ=(dlogf−αX)/3 と flatness dθ=0 は X(f)=0 を与える。
従って Xa=0、同様に Ya=0。

これで全 pool を直接再導出できる。
X=(0,H1(b,c)+a*l1)、Y=(0,H2(b,c)+a*l2)、tr Hi=0 で各五次元。
この五次元は d(P2/a²)、P2 mod Ca² の全 irregular exact terms と一致する。
実際 a*b,a*c,b²,b*c,c² の五 potentials がその全線形場を張る。
非対数成分を log characters へ置換していない。
target に u²/v² が無いので normal K は κuv のみ。

全 a-row は零である。
E∧Jac=0 の quad(a,i,j,k) は a*Jac_ijk=0 なので、P0 は strict Poisson でもある。
divergence-zero affine fields は anchor 3a²∂b∧∂c を保存し、残 Jacobi は exactly [X,Y]=0。
従って

    [H1,H2]=0,    H1*l2−H2*l1=0.

二 affine Hamiltonians の bracket は constant、又は area-form の Cartan identity から、
det(Xbc,Ybc)=d0*a²。actualG2 の全係数は γ=3κ−d0!=0。

## 2. κ!=0 の全補集合と四半単純源

全 raw36 P1、Euler gauge 六、P2 九を戻した二つの constant-root curl は記載通りである。
κ!=0 なら両 curl の a-row は C_a=C'_a=0 を強制する。
従って translations l1/l2 は残 curl kernel を変えない。
commuting sl2 pair が nilpotent/zero だけなら scalar±2κ と共通 nilpotent の和は可逆で矛盾。
一方が半単純 nonzero なら同時対角化でき、H1=diag(λ,−λ)、H2=diag(μ,−μ)。
二必要固有値式の全 signs から

    4μ−λ=±2κ,    4λ−μ=±2κ

を解くと λμ!=0、μ=λ, κ/λ=±3/2 又は μ=−λ, κ/λ=±5/2 のみ。
H1 可逆と affine commuting 式から、同一 base translation が両 l を消し、d0=0。
old base=g*new base、tensor scale σ=det g で anchor を保ち、det g=1/λ と取れば λ=1。
これで記載の四 fixed sources が全てである。
normal roots の交換や base 交換による重複があっても必要被覆は変わらない。

## 3. κ=0 の全二枝と actual inverse

γ!=0 より d0!=0。
半単純な非零 H があれば共通 affine fixed point を消した後 X/Y は proportional となり d0=0。
よって全 H は nilpotent/zero。

両 H=0 なら translation matrix L=(l1,l2) は invertible、δ=det L!=0。
σ=1/δ、old base=(σL)*new base とすれば
σ/det(σL)=1、σ(σL)^−1L=I。
従って anchor を保ち、X=a∂b,Y=a∂c の一固定源へ戻る。全 inverse は nonzero matrix の逆である。

少なくとも一 H!=0 なら根を交換して H1!=0、sl2 centralizer より H2=tH1、全 t∈C。
D=l2−t*l1∈ker H1=im H1、d0=det(l1,D)!=0 なので D!=0、l1 はその line 外。
base affine translation で l1 の image 成分を消し、basis を取って

    Xb=r*a, Xc=s*b, Y=tX+d0prime*a∂c,    r*s*d0prime!=0

とできる。
old b=βb、old c=χc、tensor scale σ=βχ とすると anchor は固定され、
Xb の係数は χr、Xc の係数は β²s。
χ=1/r、β²=1/s の両 nonzero 選択を保てば

    X=a∂b+b∂c,
    Y=tX+d*a∂c,    t∈C, d∈C*, κ=0, γ=−d.

全 shear/diagonal inverses は存在する。
d を一点へ固定した slice ではなく、新しい任意 nonzero parameter のまま保持する。
H1=0、D=0、translation rank低下、γ=0 はそれぞれ上の別枝又は対象外で、反転による未記録の除外は無い。

## 4. 正確な結論

全 first-symbol JSON の五固定源と全二parameter nil 源は上の式と一致する。
これは SF00 の whole necessary reduction であり、他 root flags、f=a²b、normal double、m三以上は別残余。
全 P1三十/P2九は未制限のままである。
後続仕事では nil の全 (t,d) domain、全 rank/denominator walls と inhomogeneous P2 を保持する必要がある。
今回の pureproof 受理だけから six higher sources の G bound 又は E membership を結論しない。
