# 既知五点の普遍補間と五生成式cofactorの厳密構成

2026-11-04。これは既知五点の標準模型に対する有限恒等式計算であり、実K5係数の評価、Gram126、normal展開、残余の空性を計算しない。

## 1. 模型とprimitive element

QQ[A,B,sigma]の五生成式を

```
f0=1+sigma-2sigma²-A²+AB-B²,
f1=A-A²+2AB, f2=B+2AB-B², f3=A*sigma, f4=B*sigma
```

とする。五点は(1,0,0),(0,1,0),(-1,-1,0),(0,0,-1/2),(0,0,1)。三つの独立形式係数alpha,beta,gammaに対しh=alpha*A+beta*B+gamma*sigmaと置く。五つのh値thetaはalpha,beta,-alpha-beta,-gamma/2,gammaで、Vandermonde D=product_(i<j)(theta_i-theta_j)は非零である。

g5(H)=product_i(H-theta_i)はmonic五次式である。三座標のinterpolationをp_A,p_B,p_sigmaと書く。計算では分数関数を展開せず、D/d_i=(-1)^i*product_(k<l,k,l!=i)(theta_k-theta_l)を使ってD*p_lを多項式のまま構成する。五点での全座標値とg5(theta_i)=0をQQで検証する。

## 2. Cofactorとcoherence

標準単項式は1,A,B,sigma,sigma²。A²,AB,B²,A*sigma,B*sigma,sigma³の六個のrewriteは、最初の三式の3×3係数行列逆と最後の二式、およびsigma³の明示恒等式から得る。graded lex A>B>sigmaの固定規則を全targetへ共通に用いる。これはparameter環上線形なreduction mapである。

保存する四つの完全恒等式は

```
D*(A-p_A(h))=sum_i C_Ai*f_i,
D*(B-p_B(h))=sum_i C_Bi*f_i,
D*(sigma-p_sigma(h))=sum_i C_sigma_i*f_i,
g5(h)=sum_i C_gi*f_i.
```

全係数をQQで再乗算して確認する。前三本のcofactorのplane次数は2以下、最後は3以下である。また、補間の一意性だけでなく、次の全五本のcofactor恒等式も実際に確認する。

```
alpha*(D*p_A)+beta*(D*p_B)+gamma*(D*p_sigma)=D*H,
alpha*C_Ai+beta*C_Bi+gamma*C_sigma_i=0 (i=0,...,4).
```

従ってh定義とnormal correctionの整合性を失うような独立cofactor選択をしていない。

## 3. 実源への適用範囲

実toric三角形を(a,b)=(xA+yB,zA+wB)へ正確に同定した場合には、alpha=x+tau*z、beta=y+tau*w、gamma=tau²を代入できる。これは標準模型に関する条件付きの適用recipeである。本計算ではそのactual GL3/ell-normalization adapterを構成・評価していない。Galois permutationに対するinterpolationの一意性は降下を与えるが、cofactorを実係数で保存するにはそのadapterまたはS3平均を明示的に束縛する必要がある。

実42法線源の5×675表示、DU=ell*e0所属、残余の空性は依然未計算である。既存のglobal monic-cover recipeに対する小さい普遍恒等式部品を保存する計算である。

## 4. 再現と制限

`code/s61_universal_five_plane_interpolation_20261104_v63.py --output work/NEW-output`。一回60秒・2GB、全保存多項式30000項・係数512bit・RAW12MiBのpreflight上界を上限にする。新しいβ点・素数・GB・field Gaussianは使わない。証拠は全多項式係数、五点の値、六rewrite、四組×五cofactorを含む。

## 計算結果

承認された一回は0.2198秒で完了し、すべてのQQ恒等式に合格した。保存多項式の延べ項数5842、最大有理係数11bit、計算途中の最大support984項。四targetの項数は605,605,370,123、固定border reductionは30,30,30,51段だった。証拠はRAW844040 bytes、gzip31927 bytesである。actual adapter、K5のscalar係数、Dの構成、所属探索は実行していない。
