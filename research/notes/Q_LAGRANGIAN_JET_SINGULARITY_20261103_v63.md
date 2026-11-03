# 平行Lagrangian fibrationの一次ジェット曲線は常に特異：追加六成分の排除

2026年11月3日。標数零での構造的な全パラメータ証明である。
標本でのc4,c6,Deltaの消滅を、全成分の恒等式へ持ち上げていない。
対象は既知追加成分zero-log、M、A0、D±、C+拡張、C−磁気の六つであり、
MPS T02,...,T38には同じ結論を仮定しない。

## 1. 一般のjetに関する補題

四次元affine空間をbase二座標(x,z)、fiber二座標(y,w)に分ける。
Poisson bivector Pについて

    P^{xz}=0,
    P^{xy},P^{xw},P^{zy},P^{zw} はx,zだけの函数

と仮定する。fiber-fiber成分P^{yw}には制限を置かない。
この条件はbaseへのaffine projectionが平行なLagrangian fiberを持ち、
Hamiltonian base方向がfiber上で一定であることを表している。

任意のrank4点を原点へ平行移動する。base-first順(x,z,y,w)の値は

    C=[[0,A],[-A^T,B]],  det C=(det A)^2

で、Aは可逆、Bはskewである。
値を標準symplectic行列へ直す有限frameを

    M=[[I,0],[-(1/2)A^{-T}BA^{-1}, A^{-T}]]

と取れる。証明書では全5変数のA,Bに対して `M C M^T=J` と
detM=1/detAを直接検査する。
Mはvertical planeを保ち、base座標を変えないので、仮定の形も保つ。

正規化後のlinear jet Lについて[J,L]=0を課す。
その全12係数から二つのcompatibility関係を除いた十変数jetを保存し、
primitive perspective補正を全項適用する。
traceの微分alphaはbase成分しか持たず、v=J*alphaはfiber成分しか持たないので、
primitive補正も同じ形を保つ。
得られる全primitive jetは、八つの自由係数a,b,c,d,e,f,g,hを用いて

    Lxy=a*x+b*z, Lxz=0, Lxw=c*x+d*z,
    Lyz=e*x+f*z,
    Lyw=g*x+h*z+(b-c)*y+(a-f)*w,
    Lzw=−a*x−b*z

と書ける。これは全十変数jetからの厳密な線形恒等式であり、
一部のjetだけのansatzを全体だと仮定していない。

## 2. canonical Pfaffian curveの普遍特異点

primitive二形式beta=−J*L*Jを取り、既存のcanonical genus-one modelを作る。
Pluecker quadric Q3の五座標を(p0,...,p4)とし、

    Pl=[[0,p0,p1,p2],[-p0,0,p3,p4],
        [-p1,-p3,0,-p0],[-p2,-p4,p0,0]]

とする。Pi=span(e_y,e_w)はp=(0,0,0,0,1)である。
beta(Pi,Pi)=beta_yw=0なので、covector a_betaはこの点で零。
curveの五つの式は正確に `(Pl*a_beta, Pf(Pl))` となる。

八つのjet係数を全て自由変数のまま、この五式と全5×5 Jacobianを再構成する。
Piでは五式は全て零、Jacobianの非零となり得る列はp1列だけであり、
最後の行のその成分は恒等的に−1である。
従って全係数環QQ[a,...,h]上でJacobian rankは正確に1である。
Q3はPiで滑らかだが、Q3上の残る四式の微分は全て零となる。
曲線のZariski tangent dimensionは少なくとも3であり、滑らかなgenus-one curveではない。

幾何的にも、LG接方向ではPlのbase-base成分の一次変化は零であり、
da_betaはPiをannihilateする。従って `Pl(Pi)*da_beta=0` となる。
保存した全係数Jacobianはこの議論を直接検査している。

標数零のdegree5 Pfaffian genus-one modelのdiscriminant criterion
（Fisher–Sadek, arXiv:1406.6953, §2）により、Delta≠0ならモデルは滑らかな
genus-one curveである。上の普遍特異点はこれを妨げるため、

    Delta(beta) = 0

が全jet係数に関する恒等式になる。
もし特殊化で曲線の次元自体が増えても、Delta≠0の結論に反するため同じ零性が成り立つ。
巨大なDelta多項式を展開したという主張はせず、普遍特異点とdiscriminant theoremによる証明である。

## 3. 六つの全symbolic sourceへの適用

次の既知成分の全族を、凍結済み `additional_rank4_families_20261103_v62.py` から読む。

* R4_ZERO_LOG_5H
* R4_MINIMAL_M
* R4_A0
* R4_DPM
* R4_CPLUS_EXT
* R4_CMINUS_MAG

各族をh=1の**実際のaffine projective bivector**へ移す。
P^{0i}を持つC+拡張も `P_aff^{ij}=P^{ij}−x_iP^{0j}+x_jP^{0i}` を使い、
単なる4×4部分行列で代用しない。
全パラメータを残したまま、P^{xz}=0と、四つのbase-fiber成分の
y,wに関する全8偏微分が恒等的に零であることを直接確認する。
DPMの係数体はQQ(r), r²=3だが、形の恒等式自体はrを自由変数にしても成立する。
Poissonであることと成分の同定は、既存の全族Jacobi・dense PGL飽和・座標catalogueの定理を保持して使う。

従って各全族の任意のrank4 affine点でcanonical first-jet curveは特異である。
この条件はframeと非零scalarに依存しない。
六族のPGL飽和が各actual既知成分に稠密である既存定理を適用すると、
Delta(beta)の零性は各S61成分の**全閉包**に延長する。
従って特別なparameter値や境界点を取りこぼしていない。
sourceの正規化chartの一つでの標本検査を全像の証明に使っていない。

## 4. a=8/9の正軸への適用

既存の全Poisson曲線のeta=1におけるprimitive betaを全61座標から再構成すると

    beta_xz=(27/8)z,
    beta_xw=(3/2)x,
    beta_yw=−(2048/27)y,
    他0

である。既存Fisher–Sadek実装を正確な有理数で評価して

    c4=0,
    c6=−949978046398464,
    Delta=−522255954073519803473068032 !=0,
    j=0

を得る。このモデルは滑らかなので、正軸のeta=1は六つの成分のいずれにも入らない。

非零etaは一つのprojective Cstar orbitである。
実際、旧物理座標のdiag(1,s,s²,s³,s⁴)のpushforwardをs^{-5}倍すると、
零section Qaは固定され、etaはs^{-6}etaへ移る。
保存する全150非零係数のweightは零sectionで5、eta項で−1であり、
この有限作用を直接確認できる。C上では任意の非零eta同士が結ばれる。
従って正軸の**全非零点**が六成分すべての外にある。

既存の別packetは同じ軸をactual T01/CA/CBから排除する。
本稿はそのpacketを変更せず、独立の普遍first-jet obstructionを加える。
MPS T02,...,T38の多くでは正確な標本ですでにDelta≠0となるため、
この六族用の補題をそれらへ拡張していない。
正軸の他の既知成分への所属、新成分、全成分数の増加は結論しない。

## 保存

* `code/q_lagrangian_jet_singularity_20261103_v63.py`
* `reports/Q_LAGRANGIAN_JET_SINGULARITY_20261103_v63/report.json`
* 同ディレクトリの `certificate.json.gz`

証明書は普遍十変数jet、八変数primitive jet、全Pfaffian model、五式、全Jacobian、
一般値の有限normalizer、六つの全symbolic affine tensors、正軸の正確なinvariantsを保存する。
全入力SHA256を固定し、元のpacket・ledger・commitは変更しない。
