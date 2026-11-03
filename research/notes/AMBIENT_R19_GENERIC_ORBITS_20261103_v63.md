# R19 の一般軌道：正確な一変数標準形と完全な不変量

2026年11月3日。対象は全射影被約交叉の19次元階数二成分Rである。
全Poisson空間の成分数を増やす結果ではない。
以下は元の八係数族で `q4*q5*q7≠0` を満たす点と、そのPGL5飽和を分類する。
これはRの稠密な部分を含むが、q4=0またはq5*q7=0の全特殊軌道は分類しない。

## 定理

上の稠密な部分の任意の点は、あるκ∈Cを用いて

    Pκ=∂y∧Vκ,
    Vκ=4x²∂x−2xz∂z+(2xw+h²+κhz+z²)∂w

と射影同値である。二つの標準形Pκ,Pκ′がPGL5同値であるための必要十分条件は

    κ²=κ′²

である。元の八係数での完全な不変量は

    I=κ²=4q6²/(3q7q5)

となる。
この部分の全軌道次元は18である。一般点のprojective stabilizer次元は6で、
Rの一般軌道には厳密に一つのmodulusがある。
特に軌道次元の上界を二つの標本微分階数から推定したものではない。

## 1. 八係数から三係数への全有限変換

t=q4と置く。実際のR族はP=∂y∧Vで

    V(h)=0,
    V(x)=4t*x²,
    V(z)=−q0*x²−2t*xz,
    V(w)=2t*xw+3q7*h²+2q3*hx+2q6*hz+q1*x²+q2*xz+q5*z².

まず

    a=q0/(6t),   z1=z+a*x

とすると、V(z1)=−2t*x*z1となる。
残るforcingのx項を

    B1=2q3−2a*q6,
    D1=q1−a*q2+a²*q5,
    E1=q2−2a*q5

と書き、次に

    w1=w+(B1/(2t))*h−(D1/(2t))*x+(E1/(4t))*z1

とする。この二つはdeterminant1の可逆線形変換で、全パラメータのまま正確に

    V(w1)=2t*x*w1+3q7*h²+2q6*h*z1+q5*z1²

を与える。
projective scaleでtを1に戻せば

    VQ=4x²∂x−2xz∂z+(2xw+Q(h,z))∂w,
    Q=A*h²+B*hz+F*z²,
    (A,B,F)=(3q7/t,2q6/t,q5/t)

となる。
三次以上を無視した正規形ではなく、元の実際の二次双ベクトルの全有限恒等式である。

## 2. intrinsicな定ベクトル因子と全invariant hyperplanes

U=C5/C∂yを(h,x,z,w)で記述する。
div_U VQ=8xなので、canonical双ベクトルは∂y∧C、ただし

    D=5C=(-8hx,12x²,-18xz,5Q+2xw)

である。

### 定ベクトル因子は一意

intrinsic condition `Eul∧v∧Pcan=0` を満たす定ベクトルvを考える。
vを∂yで割った成分をu∈Uとすると、条件は`E_U∧u∧D=0`となる。
x=0ではD=5Q∂wで、Q≠0だから、全(h,z,w)に対する外積条件から

    u_x=u_h=u_z=0

が従う。残るu_wについて、h,x,w成分は

    u_w*(h*D(x)−x*D(h))=20u_w*h*x²

となるためu_w=0でもある。
従ってvは∂yの倍数に限られ、L=C∂yは全Q≠0で一意なintrinsic factorである。

### quotient上の全invariant linear hyperplanes

ここで「invariant」はD(ell)がellで割れることを意味する。
ellのw係数が非零なら、ell=w+a h+b x+c zと正規化できる。
w=−a h−b x−c zを代入した割り算の余りは正確に

    5Q−10a*hx+10b*x²−20c*xz

となる。Qはh²,hz,z²だけを含み非零なので、この余りは零にならない。
よってw係数は零である。

ell=a h+b x+c zなら

    D(ell)=x*(-8a h+12b x−18c z).

ellがxの倍数でなければ、一次式の整除性からellはdiag(-8,12,-18)の固有covectorでなければならない。
三固有値は相異なるため、全invariant hyperplanesは正確に

    h=0, x=0, z=0

の三つだけである。一般点でのDarboux hyperplaneの尽くしを仮定せず、上の全係数計算で証明した。

そのうちxだけが`D(ell)∈C ell²`を満たす。
さらにh,zのcofactorをxのcofactorで割った比はそれぞれ

    -8/12=-2/3,    -18/12=-3/2

である。したがってhとzは入れ替えられず、三つの線はintrinsicにラベル付けされる。

## 3. 任意の射影同値を対角変換へ制限する

二つの標準形がPGL5同値なら、前節の一意なLを保つため、商Uで線形変換が誘導される。
canonical代表は一意なので、そのquotient field DはGL4線形共役と非零倍率だけで移る。
三つのラベル付きinvariant covector linesは保たれ、従って商の行列は

    h'=a*h, x'=b*x, z'=c*z,
    w'=d*w+u*h+v*x+r*z

という形でなければならない。
x成分を4x²へ戻す倍率を含めた実際の変換後のw成分は

    2xw+(b*d*A/a²)h²+(b*d*B/(a*c))hz+(b*d*F/c²)z²
                      −(2u/a)hx+(2v/b)x²−(4r/c)xz.

別の標準形には最後の三項がないので、u=v=r=0である。
canonical比較でEulerの曖昧さを残してはいない。
intrinsicにラベル付けされたhをCasimirとするraw代表をDから戻すと、上のVQが一意に戻り、
xの倍率もD(x)から定まるので、同じ計算をそのまま使える。

したがって標準形間の全同値は、quotientの対角作用で尽くされる。
y方向には独立にscaleと任意の線形shearを加えられるが、∂yとの外積で消えるのでQの比を変えない。
全GL5でのcanonical150係数の有限変換式も保存した。

## 4. 完全不変量と十八次元軌道

三係数の変換則は

    A'=b*d*A/a²,
    B'=b*d*B/(a*c),
    F'=b*d*F/c².

従ってB²/(AF)は不変である。
A,F≠0ならC上で平方根を選び、a=√A,c=√F,b=d=1としてA'=F'=1にできる。
残るB'=κは符号まで定まり、h→−hがκ→−κを実現する。
これでκ²の一致が必要十分であることが分かる。

同じ標準形を保つ対角変換は、B≠0なら

    b*d=a²=c²,  a=c

を満たす。自由対角パラメータは二つである。
これにyのscale一つと、yへのh,x,z,wの四つのshearsを加えるとGL5 stabilizer次元は7である。
中心の一次元を除いたPGL5 stabilizer次元は6となり、軌道次元は24−6=18である。
B=0ではa=−cの成分も許されるが、stabilizer次元と軌道次元は変わらない。
AF=0やt=0へこの記述を広げることはしない。

元のsource openはRへ稠密に写り、そのPGL飽和は稠密なconstructible部分を与え、従って非空の開部分を含む。
上のintrinsic flagと係数比によりIはその一般点でrationalに定まり、
完全な一般軌道不変量となる。
同値な言い方では、支配的なG×A1_κ族への制限により

    C(R)^PGL5=C(κ²)

である。特殊軌道の全分類を意味するものではない。

## 検証と未完了の範囲

* `code/ambient_R19_generic_orbits_20261103_v63.py`
* `reports/AMBIENT_R19_GENERIC_ORBITS_20261103_v63/report.json`

元の全八係数の可逆線形変換、divergence、全invariant-hyperplane割り算、対角・w-shift作用、
y-shearを含む全canonical有限作用、κ符号変換をQQ(parameters)で検査した。測定0.25秒。
t=0、A=0、F=0の特殊軌道、Rの一般点における交叉scheme厚み、R自体の正規性は未解決である。
