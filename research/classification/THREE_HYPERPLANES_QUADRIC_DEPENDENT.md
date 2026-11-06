# 従属する三超平面と既約二次因子：全支持の分岐と実変形

2026-11-07。E_reduced_nonnormal。実退化因子が
`H0+H1+H2+Q`、三つの相異なる超平面の一次式が2次元を張り、Q が既約二次式である場合を扱う。
以下の全係数比較と実変形は手計算で確認した。有限標本による一般性の推論は用いない。
第1–5節の分岐は全支持を保持する。第6節の残余全源について、第8節の有限支持計算は
その飽和閉包が24次元全成分となることを与える。
**形式切片からの論理と有限算術の独立再構成はともに受理済みである。**
従って従属する三超平面の全分岐は MPS40、CMINUS_MAG、および本24次元成分で覆われる。

## 1. 全点を保持する Poisson 方程式

一超平面を h=0 として残りを x=0,x−h=0 とする。
h-Casimir 二次代表を取ると、各因子の保存から
`{x,j}=x(x−1)c_j`、c_j は定数となる。
一般階数4なので定数ベクトルは非零であり、残る三座標を変えて

\[
 f=x(x-1),\qquad\{x,z\}=f,\quad\{x,y\}=\{x,w\}=0
\]

とできる。x を含む三つの Jacobi 式は残りの全係数の z 微分が零であることを与える。
Pfaffian から `{y,w}` は二次因子 Q の非零定数倍なので、その定数を Q に吸収する。
従って全構造は正確に

\[
 P=V\wedge\partial_z+Q\partial_y\wedge\partial_w,
 \qquad V=f\partial_x+A\partial_y+B\partial_w,             \tag{1}
\]

Q,A,B は x,y,w の次数2以下の多項式、Q は既約である。
残る一つの Jacobi は

\[
 \boxed{\ fQ_x+AQ_y+BQ_w=Q(A_y+B_w)\ }.                  \tag{2}
\]

逆に (2) を満たす任意の Q,A,B は実際の strict Poisson 構造を与える。
座標順 `(x,y,z,w)` の Pfaffian は `−fQ`、実退化五次式はその次数4斉次化に h を掛けたもの。
Q は z に依らないので、その射影二次式の階数は3または4であり、階数5は存在しない。

`z'=z+λy+μw` は A,B を `A+μQ,B−λQ` に変える。
以下で定数倍の Q を除く操作は、この実 affine shear で行う。
全二次式 Q の (y,w) に関する Hessian の階数 2,1,0 に分岐する。

## 2. Fiber Hessian の階数2

定数線形変換と x に関する affine 平行移動によって
`Q=yw+q(x)` とできる。q は次数2以下の非零多項式。
(2) の y,w 係数を全て比較すると、正確に

\[
 A=\alpha Q+a(x)y,\qquad B=\beta Q+b(x)w,
 \qquad q(a+b)=fq',                                      \tag{3}
\]

α,β は定数、a,b は affine 一次式となる。
Q 倍を shear で除く。`q|fq'` より q の有限根は0,1だけである。
斉次化した binary quadratic の根は三超平面のいずれかに限られる。
三超平面の並べ替えと尺度変更で、全支持は次の二族に入る。

\[
\begin{array}{ll}
\text{I:}&Q=yw+1,\quad V=f\partial_x+a(x)(y\partial_y-w\partial_w),\\
\text{II:}&Q=yw+f,\quad V=f\partial_x+a(x)y\partial_y+(f'-a(x))w\partial_w.
\end{array}                                               \tag{4}
\]

a(x) は任意の affine 式。I の射影 Q は階数3、II は階数4である。
両族の全支持は第7節の実変形により MPS40 成分の和に入る。

## 3. Fiber Hessian の階数1

正規化して `Q=y²+a(x)w+q(x)` とする。a は affine、q は次数2以下。

### 3.1 a=0：CMINUS_MAG の全磁気源

既約性は斉次 binary q の階数2と同値である。
全係数比較は

\[
 A=\alpha Q+u(x)y,\qquad
 B=u(x)w+B_0(x,y),\qquad 2qu=fq',                         \tag{5}
\]

B0 は任意の二次以下の式を与える。
従って q の二つの根は三超平面の中の二つである。
三超平面を並べ替えて `q=f` とし、Q 倍を除けば

\[
 Q=y^2+f,\qquad
 V=f\partial_x+{f'\over2}y\partial_y+
       \left({f'\over2}w+B_0(x,y)\right)\partial_w.       \tag{6}
\]

これは `legacy/v08/sections/cminus_component.tex` の全二次磁気項を持つ源と同一である。
同資料では binary f を `h²+X²`、Q を `h²+X²−2Y²` としている。
相異なる二根の binary 形式の座標変更、非零尺度、横変数の尺度変更によって (6)と一致する。
同源の一般点の成分認識と源の既約性から、B0 の特殊値も含めた全 (6) は
CMINUS_MAG の閉包に入る。磁気項を零にして得た結論ではない。

### 3.2 a≠0：I への実変形

Q 倍を除くと全解は

\[
 A=u y-{a c\over2},\qquad B=c y+u w+e,\qquad
 u={fa'\over a},\quad ae+fq'=2qu,                         \tag{7}
\]

c は affine、e は次数2以下である。
多項式性と次数条件から、a の斉次化は三超平面のいずれかに比例し、
q の斉次化は a で割り切れる。
実際、a が非定数なら `a|fa'` から根は0または1で、(7)のその根への代入は `a|q` を与える。
a が定数なら u=0 で、e の次数上界から q は affine、すなわち斉次化は h で割り切れる。
affine w 平行移動で q=0 としてよい。

全源は次の二次式の実変形を持つ。

\[
 \begin{gathered}
 Q_\varepsilon=y^2+a w+\varepsilon w^2,\\
 A_\varepsilon=u y-\tfrac12(a+2\varepsilon w)c,
 \qquad B_\varepsilon=c y+u w.
 \end{gathered}                                           \tag{8}
\]

`Vε(Qε)=2uQε` と `∂yAε+∂wBε=2u` は直接の恒等式である。
ε≠0 では fiber Hessian は階数2であり、平方完成後の binary q は
`−a²/(4ε)`、従って族 I に入る。
よって (7) の全支持は I と同じ MPS40 の和に含まれる。

## 4. Fiber Hessian の階数0、射影 Q の階数4

Q は `a(x)y+b(x)w+q(x)` の形である。
二つの斉次一次式 a,b が独立なら、fiber の定数変換と affine 平行移動で
`Q=y+xw` とできる。これは射影階数4である。
(2) の全係数を比較すると

\[
 A=(x+c)y+(c+1)xw-xg(x),\qquad
 B=dQ+g(x),\quad g\text{ は affine},                       \tag{9}
\]

c,d は定数となる。cQ,dQ を shear で除けば
`A=x(y+w−g), B=g` が全正常形である。

次の全パラメータの実変形がある。

\[
 \begin{gathered}
 Q_\varepsilon=y+xw+\varepsilon yw,\\
 A_\varepsilon=x(y+w-g)+\varepsilon y(w-g),\qquad
 B_\varepsilon=g(1+\varepsilon w).
 \end{gathered}                                           \tag{10}
\]

ε≠0 で `Y=y+x/ε,W=w+1/ε` とすれば `Qε=εYW−x/ε`。
さらに Aε から Qε 倍を除いた Vε の Y,W 方向係数は
`(x−1−εg)Y, εgW` で、和の係数は `x−1=f/x`。
従って (3) の全 Poisson 条件を満たす。
(10)は ε=0 でも正則なので、全正常形は族 II の閉包、従って MPS40 の和に入る。

## 5. 第2–4節で係数比較の例外を落としていない理由

Q の Hessian 階数は定数行列の階数なので、0,1,2で全てを尽くす。
階数1の a=0 と a≠0 をともに保持し、階数0では a,b の独立・従属を分けた。
階数2の q=0 は Q=yw が可約となるため、対象の既約性から正当に除かれる。
階数1の a=q=0 や binary q の平方の場合も、Q の階数≤2で可約なので対象外である。
次節は階数0で a,b が従属する全残余であり、省略しない。

## 6. 残る全源：Fiber Hessian 0、射影 Q 階数3

a,b が従属で零でない場合、fiber 座標を変えて
`Q=a_h(h,x)y+q(h,x)` とできる。既約性は `gcd(a_h,q)=1` と同値である。
三超平面のうち、a_h の根で消えない二因子の積を選ぶ。
その積の非零定数倍と q の差は a_h で割り切れるため、y の実線形平行移動で
その積を q としてよい。三超平面の並べ替えと尺度変更の後、全点を

\[
 Q=a(x)y+f(x),\qquad a(0)a(1)\ne0                         \tag{11}
\]

へ送れる。a は零でない affine 式で、定数の場合も含む。
(2) の全係数比較と Q 倍の shear により、正確に

\[
 \boxed{\quad
 A=y(a'y+f'),\qquad B=-a'yw+R(x,y),
 \qquad R\in\mathbb C[x,y]_{\le2}\quad}                    \tag{12}
\]

を得る。
`V(Q)=(a'y+f')Q`、`A_y+B_w=a'y+f'` は全係数の恒等式である。
従って a の二係数と R の六係数を持つ実際の既約な源であり、全磁気項を保持する。
射影退化式は

\[
 -h x(x-h)\{a_h(h,x)y+x(x-h)\}.
\]

この全源は第8節の支持計算と独立再構成により、新しい24次元全成分を与える。
第2–4節で扱った族への包含を仮定していない。
固定した三超平面のまま Q に任意の二乗項を加えるだけでは、
第2節の binary 根条件が保たれるとは限らない。
第8節では全実源の次元下界と全周囲の局所次元上界を合わせる。

## 7. 族 I と II の全 MPS 被覆

ここでは二つの実変形を独立に全パラメータで確認する。
結論は MPS40 の**和**への包含であり、T01 単独への包含ではない。
元の既約二次因子は T01 全閉包の五一次因子条件に反するため、両者を混同しない。

### 7.1 族 I

a(x)=Ax+B、`p0=−A−B,q0=−B`、`l_t=x−1+tz,F_t=xl_t` とし、

\[
 X_t=tp_0x\partial_x+
       \{(q_0-p_0)x-q_0+tq_0z\}\partial_z,
 \qquad U=y\partial_y-w\partial_w.
\]

次の実二次族を用いる。

\[
 P_t=F_t\partial_x\wedge\partial_z
       +(yw+1)\partial_y\wedge\partial_w+U\wedge X_t.     \tag{13}
\]

`X_t(x)=tp0x,X_t(l_t)=tq0l_t`、`div X_t=t(p0+q0)` より
X_t は F_t∂xz を保つ。U は (yw+1)∂yw を保ち、異なる二ブロックの場は可換。
従って (13) は全 t で Poisson であり、t=0 は族 I の任意の a に等しい。

t≠0 で (h,x,l_t,y,w) を座標に取ると、全ての交差項は logcanonical で、
唯一の追加項は `h²∂y∧∂w`。
h の重みを1、他を0とすれば、この追加項だけが消える。
極限の Pfaffian は非零な `±t x l_t y w` で、実五次式は五独立超平面の SNC 因子となる。
従って各 t≠0 の構造は MPS40 のいずれかに入り、和の閉性から t=0 の全族 I も入る。

### 7.2 族 II

a(x)=Ax+B、`r=−B,s=A+B,l=x−1+tz,F=xl` とする。
次の六括弧を順序 `(xy,xz,xw,yz,yw,zw)` で定める。

\[
 \left(
 \begin{array}{c}
 tsxy\\
 F\\
 t(1-s)xw\\
 (rl+sx)y\\
 F+[1+t(r-s)]yw\\
 -[(1-r)l+(1-s)x]w
 \end{array}\right).                                    \tag{14}
\]

t=0 で正確に族 II となり、Pfaffian は恒等的に `−F(F+yw)`。
t≠0 の (x,l,y,w) 座標では、定係数 logcanonical bivector の六係数は

\[
 (t,ts,t(1-s),-tr,-t(1-r),1+t(r-s))
\]

で、これに `xl∂y∧∂w` を加えたものが (14)である。
二つの互換条件は `t=ts+t(1−s)` と `−t=−tr−t(1−r)` であり恒等的に成立する。
追加項の自己 Jacobi も零なので全 (14) は Poisson である。
別確認として、その逆形式は

\[
 {1\over t}X\wedge L+
 (Y-rX-sL)\wedge\{d\log(F+yw)-X-L\},
 \quad X=d\log x,L=d\log l,Y=d\log y
\]

であり、全六 cofactor が (14)に一致する。
斉次座標 (h,x,l,y,w) に重み `(0,1,1,0,0)` を付けると追加項だけが消え、
logcanonical 部分の Pfaffian は `t x l y w≠0`。
従って実 SNC 極限を持ち、全 t≠0、さらに閉性から t=0 の全族 II が MPS40 の和に入る。

以上により、従属三超平面分岐のうち第2–4節で覆われない支持は第6節の全源へ絞られる。

## 8. 残余源の24次元全成分性と既知60との相違

本節の有限算術 producer は exact QQ で成功し、別の Fraction 疎多項式実装が
producer や SymPy を import せず独立再構成した。
形式切片と Krull 次元不等式の論理も、算術とは別の担当が独立に確認した。

### 8.1 実源の下界と全周囲の二法線

`a=a0+a1x` と R の六係数を全て自由にした源は、
`a0(a0+a1)≠0` 上の既約な affine 八変数源である。
PGL5 で飽和した射影像閉包を Z とする。
次の点 P* を用いる。

\[
 a=x+1,\qquad R=1+2x+3y+5x^2+7xy+11y^2.
\]

全パラメータの十個の strict Jacobi 恒等式と Pfaffian 恒等式を保存した。
P* の**全150係数・全350射影 Jacobi**に対する QQ 階数は

| 行列 | 階数 |
|---|---:|
| 全微分 D | 99 |
| Euler-vertical | 24 |
| 全GL・vertical・scalar | 46 |
| 上記と全八源微分 F | 49 |

であり、`D F=0` を全350行で確認している。
従って `dim Z≥49−24−1=24`、全射影接空間次元は `125−99=26`。
実族接49方向を、群接46方向だけに置き換えていない。
全微分の核は51次元であり、実族の49次元像を補う二つの literal bivector R0,R1 を保存した。

二法線上の Jacobi の三二次単項式を全 D 像で割った投影は、正確に

\[
 \langle n_0^2,n_0n_1,n_1^2\rangle                       \tag{17}
\]

を張る。非零99角 Jacobian 小行列式は `34626699791224012800`。
全 D を消す三つの疎な左余核 functional の支持は36,60,58項であり、
`Jac(n0R0+n1R1)` の像は、単項式順 `(n0²,n0n1,n1²)` で

\[
 \begin{pmatrix}
 13320&70200&-58320\\
 22770&69930&-210600\\
 -270&28890&119880
 \end{pmatrix}.                                         \tag{18}
\]

この三行は独立であり、(17) の全二次空間を得る。
functional の全係数と `λD=0`、(18)の全係数を report に保持する。

### 8.2 真の線形切片と Krull 上界

全150次元係数空間内で、実族接49方向を切る49本の線形式を取り、
その切片の接空間を99個の D 可逆方向と二法線で張る101次元空間とする。
選択した99個の全 Jacobi 式を形式陰関数定理で解くと、
残る**全て**の式は二変数の形式的方程式となる。
左余核 functional が D を消すため、その真正二次初期式は (17)を含む。
したがってこの切片の完成局所環は零次元である。

49本の線形切断に対する Krull の次元不等式から、150係数上の全射影 Jacobi 局所次元は高々49。
実際の既約像は微分階数49を持つので下界も49であり、射影化・vertical 商後は24となる。
Z を含む全 Poisson 成分 C は P* を含むため

\[
 24\le\dim Z\le\dim C\le24.
\]

従って Z=C は24次元の全空間既約成分である。
局所成分の一意性や正確な冪零厚みは、この推論に不要であり主張しない。
全源が Z に含まれることはその定義から成立するため、
唯一性を追加しなくても第6節の全支持がこの成分へ入る。

### 8.3 相違を支える閉条件

Z の退化五次式は、三つの一次因子が2次元以下を張る
`L0 L1 L2 Q` という条件を全閉包で満たす。
これは射影的な因子パラメータ空間の閉 rank 条件の積写像による像なので閉である。
零五次式を含む cone を取れば rank2 境界も保持できる。
この閉条件は被約 SNC 五次因子を一つも含まない。
三因子が重なれば非被約で、相異なれば三超平面が余次元2で交わり SNC に反するためである。
従って SNC 点を持つ全 MPS40 成分とは相異なる。

他の既知成分のうち射影次元24なのは CPLUS_EXT と新しい Tw0 成分だけである。
CPLUS_EXT の一般退化因子は一超平面と相異なる二つの既約二次因子であり、
三一次因子を持つ上記閉条件に反する。
Tw0 の全閉包は `LQ²` という閉非被約因子条件を満たすが、P* の実因子は被約である。
従って両者とも Z と異なる。残る既知成分とは次元で区別され、
六つの rank2 成分は閉 rank 条件によっても区別される。

独立算術の確認が完了したため、Z は既知60に加わる一成分となり、
個別成分数下界61、一般 rank4 成分下界55を与える。
この追加は全成分数の上界や他の因子型の全被覆を証明しない。

### 8.4 有限証拠と再現

* `code/classification_three_hyperplanes_quadric_dependent.py`：
  `--mode ranks` と `--mode normal` の一つの再利用可能な runner。
* `reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_20261107.json`：
  全源恒等式、全八微分、四行列の QQ 階数・小行列式・hash。
* `reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_NORMALS_20261107.json`：
  全四行列の hash 再現、二法線、選択99行列、全左余核 functional、二次初期式。
* `reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_REVIEW_20261107.json`：
  Fraction 疎整数実装による全 D99、実源49、核51、全八微分、二法線、
  三 functional の全150列零化の独立再構成。
  (18) の行列式は `−34390137600000` と独立に一致した。
* `reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_DEPENDENT_SOURCE_20261107.json`：
  最初の rank 証明書が参照した runner と本ノートの旧 bytes。
  mode 追加・本節統合の前の hash を保持し、旧証明書を書き換えていない。

階数実行は約1.98秒、二法線実行は約2.36秒で、どちらも active
`E_reduced_nonnormal` の事前 gate と60秒・1.5GiB制限の下で完了した。
normal report は現在の producer 全文も内包する。
旧 rank 実行を再現する場合は、一時 root 内に SOURCE report の保存 bytes を復元する。
独立算術は約0.52秒で完了した。
