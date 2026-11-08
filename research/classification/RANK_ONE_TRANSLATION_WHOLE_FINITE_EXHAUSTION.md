# rank一・c非imageの全long／short源：十七whole源への有限被覆

2026-11-08、PRIMARY E_multiple_hyperplane。**全源routingと三つの追加既知認識の純証明。全支持の完成はまだ主張しない。** 新ambient点、M1、minor又はjetを実行しない。long非定数gの一既実行point／既別作者receiptはそのまま保持する。

## 1. 全入口と現在閉じたg枝

`MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §5 の全六括弧(19)を用いる。X=∂x+x∂y、r=y−x²/2、q=(z,w)、d=div Uとして

```
U(g)=d*g, (U+1)H=F*g,
deg g≤1, deg U≤1, deg F≤2, deg H≤3.
```

全short源は `RANK_ONE_SHORT_TO_LONG_WHOLE_POLYNOMIAL_CURVE.md` のactual曲線で同じg分岐の全long源へ入る。従って以下のlong被覆でshortの補集合も全て保持できる。

* 非定数gの全源：新25成分 `RANK_ONE_LONG_NONCONSTANT_25_COMPONENT.md` へ帰属済み。
* g0≠0、H3=0：旧全源Gm×A5×A6から既存positive-weight NEW25へ帰属済み。
* g0≠0、H3≠0：`RANK_ONE_LONG_RESONANT_AND_BINARY_CUBIC_ENTRANCE.md` §2 のμ=1又は1/3の**二whole源**。全H≤2六係数、全translation、全非零H3及びg0を保持し、各像は25次元だが成分性／帰属は未解決。
* g=0：Fは**任意二次式六係数**、Hは非零で(U+1)H=0。以下で全(U,H)を十三whole源へ尽くす。

g零かつH零はPf零で一般rank四でない。非零定数Hは(U+1)H=H≠0なので起こらない。

## 2. g零の全十三源

以下のnormalizationはqのactual affine GL2だけを使い、元(x,y)を非線形変更して帰属を推論しない。Fは各normalization後も全二次式のまま保持する。Hの全表示係数は非零の場合だけunitと書き、零は後述の別行へ戻す。

### 2.1 五つの単項式型

| ID | H | 全U | 現在の帰属 |
| --- | --- | --- | --- |
| M1 | t z | −z∂z+v(z,w)∂w、全affine v | 未解決 |
| M2 | t z² | −z∂z/2+v(z,w)∂w、全affine v | nil31-third25、§4 |
| M11 | t zw | λz∂z+(−1−λ)w∂w、全λ | 未解決 |
| M3 | t z³ | −z∂z/3+v(z,w)∂w、全affine v | nil31-half25、§4 |
| M21 | t z²w | λz∂z+(−1−2λ)w∂w、全λ | nil31-DDH25、§5 |

t≠0。M1/M2/M3のvはconstant／cross／diagonalの全三係数を持ち、Jordan、零固有方向の非零translationを落とさない。M11/M21では二つの相異なるfactor idealがUで保たれるので、中心を取った後は表示のdiagonal Uが全解である。λ零／scalar値等も表示の全sourceに含む。

### 2.2 相異なる非零二固有値の全六binomial型

U=λz∂z+μw∂wとし、Hの二係数A,Bはどちらもunit。

| ID | (λ,μ) | H |
| --- | --- | --- |
| B1 | (−1,−1/2) | A z+B w² |
| B2 | (−1,−1/3) | A z+B w³ |
| B3 | (−1,1) | A z+B z²w |
| B4 | (−2/3,−1/3) | A zw+B w³ |
| B5 | (−1/2,−1/3) | A z²+B w³ |
| B6 | (−1/2,−1/4) | A z²+B zw² |

全q-affine frameと中心を戻したwhole源として扱い、A/Bを根拠なく一へ固定しない。全六源の全Fと帰属は現在未解決。A又はBが零なら§2.1の該当単項式型へ戻る。

### 2.3 二つの追加vertical源

* **V1（非零固有値と零固有値）**：U=−z∂z、H=z(Aw²+Bw+C)。新たに残す開はA(B²−4AC)≠0。A零、discriminant零、H零の全補集合は§3で単項式型へ戻す。全F六係数を保持する。全affine U/frameを戻した像は25次元であり、帰属は未解決。
* **SC（scalar cubic）**：U=−Euler(q−q0)/3、H=T(q−q0)、Tは任意三次binary式。新たに残す開はT squarefree。重根境界はM3/M21へ戻す。全F六係数を保持する。whole源とその25次元像の証明は前掲新三入口note §3。帰属は未解決。

これはfinite sourceの被覆であり、十三の像閉包が全て異なる、又は全て成分であるとはまだ言わない。

## 3. 全(U,H)を尽くす証明と零係数補集合

U1が可逆で相異なる固有値を持つ場合、全translationを中心へ移し、Uをdiagonalにする。Hの支持は

```
W={(i,j): i,j≥0, 1≤i+j≤3},  iλ+jμ=−1
```

に限られる。λμ(λ−μ)≠0の下では同時に許される格子点は高々二つ。九点Wの二点を選び二元一次式を解き、λμ=0又はλ=μを除くと、z/w交換を除いて§2.2の六組だけを得る。support一点なら、交換を除いて(1,0),(2,0),(1,1),(3,0),(2,1)の五単項式型となる。

U1=λI+N、λ≠0、N≠0がJordanなら中心を移し、degree kの固有値はkλ。固有値−1にはλ=−1/k、k=1,2,3だけが可能で、Nのbinary degree k上のkernelは一つのlinear formのk乗。従ってM1/M2/M3の全affine-v表示に既に入る。

U1=λI、λ≠0ならHは一つのdegree kで斉次、λ=−1/k。k一はM1。k二のsquare／相異二factorはM2／M11。k三のtriple／double+simple／相異三factorはM3／M21／SC。全binary式は複素数体上でこの有限factor分割に尽くされる。

U1の固有値がλ≠0と零の場合、非零方向のtranslationを消し、U=λz∂z+b∂wとする。

* b≠0では各z次数のpolynomial ODEからH=t z^k、λ=−1/k。M1/M2/M3の全vがこのbも保持する。
* b=0ではH=z^k p(w)、λ=−1/k、deg p≤3−k。k三はM3。k二でp定数ならM2、p一次unitならw平行移動でM21。k一でp定数／一次ならM1／M11。p二次でdiscriminant零ならH=z*(linear w)²なのでz/w交換後M21。p二次でdiscriminant非零だけがV1の追加開である。p零ならH零でrank四の外。

U1の全固有値が零なら、任意のconstant項を含むaffine Uは有限次数polynomial上でnilpotentである。従ってU+1は可逆で、非零Hは存在しない。この空枝も落としていない。

以上で全十三源と全zero/unit補集合が尽くされた。個別点や一般periodだけの支持推測ではない。

## 4. M2とM3の全六係数による既存成分認識

旧 `NIL31_IMAGE_ZERO_COMPONENTS.md` §1 のcoords(h,a,b,c,d)、全磁気quadratic A(a,h,d)、S=u_old*a+v_old*h+(τ_old−r_old)dを使う。元の全十変数で

```
Pab=r_old*a*c, Pac=r_old*h*a, Pad=0,
Pbc=A+(3r_old+τ_old+1)(hb−c²/2),
Pbd=g_old+cS, Pcd=hS.
```

実PGL置換(long x,y,z,w)=(c,b,a,d)の後、次の**非零tensor scalar**を掛けるだけでよく、Euler補正は不要である。

* **M2**：third、r_old=−1/3、g_old=q_old*h*a、tensor −3/2。得られるlong sourceは
  
  `H=(3q_old/4)z², F=−3A/2, U=−z∂z/2−(3S/2)∂w`。
  
  任意のM2係数t,v=v0+v1z+v2w,Fに対し、逆は
  
  `q_old=4t/3, A=−2F/3, u_old=−2v1/3, v_old=−2v0/3, τ_old=−2v2/3−1/3`。
* **M3**：half、r_old=−1/2、g_old=p_old*a²、tensor −2/3。得られるlong sourceは
  
  `H=(2p_old/9)z³, F=−2A/3, U=−z∂z/3−(2S/3)∂w`。
  
  任意のM3係数に対し、逆は
  
  `p_old=9t/2, A=−3F/2, u_old=−3v1/2, v_old=−3v0/2, τ_old=−3v2/2−1/2`。

各式は全六括弧へ直接代入して一致する。逆の分母は非零有理定数だけで、全F六係数と全affine v三係数を保持する。t非零は旧p/q非零のrank四条件と一致する。一般affine linear formの平方／立方はq-affine逆を戻して全て含まれる。従ってM2は既nil31-third25、M3は既nil31-half25の全閉成分に属する。

## 5. M21の全sourceは既nil31-DDH25へ入る

H=t z²w、t≠0、U=λz∂z+μw∂w、μ=−1−2λ。斉次liftのz-Hamiltonian cofactorは全sourceで

```
Yz=P#dz/z=λh∂x+(λx+t z)∂y.
```

P+E∧Yzはz-Casimir liftとなり、その全divergenceは

```
Xnew=(1+5λ)h∂x+[(1+5λ)x+5t z]∂y.
```

a=1+5λ≠0では、z=1の新affine chartにおける線形部はnilpotent Jordan(3,1)、定数部はそのimage内。x′=x+5t/aで定数部を実affine translationにより消せる。その後は既 `MODULAR_NIL31_COMPONENT.md` の全15sourceへactual linear normalizationで入る。

新affine Pfaffianは非零定数倍h²w。chain末端一次式im(A*)²はhで、wは独立なinvariant covector。旧n7零の全rank四枝のPfaffianは末端一次式だけのpolynomialであるから、ここではn7≠0が強制される。従って元の全Fも既nil31-DDH25成分へ入る。

λ=−1/5を除去していない。全M21のλ,t,F源は既約で、a≠0は稠密開であり、元の全六係数はλの多項式。既成分の閉性によりλ=−1/5の全Fも含む。q-affineの逆で全double+simple binary-cubic monomial源へ戻る。これは旧成分への一点の接触やβ上界ではない。

## 6. 全statusと残る十二whole源

合計は、非定数g一、定数g/H3零一、定数g/H3非零二、g零十三の**十七whole源**。このうち非定数g、定数g/H3零、M2/M3/M21の五源を既成分又は新受理成分へ全点帰属した。

残るのは、定数gの二共鳴源、g零のM1/M11、B1–B6、V1、SC、計十二whole源。各源はfullF／H／translationと全元frameを持つ。この有限リストを、三つの先行candidate入口だけに縮めて他のsourceを落とさない。残りの成分性、実親、全Eの完成は未証明である。
