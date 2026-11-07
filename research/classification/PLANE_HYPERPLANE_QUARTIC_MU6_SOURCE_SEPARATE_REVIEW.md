# μ6 marked quartic：全源・全配置・旧成分比較の別作者検証

2026-11-07、独立 reviewer `/root/quartic_review`。**PASS**。
対象は type6 nil22-image、base 括弧零、実被約五次 `D=H F`、`F` 絶対既約四次、
`μ(F)=6`、`#(H∩F)=2` の全 cotangent 源である。
入口 note の必要条件と全 marked list、全十一源及び二 E6 control を査読した。
全150/350 ambient 成分性は root の別作者計算
`CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU6_AMBIENT_REVIEW_20261107.json`
に分離する。本検証はその source・既約性・相違の仮定を供給し、ambient 計算を重複しない。
全 μ5、全 H+quartic 又は E 全体の尽くしはここから結論しない。

## 再現できる有限算術

`classification_marked_quartic_mu6_review.py` は producer/CAS を import しない。
既存の別作者 Fraction・疎多項式演算だけを使い、gate PASS 後 0.205 秒、18,640 KiB で完了した。
公開 report `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_MU6_SOURCE_REVIEW_20261107.json`
には実行 reviewer の全文も保存した。producer の可変現行 runner には依存しない。

全十三行で actual H→h 行列と逆、元 F の可逆 pullback、全28×18 gradient 行列、
全 kernel、非零小行列式を再構成した。degree6 行順は `h^(6-i-j)y^i w^j` であり、
source 列順の十八係数を全て保持した。十一候補及び非零 E6 control は rank16/kernel2。
E6 a0 transverse は rank15/kernel3 だが、**全 kernel anchor の任意の二本の determinant が零**。
一つの失敗 frame だけによる除外ではない。

全 gradient syzygy W は `W_h=h ℓ` を満たす。`W−ℓ E` の二成分が完全な h-Casimir anchor であり、
保存 U,V と一致する。det(U,V)=cF、c≠0、[U,V]=αU+βV を全係数で検証した。
全六磁気 R を保持し、更に**任意の GL2 frame 四係数**を独立不定元にした。

```
U'=aU+bV,  V'=cU+dV,
α'=dα−cβ, β'=−bα+aβ,
Π=∂x∧U'+∂z∧V'+(α'x+β'z+R)∂x∧∂z.
```

全十 strict Jacobi と従って全350 projective Jacobi が恒等零で、
実 Pfaffian は `(ad−bc)cF`。det frame=0 の補集合は rank4 の対象でなく Pf=0 として保持する。
全既約 rank4 frame はこの定数 GL2 source に尽くされる。全 six magnetic form を捨てていない。
発散 N は全 frame・全磁気で fiber 成分だけを持ち、その係数は base の一次式である。
従って **N²=0 は全源で構造的恒等式**であって、一点から開条件として推測したものではない。
canonical 十一行は rankN=2、h は双対像外、非零 E6 transverse control は rankN=1。

## 全六磁気の正則 inverse と一軌道閉包

q'=q+L、L の各成分は base の全一次式、という実 PGL shear に対し

```
δR=−αL1−βL2+U(L2)−V(L1)=T L.
```

全十二 frame で T は rank5。保存全6×6行列、小行列式、左 kernel v を独立再生した。
最初の v_j≠0 に対する projection `N(m)=(v m/v_j)e_j` を使い、保存5×5 minor の逆から
**全6×6の定数行列 B** を作って `T B=N−Id` を全係数で検証した。
よって未知磁気係数を一つも反転せず、全 R を `(I/v_j) mon_j`、I=v m へ送れる。
GL2 frame を変えても二つの L の定義域を adjugate で可逆に変えるだけで、rank と商次元は同じ。
I≠0 なら q を I/v_j で割り tensor を掛けて canonical R=mon_j へ移る。
I=0 は実多項式曲線 R=t mon_j の t=0 で保持する。
従って**各全六磁気 source の PGL 飽和閉包は canonical 一軌道閉包に等しい**。
零磁気境界を反転により捨てていない。

十一 underlying marked curve の連結 automorphism は零。
A6 は normalization の係数比較で r→ζr、ζ³=1、A4+A2 は異種二 cusp を固定して残 scale も固定、
3A2 は三 cusp の S3、非零 E6 は r→±r だけである。E6 a0 はこの有限性を仮定せず no-frame として除いた。
D の directional kernel は intrinsic な fiber 二平面、唯一の一次因子は H なので、
連結 tensor stabilizer は base 恒等の `q'=Aq+L` へ落ちる。
projective tensor の anchor 比較が A=σId、非零 magnetic 商比較が σ²=σ を与え、σ=1。
残る stabilizer は kerT の一次元だけ。従って canonical orbit は正確に23次元。
I=0 では fiber homothety が追加され orbit≤22 であり、非零商との区別を保持する。
この議論には ambient smoothness・局所成分一意性を用いない。

## 全 marked list の必要性を独立に読解した範囲

元の全 cotangent source と Saito determinant より D は free exponents(2,2)、τ(D)=12。
Milnor 積公式 `μ(D)=μ(F)+8−r`、四次の δ≤3 より r≤2、正規化は P1、μ(F)=5 又は6。
局所全枝数を保持した δ≤3 列挙は入口 note の九 profile を与える。
μ6 は A6、A4+A2、3A2、E6 の四つ。r2 は μ=τ=12 を強制する。
r1 と μ5 は別枝として保持され、ここで二点配置へ無断に移さない。

四つの全 normalization は次であり、target shear・P1 座標変更に使った係数の零補集合を読解した。

```
A6       (s4,s²t²,8s³t+t4)
A4+A2    (s4,s²t²,12st³−6t4)
3A2      (s²(s−t)²,t²(s−t)²,s²t²)
E6       (s³t,s4,t4+a s²t²),  a=0 又は a≠0→−6.
```

A6 の一般第三 section `a s³t+b st³+c t4` は c≠0。
A6 条件が b=0,a≠0 を強制するので非零係数だけで正規化できる。
A4+A2 は異種 cusp の joining section と tangent section を固定し、A2 条件が s³t 項を消す。
3A2 の三 joining line は cusp の非共線性と Bezout により上の三 section そのもの。
E6 は t-translation と target shear で st³ を消すが、残係数 a=0 は別 control に保持する。

一般 line の binary quartic pullback の二 root 条件を全て使う。
A6 は cusp を含む C=0 と含まない C≠0 を分け、後者で A=0 の infinity 単根を処理。
A,C≠0 の (3,1) は v=−u/3,u³=1、(2,2) は矛盾。
A4+A2 は C=0、A=0、A C≠0 を全て分ける。
A=0 の残 quadratic discriminant が B=−6C、cusp 回避の triple root 比較が u=4/3、
double-double は C≠0 と矛盾する。3A2 は S3 で全 cusp を覆い、
cusp 回避 double-double が唯一 H=h+y+w、triple-single は u²(u−1)²=0 により cusp へ戻る。
E6 は cusp 回避で B=0 の infinity 補集合も保持し、triple-single は a u²=−6、
double-double は B/C=a²/4。a=0 はいずれも解を持たず、残 smooth hyperflex は τ13 で対象外。
従って二・四・三・二、計十一の全 r2 marked orbit が得られる。

独立算術は四 parameter map の basepoint 不在、**全 degree3 implicit map の rank10**、
各 F の代入零を検証した。よって image は絶対既約四次で、parametrization は birational。
局所有限冪級数で A6 の first odd7、A4+A2 の odd5/3、3A2 の三 odd3、E6 の orders3/4 を再生した。
δ の和は3なので他の特異点・複数 image branch はない。
十一 H の pullback も報告した binary 因子と全係数で一致し、各接触 partition を直接再生する。
これは全 list 必要性の手証明と、個々の presence/型認識の算術を別に検証したものである。

E6 cusp transverse の全 marked family は一般四次の h³y 係数 b を保持する。
局所 hF の weights(3,4) の Milnor 代数は長さ11、top16。
b h4y が非零 socle、c h3y² は高 weight で零だから τ12 条件は b=0。
a0 no-frame と非零 a rank1 cofactor を上の全 kernel 計算で両方保持した。
非零 control は accepted type8 全支持定理により **MINIMAL_M** へ入る。
単なる type6 の一点との見かけ上の類似による認識ではない。

## 旧 printed source と positive 全源認識

旧六括弧は独立 checker に再入力して strict Jacobi、実 Euler/Pf、発散を直接展開した。

* ESSENTIAL5 の F は `h4+h²xz+x³y−x³z+x²z²/2−x²zw`。D=hF の gradient rank は5。
* MINIMAL_M は F=`h4+3h²x²+12x³z`、N=`6x∂y`。新十一 dense orbit の rank2 と異なる。
* A0 は F=`h4+h²xz−z4/4`、N=`z∂y+x∂w`。
* DPM は F=`h4+h²xz+p(x,z)`、同じ N。QQ(√3) の全厳密演算で確認した。
  p(t,1) と p(t,1)−t²/4 はそれぞれ triple+simple の四次で、triple root は
  `−3+√3` と `3+√3`。いずれも非零である。
  前者は H=h 上の A2 と contact3+1、後者は h²=−t/2 の二つの A2 を与える。
  δ和3により正確に3A2。新九中の3A2は H が二 cusp を結ぶか cusp を避けるので異なる。
* Ga9constant の受理済み整数 anchor 八係数括弧を直接再生すると、**D の gradient rank は5**。
  新成分全体は closed essential≤3 に入るので、この一 actual witness だけで相違する。
  cofactor-image-membership 自体の closed 性や dense source overlap を仮定する必要はない。

A0 と DPM の旧 printed tensor は base=(h,x,z)、fiber=(y,w) の全 cotangent 形である。
全磁気 T をここでも独立に解き、rank5・左 cokernel を再生した。
旧 actual R の quotient は、checker の basis に対してそれぞれ
**I(A0)=2**、**I(DPM)=177+100√3≠0**。旧 orbit dimension だけから非零性を推論する必要もない。

D5 cusp-tangent marked pair の全正規化では F=h²yw+a h4+d h y3+e y4、e≠0。
weights(3,2) の hF は長さ12/top12の Milnor 代数を持ち、d h²y³ が非零 socle である。
τ12 が d=0 を強制し、既約性が a≠0 を強制する。
y=βY,w=(a/β)W,β4=a/e で全 pair を `h²YW+h4+Y4` へ可逆に送る。
旧 A0 の a=1,e=−1/4 はこの unit chart に入り、全二次 Saito space は二次元なので、
旧 anchor はその全 GL2 frame、旧磁気は全 six source 内に入る。
I≠0 の上の実 normalizer により canonical orbit そのものに一致し、
**D5 cusp-tangent 全 source closure は旧 A0**。

同様に 3A2 cusp-tangent pair は上の全 normalization/S3 により一 orbit。
完全二次 Saito space・全六磁気・旧 I≠0 を合わせ、
**3A2 cusp-tangent 全 source closure は旧 DPM**。
Qcoker零から旧親 membership を推測したものではなく、actual full-source 逆変換による認識である。

## 全 prior111 と九候補36組

公開 prior111 は commit37d56e2afb3fd467401d395e39f9c34f32bae7b3、
blob dfc7c037b2050012fbe65c6c34b0e1ad53d3ff59 の全111行と byte hash を Git から照合した。
rank/dimension の異なる行、旧 proper multiplication/nonreduced closed locus、
旧 irreducible quintic witness、新 whole essential≤3 を個別に適用し、全111行の理由を保存した。
同23 MPS29 の actual full-period toric-SNC witness は既別作者確認済み
`CLASSIFICATION_ESSENTIAL_THREE_QUINTIC_FIRST_REVIEW_20261107.json`
の全29 named IDs と照合して再利用する。その実 D=κhxyzw、κ≠0 は essential rank5。
一般因子次数の名称だけからの比較ではない。旧五つの1+4成分は上記 actual witness で処理する。

九候補相互の全36組も保存した。異 curve profile は intrinsic cusp 型、
同 A6 は H 上の cusp 有無、同 A4+A2 は cusp 数及び A2 contact3/2、
同3A2は cusp数2/0で相違する。唯一の一次因子 H と essential三次元商は intrinsic なので、
PGL5 同値から marked PGL3 同値が必要である。
各 source は dense23 orbit を持つから、同じ23 component なら二つの dense orbit は交わらなければならず矛盾する。
これと root の独立九 ambient PASS を合わせれば、九つの成分採録を支持する。

**E6 smooth-flex は D100/two normals の未知として保持する。**
μ5及び新十三候補を全て既知とする推論はしない。
