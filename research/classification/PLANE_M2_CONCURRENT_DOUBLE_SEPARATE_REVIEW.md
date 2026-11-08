# concurrent 三直線・m二 double：全十七族と全複素壁の独立確認

2026-11-08、PRIMARY G、別作者 `/root/quartic_review`。**全源 G 上界23を独立受理。**
全 Poisson 分類の完成又は下記三残源の E 帰属・存在は主張しない。

producer 凍結 manifest は `inputs/plane_m2_concurrent_double_manifest_20261108.json`
(`6bce852ec8812e0cb9f861002293f9bda18fd1081fa4089f453e2d0ba139aac1`)、
report は `CLASSIFICATION_PLANE_M2_CONCURRENT_DOUBLE_20261108.json`
(`5b39d3b708e4655453bef40f6ff61c9c14424894261e0a0b7605ff6cef0f08a3`)、
実行 producer bytes は同 report 内の `448c081c5110dc7f99e53f29631c4d655d16b29cd5597a35a66e490b601e70b8`。
現在の共通 runner を import 又は再実行していない。

## 1. 全 source の入口と正常化の逆

仮定は rank 四 projective Poisson、被約既約 actual quintic、restriction rank 二 plane、
被約 concurrent cubic f=xy(x−y)、normal 位数二、generic double normal quadratic。
primitive factorization G2=h*ell² の全 degree は0/1である。
degree 一の moving root は既全 connection 入口と split-q の actual Pf coefficient 障害で空。
degree 零では constant normal GL2 により ell=u、全 marked h は f、L²M、L³。
三 concurrent components の全 marked 選択は actual base PGL により戻す。

原全 normalform note `PLANE_M2_DOUBLE_ROOT_CONSTANT_NORMAL_FORMS.md` §§4–7 を独立読解した。
bad u-row v²C の base x/y components (c,d) には、記載の四必要式がある。
c のみ非零なら r=−6、d のみなら s=−6、両方非零なら determinant 6(6−r−s)=0 が必要。
全 Δ 型 (r,s)=(0,0),(0,3),(−3,6) はいずれにも入らず、全 period で C∥∂z である。
off-field、C の零/nonzero 境界、normal K の零境界を chosen slice にしていない。

h=f の全 off-field Z=(cx+dy)∂z と全 X=Y=(ax+by)∂z を保持する。
det(E,X,Z)=Zf=0 より actualG2=3kfu²、従って k!=0。
h=L²M では全 Z=A0*y∂z+B0*diag(−1,2,−1)+C0*x∂z を戻す。
actual Pf の全係数から B0!=0、a+b+4!=0、k=B0(b+4)。
normal shear が A0 を消し、relative normal scale が B0 を一にする。
全 inverse はそれぞれ逆 shear と B0 倍であり、c=C0/B0 の零/非零を保持する。
actual diagonal weights (0,0,1,0,1) と tensor scale で c だけが非zero character で動くため、
c=0/1 が全 source を覆う。

h=L³ では target 非零から原 A0!=0、a!=4。
残式 A0(−b−8)+3E0(−a−b−4)=0 は a+b+4=0 の chart を矛盾で除き、
A0=1 後 E=−(b+8)/(3(a+b+4)) を与える。
base shear z'=z+p x+q y の C0/D0 への作用は

    (C0,D0) -> (C0+3Ep+q, D0+6Eq).

E!=0 で determinant 18E² は非零。
E=0 は b=−8 の whole 境界で、base q-shear と normal shear Δ=(−3x+6y)∂z が独立に C0,D0 を消す。
この閉補集合を E で割って落とさない。
全 connection の非対数 closed-form module は原 note の全 Z に残されている。

## 2. equal-diagonal nilpotent off-field の全空性

`PLANE_M2_DOUBLE_EQUAL_DIAGONAL_NIL_OFF_EXCLUSION.md`
(`cc1c403472880b2a93fc746663d052acfed29ff399eb3f70523aa0ad740546a6`) を独立再導出した。
全 raw36 P1 と全9 P2 に対して row divergence は Ru=0、Rv=−2ku。
canonical strict lift Π=P+E∧R/5 と actual Casimir F を使うと

    M=Π#d(v/u)=P#d(v/u)−(2k/5)E,
    M(u,v)=(3k/5)(u,v),
    M(b)=(Z−2kI/5)b+B(u,v),   M(F)=0.

B は全 higher mixed source から来る base-independent rational homogeneous degree 一 vector。
Z nilpotent、k!=0 より T=−(kI−Z)^−1B は全 parameter で定義され、b'=b+T によって translation が消える。
逆 b=b'−T は fraction field C(u,v) 上の affine map であり、global PGL と呼ばない。
base degree d、normal total degree m=5−d の summand に対する operator は
Z|Sym^d+k(m−2)I。m!=2 で可逆なので、全 F は exactly u²f(b')、全体 scalar を除く。
normal order 二の leading は T の degree 一性で変わらない。
reducible f の全 factors はこの affine change 後も正の base degree を持つ。
従って fraction field 上の分解は Gauss lemma により実 polynomial F の既約性に反する。
これで tangent h=f と concurrent h=f の全源が空。
この証明に bad C の nonzero normalizer は不要であり、旧 tangent Z0 の上界22を空性へ強化する。

## 3. 全二変数 source を十七 a-lines へ制限する証明

一般 Δ を戻すと、同じ ratio-Hamiltonian の base field は
A(q)=Z−qΔ、q=v/u invariant、translation equation は (kI−A(q))T=−B。
base degree d=5−m の operator は

    A(q)|Sym^d+k(m−2)I.

この式は全 P1/P2 を保持した rational identity である。
translation が存在し m=3/4/5 の全 operator が可逆なら、全 F=u²h(b') となり reducible。
従って translation 又はその spectral wall の全和だけが必要残余。
normal order 二から元 m=0/1 summands は無く、degree 一 translation 後も生じない。

h21 では A(q) は全 q,c で triangular、eigenvalues (−1,2,−1)。
Sym^d spectrum 3j−d より、m=3/4/5 の全 k は
{2,−1,−4} / {1/2,−1} / {0}。
translation determinant (k+1)²(k−2) の失敗も同じ和に入る。
k=b+4 より全 b は −2,−5,−8,−7/2,−4、c=0/1 と合わせて十 whole a-lines。

h3 では全 q の eigenvalues (0,3E,−3E)、k=E(b+8)。
E!=0 で全 spectral walls は b+8∈{−6,−3,−3/2,0,3/2,3,6}。
translation failure k(k−3E)(k+3E)=0 もこの和に入る。
E=0 の b=−8 は原 source へ直接戻す。
従って七 b 値 −14,−11,−19/2,−8,−13/2,−5,−2 が全てである。
全 a の source domains と leading/order-change complements をそれぞれ保持する。
これは normalizer で period を一点へ置換した有限化ではない。

## 4. 独立全行列・全複素壁の検証

新 checker `classification_plane_m2_concurrent_double_review.py` は、成功済み独立 sparse/Fraction helpers のみを import。
producer/CAS は import しない。gate PASS 後の一回の run は **16.787秒 / 95,692KiB** で全 PASS。
新 report は `CLASSIFICATION_PLANE_M2_CONCURRENT_DOUBLE_SEPARATE_REVIEW_20261108.json`。

全17 generic QQ(a) source の full150 P0、full350 Jacobi、actualG2、350×30/350×9 operators、
150×39 Euler/Levi/scale matrix、complete kernels と非零 minors、全39 actualD linear terms を別に再構成した。
full30+Euler6=raw36、P2全9 と Euler rank24 を再確認した。
普遍39 quadratic table は旧独立再構成済み tangent packet との全データ一致で再利用し、旧 suite は再実行しない。

全 generic determinant numerator/denominator と complete-kernel numerator/denominator が
記載の全 wall polynomials の積で割り尽くされることを exact polynomial division で確認した。
全 generic operators の pole はこの和又は明示 source-domain 除外点だけである。
wall と除外点の和がその全 factor 和に一致する。
全84 allowed wall algebras を原 source/operators へ specialize して全証明書と照合した。
うち五つは非有理二次/四次 quotient algebras。

数体部分は producer と異なる multiplication-matrix inversion により QQ[a]/g を実装した。
選択 minor は quotient 内で **unit inverse** を直接確認し、kernel の free-coordinate submatrix は identity。
従って g の一 complex embedding だけの nonzero 確認ではない。
g の既約性 oracle を使わなくても、全 complex roots に同じ rank/kernel upper が成り立つ。
複素壁を有理点だけで置換していない。

## 5. actual u 因子と正しい全像上界

全 generic17 fields 及び81 wall algebras では、complete K1 の全 bad u-row maps が零。
P0 の u-row は u 倍、全 P2 に normal columns は無いので、全 full solution の u-row が u 倍となる。
actual D は u 因子を持ち、指定既約 source は空である。
K1 の一般元だけの因子又は chosen P1 slice ではない。
残りは次の三 fixed P0 の全 higher source だけである。

| whole family / a | K1 | K2 | filtered group fiber | projective image upper |
|---|---:|---:|---:|---:|
| h21,b=−2,c=0 / −4/11 | 9 | 0 | 10 | 23 |
| h3,b=−2 / 5/2 | 8 | 0 | 9 | 23 |
| h3,b=−2 / −1/2 | 8 | 0 | 9 | 23 |

全 field で normalizer rank を exact に確認した。
Euler rank24、effective projective Levi dim12 から graded stabilizer は s=37−rankG。
全六 base←normal shears の graph fiber が加わり filtered fiber は s+6。
各 fixed P1 に対する P2 は、空又は complete K2 torsor の部分集合である。
従って全非線形 compatibility を満たす部分の upper は

    24−(s+6)+parameter_dimension+dim K1+dim K2.

残三点は parameter dimension 零なので23。
generic 粗線形上界25/28等を実成分の存在とは解釈せず、全 u 因子による空性を使う。
これに全 pureproof complements を合わせて、**concurrent reduced R の全 m二 normal-double 源は G 像23以下**。
残三源の存在、空性の追加強化又は E 帰属、nonreduced R、m三以上及び他 surfaces は未結論である。
