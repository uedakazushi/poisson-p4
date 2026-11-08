# 二つの E1DX weak wall：全 W六源を二つの実軌道へ還元

2026-11-08、PRIMARY E、producer 純証明・独立 review 前。
新 Jacobi／operator／rank／jet 計算は用いない。
入力は十八 module report の `E1DX_1`、τ=4/21,1/12 の二つの
**完全 affine W-kernel 六**と actual divisor の全 affine coefficient。
proper factor 未証明という先の weak-wall erratum はその時点の正しい限定であり、
本稿は実 normalizer により全点の認識義務を二つの actual points へ縮約する。
既 component への包含又は新 component はまだ主張しない。

## 1. 同じ literal Q と二つの canonical divisor

座標(a,b,c,z,w)、y=(a,b,z,w)とする。両 wall の全 tensor は
`Π=Qτ+∂c∧W`、W=W0τ+全十八の中の六次元 affine solution。

```
Qab=−(3τ+2)az+w², Qaz=0, Qaw=−3z²,
Qbz=(6τ+1)z², Qbw=9τzw, Qzw=0,
W0=(a²,−2ab,3τaz,(6τ−1)aw−3bz).
```

Euler normalization `Wz[z²]=Ww[w²]=0` を用いる。
元 full350 affine system と完全 kernel は report に保存済みであり、W六は homogeneous P2 の制限ではない。

particular solution の actual divisor を

`D0=K ab z³+L a w³z+M b w²z²−K c z⁴+E w⁵`

と書く。定数は次の通り。

| τ | K | L | M | E | η=1−3τ |
|---:|---:|---:|---:|---:|---:|
| 4/21 | 45/7 | 0 | 15/2 | 1/2 | 3/7 |
| 1/12 | 9/2 | 9/4 | −15 | −1 | 3/4 |

完全 W六の actual D image は正確に

```
D0+R(awz³+κw³z²)+A0 az⁴+B0 bz⁴+C0 z⁵+C1 wz⁴+C2 w²z³,
κ=7/6 (τ4/21), κ=−4/3 (τ1/12).
```

六方向は monomials awz³,az⁴,bz⁴,z⁵,wz⁴,w²z³ の係数から独立。
従って fixed Qτ／fixed Euler-normalization の W六から actual D への affine map は injective。
六次元 W の chosen point だけの式ではない。

## 2. 全 Qτ を保つ六次元実 unipotent 群

actual coordinate changes を

```
a'=a+A z,
b'=b+C z+(ηA/3)w,
c'=c+ell(a,b,z,w),
z'=z,w'=w
```

とする。ell は全 linear form、A,C は全 C。
全 map は determinant一、inverse は negative triangular substitution。

一般の a-shift A z、b-shift D w が Qab に加える central coefficient は
`(ηA−3D)z²`。D=ηA/3 により正確に零。
他の mixed entries はこの shear で変わらない。
c-shift は Q の core を変えず W のみを変える。
従って **Qτ は全群で literal に保たれる**。
normalweight零の P0 も変わらず、必要なら正の Euler redundancy で上の W gauge へ戻る。
元全 source／fullJacobi を保つ実 PGL 作用である。

## 3. actual D image 上の全 transitivity：global triangular inverse

inverse substitution を a→a−Az、b→b−Cz−ηAw/3、c→c+ell として書く。
particular D0 の new awz³ coefficient は `−KηA/3`。
w³z² coefficient は `−LA−MηA/3` なので、その比は

`κ=(L+Mη/3)/(Kη/3)`。

表の二 κ と一致する。
さらに w²z³ coefficient は `−M C`。
Kη と M は両 wall で非零定数なので、任意の R と C2 に対して A,C は一意に戻せる。
残 az⁴,bz⁴,z⁵,wz⁴ coefficient は、既知の A,C の polynomial を引いた後に
ell の a,b,z,w coefficients が各 −K 倍で独立に動かす。
ゆえに任意の六 parameter を global に消せ、逆も全 parameter 上の polynomial map。
generic dense orbit 又は小さい近傍だけの transitivity ではない。

群作用は Qτ/P0 を保ち fullJacobi source を保つので、像 tensor は同じ W六に入る。
§1 の actual D map の injectivity から、D を D0 に戻した後の tensor は
Euler redundancy を除いて report の particular tensor そのもの。
したがって **各 wall の全 W六 source は一つの actual projective PGL5 orbit**。
以前の cohom／first-symbol normalizers の逆も actual なので、元二つの全弱境界源も
この二 actual orbit の saturation に入る。

## 4. canonical points は実際に rank四・既約・本質五

particular W は元全三百五十式を満たす exact affine solution。
D0 は c一次、c係数−Kz⁴と残項の gcd は E w⁵ のため一。
Gauss lemma により C 上絶対既約、かつ非零なので tensor は generic rank四。
ゆえにこの二弱境界の rank四既約 source の存在は、未知としたままではなく確認出来る。

`∂cD0=−Kz⁴` から SingD0 は z零に入る。
そこで `∂wD0=5Ew⁴` なので wも零。
逆に z=w零では全 partial が零である。
従って reduced singular support は正確に一つの plane S={z=w=0}。

全五 partial derivatives は C-linear independent。
まず ∂z の ab z² monomial、次いで ∂w の w⁴、次に ∂a の b z³、
∂b の a z³、最後に ∂c の z⁴の係数を比較すればよい。
したがって D0 は essential-five であり、essential-three 八成分としての認識は不可能。

## 5. 次の正確な E 義務

全 source の E 認識は二つの particular tensor の既知 parent containment へ縮約された。
K22/nil22-long 既23成分、K4a smooth24親、他の既 MPS closure との
actual curve／actual linear normalization を優先して探す。
unique singular plane 又は normal cubic の型が他の一般元と違うことだけから
高次元 parent の boundary containment を排除しない。
新 componenthood、全 prior distinctness、orbit の正確 dimension 又は新成分数は本稿の結論ではない。
