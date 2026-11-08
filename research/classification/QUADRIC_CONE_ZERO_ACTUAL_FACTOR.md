# Quadric cone rank零の全 source：実五次は必ず一次因子を持つ

2026-11-08、PRIMARY `E_reduced_nonnormal`。既受理の全 Koszul/Euler・Ore source からの純証明。
**quadric cone 上で零となる指定全 source に、actual D が被約既約 quintic である rank四点は存在しない。**
旧 G像≤23 を E の全空性へ強める。新 CAS、operator、rank、minor、jet 又は parameter sample は使わない。
可約／非被約 D を持つ rank四構造全体の不存在は主張しない。

## 1. 全 source と全逆を保持した入口

`QUADRIC_CONE_ZERO_REVIEW.md` と正式 independent receipt の全源定理を入力とする。
座標 `(u0,u1,u2,ζ,w)`、q=u0u2−u1²、S=(w=q=0)。
ここで ζ は旧 note の vertex 座標 z であり、以下の base 座標 z とは別である。
全原形は Π=wL+qA、全 L五十／A十係数と Euler 冗長性を残す。
conormal 必要条件、全mod-w Jacobi、全 thirty-kernel と twenty-five Ore variables の可逆橋は既独立確認済み。
clean、ordinary node、vertical divisor 無し等の追加仮定は置かない。

旧記号 λ=A^{w,ζ}、base skew C に対し、λ=0 又は C=0 は既に H=(w=0) 因子を持つ。
従って既約 D の仮定の残余では λC非零。
全 normal shear、SO3 の非等方／等方両 orbit、tensor scalar と vertex scale の実逆により

`A=∂w∧∂ζ+C, C=∂0∧∂2` 又は `C=∂0∧∂1`

の二全型へ移る。
Δ非零を仮定して nil 型を落とさず、Δ零／C非零も全て後者へ戻す。
有限 complex root の選択を含む全逆は既受理 source theorem のものを使う。

chart w=1 で **projective Euler 補正まで含む全 affine tensor** は

`P=B(u)+(Y(u)+ζT(u))∧∂ζ`,

`B=qC+u∧X+D−u∧a, X=C*dq, T=t−ρu−X, t=a−d_vector`,

`Y=−u q−u(b·u)−M u−e`。

ρ は旧 Ore source の scalar c。b三、M九、e三の全十五 coefficient を保持する。
Poisson 条件は全源上で

`Jac(B)=0, Lie_TB=0, Lie_YB=Y∧T`

に同値である。Lie の符号は

`(Lie_YB)_ij=Σ_l[Y_l∂_lB_ij−B_lj∂_lY_i−B_il∂_lY_j]`

と固定する。以下はこの full equality の元一／二 component だけを使う必要条件であり、
chosen parameter、homogeneous particular 又は subset module へ制限しない。

## 2. Semisimple 型：全 ρ／全境界で不変 hyperplane

既受理の全 Jac(B) support は a=(p0,0,p1)、D=δC。
全 Lie_TB 条件は

`t=((ρ+1)p0,0,(1−ρ)p1), ρ(δ−p0p1)=0`。

全 parameter に対して実 affine translation

`x=u0−p0, y=u1, z=u2+p1, d=δ−p0p1`

を行う。これは homogeneous PGL では u0→u0−p0w、u2→u2+p1w、
逆は反対の translation であり、p0、p1 又は d を反転しない。
全 base tensor と T は正確に

`Bxy=−xy, Byz=−yz, Bxz=d−xz−y²`,

`T=−(ρ+1)x∂x−ρy∂y+(1−ρ)z∂z`。

g(x,z)=Y_y(x,0,z) と置く。
元 Y の fixed cubic と radial quadratic の y-component はそれぞれ−yq、−y(b·u)なので、
y=0 上で消える。従って **g=Ax+Bz+C は affine linear**。
全 M、e は自由なままで、この次数制約は source の構造だけから従う。

full `Lie_YB=Y∧T` の xy、yz 成分を y=0 へ制限する。
J=d−xz と書くと、正確に

`J*g_z=−(ρ+2)xg, J*g_x=(ρ−2)zg`。

第一式の x²,xz,x 係数と第二式の xz,z²,z 係数は

```
(ρ+2)A=(ρ−1)A=0,
(ρ+1)B=(ρ−2)B=0,
(ρ+2)C=(ρ−2)C=0.
```

差を取ると3A=3B=4C=0、従って g=0。
ρ、d を反転しておらず、ρ=0,±1,±2、d=0／非零の全補集合を同時に含む。
旧 generic matrix minor の exceptional values を個別 profile に戻す必要は無い。

よって Y_y は yで割れ、Bの y-row、T_y も yで割れる。
**full affine P の y-row 全体**が yで割れるので、H=(y=0) は Poisson invariant である。
core B だけの不変性ではない。

## 3. Nilpotent 型：一つの full Jacobi row で全 affine bad component が消える

全 source は a=(p0,p1,0)、D=δC。
全 Lie_TB 条件は

`t1=p1−ρp0/2, t2=0, ρp1=0, −4ρδ+ρp0²−2p0p1+2p1t0=0`。

全 parameter に対し

`x=u0, y=u1+p0/2, z=u2, p=p1, d=δ−p0²/4, k=t0−p0`

とする。唯一の分母二は constant unit、inverse は全 parameter 上の translation。
このとき正確に

```
Bxy=y²−p*x+d, Bxz=2yz, Byz=z(z+p),
T=(k−ρx+2y)∂x+(p−ρy+z)∂y−ρz∂z,
ρp=0, 2ρd=pk.
```

g(x,y)=Y_z(x,y,0) と置く。
元の radial cubic/quadratic は z=0 上で消えるので、ここも全 source に対し g=Ax+By+C。
full xz 成分を z=0 に制限すると

`(y²−p*x+d)*g_y=(k−ρx+4y)*g`。

y²、xy、y 係数を順に比較すると

`B=4B, 0=4A−ρB, 0=kB+4C`。

従って B=A=C=0。
この結論には ρp=0、2ρd=pk の反転も零分岐の削除も不要である。
全ρ=0、p=0、k=p0由来の壁、δ任意を含む。
よって full affine P の z-row は zで割れ、H=(z=0) は Poisson invariant。
旧 c=0 の二閉枝を一つの generic source に取り違えていない。

## 4. Full projective actual D の一次因子へ戻す

四次元 affine Poisson matrix で一つの row が hyperplane h=0 上で零なら、
同所の rank は高々二、全 symplectic Pfaffian は零。
従って full affine Pfaffian は hで割れる。
これは core B の三変数 Casimir だけの命題ではなく、
§2／§3が全 vertex entry `Y_i+ζT_i` まで含むために適用出来る。

w=1 の actual quintic は、Euler 補正した full affine matrix の Pfaffian に非零符号で一致する。
semisimple 型の h=y は元の homogeneous u1、nil 型の h=z は元の homogeneous u2である。
translation はこの bad coordinate を変えていない。
従って homogeneous actual D は、それぞれ u1、u2を因子に持つ。
同値に H∩(w≠0) の dense open で D=0 なので、H の閉性から全 homogeneous H 上で D=0。
Euler 項を戻しても actual D=E∧Π² は不変であり、full projective source への因子の復元は正当である。

元の全 normalizer の actual inverse を戻すと、この非零一次因子も非零一次因子へ戻る。
λ=0／C=0の先行 H=(w=0) 因子枝と合わせ、全 quadric-cone rank零 source を尽くす。
actual D が quintic で被約既約という対象仮定に矛盾するので、その whole branch は空である。

これは既知／新 component への membership 推測ではなく、全元式からの actual factor 恒等的必要条件である。
新しい成分数は零。quadric-cone rank二、他の導手曲面、可約 D の source は本結論へ混ぜない。
