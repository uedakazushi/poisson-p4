# rank一long未解決五源：真正F50切片と一回のquadratic余核支持判定

2026-11-08、PRIMARY E_multiple_hyperplane。**有限計算前の数学的契約。新点・再ambient profile・高次jetは作らない。** 対象は元指定点のM1,M11,B1,B3,V1だけ。各whole像は既約25次元、保存Dのrankは97,98,98,99,99で、全接dimensionが下界より大きい。

## 1. 実全sourceのcanonical parameter方向

全g零源の元六括弧(19)を使う。全二次Fの基底は1,z,w,z²,zw,w²。次の全canonical源と全GL5／Euler／scalar作用を保持する。

| source | canonical U,Hの全parameter |
| --- | --- |
| M1 | U=−z∂z+(v0+vz*z+vw*w)∂w、H=t z、F全六 |
| M11 | U=λz∂z+(−1−λ)w∂w、H=t zw、F全六 |
| B1 | U=−z∂z−w∂w/2、H=A z+B w²、F全六 |
| B3 | U=−z∂z+w∂w、H=A z+B z²w、F全六 |
| V1 | U=−z∂z、H=z(Aw²+Bw+C)、F全六 |

全parameterで(U+1)H=0、degree制限を保ち、全sourceの実Poisson方向である。元actual q-affine／固有frameの逆は既whole源に戻す。V1のB方向やM1の全v方向を捨てない。固定点、整数scale十二は既nine packetと完全に同一。

保存M51にこれらcanonical parameterの実微分を加えた行列をFとする。新しいDを構成せず、保存Dで**全DF=0、rank F=50**を確認する。50に届かなければそのcaseのこの経路を停止し、方向又は点を追加して追わない。

## 2. 真のnormalと実50線形cut

保存Dの全kernel Kはn=150−r次元、free coordinate rowでK_free=I。全Fのfree-row行列で50列を選び、その実像をF50とする。K-coordinateでF50を補う標準basisを選び、元150係数方向へKで戻したnormal Nはd=n−50=3,2,2,1,1本。

保存Dのr独立標準座標列をUとすれば、D_Uはinjectiveで、[F50,N]はker Dを尽くす。従ってUとNを合わせた100次元LはF50の真の補空間である。free-rowのn×n basis inverseから**実50線形cut**を明示し、cut F50=I、cut U=cut N=0を確認する。固定divisor源だけのnormalではなく、元150係数空間における切片である。

## 3. 全quadratic Jacobiと真正余核

全350元式Jのnormal quadraticはJ(Ni)、J(Ni+Nj)−J(Ni)−J(Nj)の全組で尽くす。係数行列Qのcolumn数はd(d+1)/2のみ。

保存nonzero minor A=D[I,J]を使い、

```
C=Q−D[:,J] A^{-1} Q[I,:]
```

を計算する。全normal coefficientを残し、selected二次式の都合のよい部分だけを使わない。選択されたCの独立rowにはfull Dを零化するfunctional

```
φ_i=e_i−D[i,J]A^{-1}e_I
```

を保存し、φ_i D=0を全150列で確認する。全C rowがこのbasisの線形結合であることも保存する。

selected r equationsのU微分は可逆なので、形式陰関数定理でU座標をnormal tから消せる。DN=0なのでU(t)は二次以上から始まり、φは全Dを消すから、上のquadratic formsは**真の切片idealの二次初期式**である。auxiliary方向の補正任意性や未解決P2を無視した余核ではない。

## 4. 一回だけの有限m-primary証明書

normal変数d≤3のhomogeneous quadratic idealI2に対し、**最初からN=d+1を固定**する。全quadratic generator×全degree N−2 monomialのMacaulay係数行列を作り、degree Nの全monomialがspanされるかを一回だけ確認する。

full row rankなら、各degree N monomialの明示QQ liftを保存する。従ってm^N⊂I2でradicalはnormal maximal ideal。full row rankでなければ非零left functionalを保存する。homogeneous m-primary quadratic idealならdegree二のregular sequence d本を内部に選べ、そのcomplete intersectionのsocle degreeはdなので必ずm^{d+1}⊂I2。従ってこの固定degree testはd≤3の全quadratic idealについて十分で、failure時にdegreeやjetを延ばす必要はない。

正dimensionのquadratic coneが残った場合、**元のPoisson schemeも正dimensionの余剰枝を持つとは結論しない**。より高い障害があり得るため、今回の最弱経路が非閉鎖という意味であり、別点／高次jetの反復は行わず数学的routeを変更する。

## 5. 成功時の結論は成分性まで

quadratic初期式のradicalが原点なら、真の100次元切片のKrull dimensionは零。実50線形cutから、元150係数Jacobi局所環Rはdim R≤50。一方全source ZのEuler逆像は50次元なので下界と一致し、各閉既約Zが25次元全成分となる。

**smooth性、局所成分の唯一性、全nonreduced横断代数又は正確nil指数は主張しない。** 例えばy(y−x²)=0でx方向をcutするとy²がm-primaryでも二枝がある。同じ理由で、今回の切片支持だけからuniqueを付加しない。D100だった四positive例のsmooth unique証明とは射程が異なる。
