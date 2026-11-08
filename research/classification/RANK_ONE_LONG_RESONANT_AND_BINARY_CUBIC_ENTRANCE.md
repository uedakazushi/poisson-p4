# rank一longの残余：二つの定数g共鳴源とg零scalar-cubic源の全25次元像

2026-11-08、PRIMARY E_multiple_hyperplane。**純全源入口。以下の三閉族の全空間成分性・成分数は未証明で、新ambient点は選択も実行もしていない。** 既知源への包含を源次元やβ上界だけで推論しない。

## 1. 共通の係数逆と飽和次元

`MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §5 のlong全源(19)を使う。

```
X=∂x+x∂y, q=(z,w), r=y−x²/2,
Pxy=−F−(1+d)r, Pxz=u, Pxw=v,
Pyz=xu+r*g_w+H_w,
Pyw=xv−r*g_z−H_z, Pzw=g,
U=u∂z+v∂w, d=div U,
U(g)=d*g, (U+1)H=F*g.
```

全degreeはg≤1,U≤1,F≤2,H≤3。Jacobi・divergence・Pfは前の全源証明と同じく

```
div P=X, Pf P=−(rg+H).
```

係数からg,U,FとHの両偏微分が復元され、(U+1)H=Fgが定数差も零にするから、全parameter mapはinjectiveである。

`RANK_ONE_LONG_NONCONSTANT_COMPONENT_ENTRANCE.md` §§3–4 の**全Aff4 stabilizer十次元**はgの型に依らない。以下の各SがStab(X)不変で既約12次元なら、Aff4飽和は正確に22次元。非零scalar stable、一般Dの唯一double-Hにより、全PGL飽和射影像閉包は正確に25次元となる。同次元の根拠を新しい標本から得るのではない。

## 2. 非零定数g・三次共鳴の全二源

g=g0≠0ならd=0、F=(U+1)H/g0。H3≠0では(U1+1)H3=0。trace零二変数場の四つの三次重みから、U1は相異なる固有値±μを持ち、μ=1又は1/3（符号は固有方向交換で同一）だけである。nilpotent又は零固有値ではH3=0なので、ここへ未確認に混ぜない。

各μでUはtrace零semisimple GL2軌道二次元と任意translation二次元からなる既約四次元源。U1は可逆なので中心q0が一意。q−q0の三次式上の固有値−1部分は一次元である。

```
μ=1:   U1=−z∂z+w∂w,   H3=t*z²*w,
μ=1/3: U1=−z∂z/3+w∂w/3, H3=t*z³,
```

ただし表示は固有frame／中心で行い、元全U、全translationと非零tを保持する。Hの二次以下の六係数は**自由**であり、Fだけが一意に決まる。従って各whole sourceは既約で

```
dim Sμ=4(U)+1(H3)+6(H≤2)+1(g0)=12.
```

sourceのStab(X)不変性も全域で成り立つ。qのaffine変更はU1を共役し、rへのaffine q式の加算はHへg0倍のaffine式を加えるだけなので、μと非零H3を保つ。

scalar stabilityは明示できる。aλ:(x,y,q)↦(λx,λ²y,q)とtensor λ^{-1}の合成は

```
g′=g/λ, U′=U, H′=λH, F′=λ²F
```

を与える。μも不変である。

全H3非零で実因子は

```
D=−h²*C,
C=g0*h²*y−(g0/2)h*x²+H3(z,w)+h*H2(z,w)+h²*H1(z,w)+h³H0.
```

Cはyについてprimitiveな一次式である。leading coefficient g0h²と定数部は、H3≠0によりgcd一なのでCは絶対既約。従ってhは唯一の重複因子、重複度二であり、射影像閉包Zμは正確に25次元となる。H3零境界は既存NEW25の源だが、それを理由にH3非零源をNEW25へ含めない。

## 3. g零・scalar U・全squarefree binary cubicの源

g=0、q0任意、

```
U=−[(z−z0)∂z+(w−w0)∂w]/3,
H=T(z−z0,w−w0), T任意の三次斉次binary式,
F任意の二次以下の式
```

は全域で(U+1)H=0を満たす。Fは全六係数を保持する。中心二、T四、F六なので、sourceSBは既約12次元で、係数逆から像も正確に12次元。Tのsquarefree非零開で

```
D=−h²*T(z−z0h,w−w0h)
```

となり、h²と相異なる三つの一次因子、計四つの異なる超平面を持つ。hが唯一の重複因子である。三つの単純因子のspanは二、hを加えたspanは三。

Stab(X)のq-affine変更はscalar Uの中心とbinary式を移し、rへのaffine式の加算はFを別の許容二次式へ移すだけ。全SBを保つ。上記aλ/tensor作用はg零、U不変、H′=λH,F′=λ²Fで同じ源内にある。従って全飽和射影像閉包ZBも正確に25次元である。

squarefreeは次元／唯一Hを確認する非空開であって、Tの重根境界を除去する仮定ではない。全Tが同じ既約源の閉境界として保持される。

## 4. 既知129及び非定数g閉族との分離

prior129 projectionは `rank_one_long_nonconstant_prior129_20261108.json` に固定してある。rank二全成分とdimension24以下の全成分は前と同じ理由で各25次元像を含まない。

Zμの一般Dは2Hと絶対既約cubicであり、CA/CB及びnil31三25成分のall-linear閉条件、positiveNEW25のtriplelinear閉条件の外である。T01/T39/T40への包含も、同次元なら一致する一方でZμは全てmultiple-Hの閉条件内、各MPS25には被約Dの開があるため不可能。

ZBの一般Dは四つの相異なる一次因子で重複度2,1,1,1。CAとpositiveNEW25は三乗因子を必要とし、CB及びnil31-DDHは二つの二乗一次因子を必要とする。half/thirdも三乗因子を持つ。これらはいずれも適切なprojective multiplicationのproper像という閉条件なので、generic ZBはその外。MPS25への包含は同じmultiple-H／同次元の矛盾で排除する。

さらに、非定数gの閉族Zncは全rank四部分のDがH²*L*Qというclosed factor incidenceに属す。絶対既約cubicを持つZμはこの外。ZBはこのfactor incidenceには属すが、もしZB⊂Zncなら同次元で等しいはず。一方ZBは全D all-linearの閉条件内、Zncはgeneric rank五Qを持つため一致しない。

μ=1と1/3の二像は、一般点で互いに同値にならない。唯一double-Hを固定したprojective同値はaffine同値とtensor scalarになり、modular Xを正規化したその全群はStab(X)と上のaλ/tensor作用で生成される。前者はU1を共役し、後者はUを保つ。従って固有値集合{−μ,μ}は全域不変である。各像が同じ25次元閉既約族に稠密なら、二つのconstructible像の非空開は交わるはずなので矛盾する。

以上は三つの**相異なる25次元像閉包**及び全priorからの分離まで。未知のより大きな全空間成分がこれらを含む可能性は、全ambient局所証明前には排除していない。

## 5. 残余と次の最弱判定

定数gのH3非零全二源は§2で尽くした。shortの同源は既actual曲線でlongへ入る。g零の§3は全incidenceの一枝であって、任意のU,Hをscalarへ正規化してはいない。一般の固有値共鳴、二つの単項式共鳴、零固有方向、Jordan境界、全Fは別残余として保持する。

各Zの全成分性には、もし次に許可されれば、各全源内の一つの構造的に許容される固定有理点でD100を得ることが最弱十分である。現在はその点も新計算も用意しておらず、有限contract／独立source reviewの前に追加profileを始めない。
