# 定数modular source IIの全点は既NEW25に入る

2026-11-08、PRIMARY E、producer plane_collision_review。純手証明。独立査読前。
**f=u+zwの完全A8 source、全frame／全modular c／全Hの全点は、既閉成分 R4_POSITIVE_WEIGHT6_25 に入る。**
新成分増分は零。先のD98診断をcomponentと呼ばず、実既知親への全係数同定を与える。
新operator、ambient point、kernel、profile、GB、jetは実行していない。

## 1. 全八sourceとactual normalizer

`CONSTANT_MODULAR_DEGREE_TWO_WHOLE_REDUCTION.md` のsource IIを使う。

```
f=u+zw, L=αu+γz-βw,
Yz=zL+βu+λz+ρw+b,
Yw=-wL+γu+σz-λw+c,
Yu=f-βuw-γuz-ρw²-σz²-bw-cz.
```

P=∂t∧Y+Jac(f)、全八係数は自由。全値がPoissonでPf=f、κ=1。

actual t shear

```
t'=t-λu-cz+bw
```

はλ,b,cを同時に零へ送る。Jac(f)のHamiltonian of a linear formをLeibnizで戻せば、
Yzのz/constant項とYwのw/constant項がちょうど消え、Yuの-bw-czも消える。
逆は同じlinear formを引く。quadratic coefficientsは変えない。

次にactual paraboloid translation

```
z'=z+p, w'=w+q, u'=u-qz-pw-pq, t'=t
```

を使う。det一でfをexactに保存し、逆は(p,q)の符号反転。
全parameter作用は

```
α'=α, β'=β-αp, γ'=γ+αq,
ρ'=ρ-αp²+2βp, σ'=σ+αq²+2γq,
λ'=λ-2αpq-2γp+2βq,
b'=b+αp²q+γp²-2βpq-λp-ρq,
c'=c-αpq²-2γpq+βq²-σp+λq.
```

元の全括弧に可逆代入して得る全式であり、Y2だけの同値ではない。
α≠0でp=β/α,q=-γ/αを取り、β'=γ'=0となる。
残るλ',b',c'は上のt shearで除く。このとき

```
R=ρ+β²/α, S=σ-γ²/α
```

が残quadratic coefficientsである。

RS≠0のopenでは

```
u'=a*b*u, z'=a*z, w'=b*w, t'=a*b*t,
tensor scale=1/(a*b),  a*b=α, a/b=1/R
```

を使う。a²=α/R,b²=αRの有限平方根を選ぶだけで、全変換はactual GL/tensor scale。
結果は

```
α=ρ=1, β=γ=λ=b=c=0, σ=m=R*S.
```

従って全generic source IIは一つのmodulus m≠0の実正常形へ戻る。
αRS零は現段階で除外したままにせず、§4で全A8へ閉性を戻す。

## 2. 既positive-weight6全三parameter sourceのdense boundary

旧componentのfull raw affine sourceをそのまま使う。
`Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz` の
`full_affine_P` に n0=0,n1=s,n2=r を代入する。s≠0,r≠0を保持する。
source座標を(x,y,Z,W)と書くと全非零括弧は

```
Pxy=s/36-256W/3,
PxZ=(s/48-64W)*(r*y+x), PxW=8/9,
PyZ=-s²W/1536+sW²-sx²/64+64Wy-2/9,
PyW=4x/3, PZW=rxy+x²+2y/3.
```

全parameterでのPoisson性と、この全sourceの既閉25成分所属は旧独立証明を再使用する。
canonical homogeneous coefficientsをrawと取り違えず、上のaffine/raw表示から変換する。
modular vectorは正確に(s/48)∂Zである。

k=s²/221184+256/(27s)と置き、sourceからtargetへのactual affine mapを

```
t=4096Z/(9s),
u=-8y/9+k,
z=x+8W-s/384,
w=2x/3-16W/3+s/576,
tensor scale=-27/256
```

とする。全inverseは

```
x=z/2+3w/4,
y=(9/8)(k-u),
Z=(9s/4096)t,
W=s/3072+z/16-3w/32.
```

対応するhomogeneous GL5はhを固定し、detは-1048576/(243s)、非零。
新modular vectorは-∂t、従ってκ=1。

base (x,y,W)のJacobian potentialは

```
f0=2x²/3-8y/9+sW/36-128W²/3+256/(27s).
```

上のbase mapでexactにu+zwとなる。
base determinant=-256/27とtensor scale=-27/256の積は一なので、
base三括弧はexactにJac(u+zw)へ移る。
全四変数Pfaffianのscaleも(-27/256)² det(M)=-48/sで、同じfへ一致する。

## 3. 残る三括弧の全parameter同定

mapとtensor scaleのt係数の積は-48/s。
W0=W-s/3072とすると、Leibnizから

```
PZx=64W0*(r*y+x), PZW=r*x*y+x²+2y/3
```

を使う。xその他のspatial coordinateを反転していない。
これらを足すとtargetのYz,Ywは次のexact式になる。

```
Yz=-(384/s)*[z*(r*y+x)+2y/3],
Yw=(384/s)*[w*(r*y+x)+4y/9].
```

x=z/2+3w/4,y=-9u/8+9k/8を戻すと、source IIの係数は

```
α=432r/s, β=288/s, γ=-192/s,
λ=-432r*k/s, ρ=σ=0,
b=-288k/s, c=192k/s.
```

Yuも直接変換すると

```
Yu=f+(192z/s-288w/s)u-192kz/s+288kw/s
```

であり、同じ全source II式と一致する。
したがって全六raw括弧／全150 homogeneous係数のactual等式である。
Euler gaugeは零。選択Pfaffianや一つのcomponentだけを一致させたものではない。

この既sourceでは

```
R=192/(r*s), S=-256/(3*r*s),
m=R*S=-16384/(r²s²).
```

r,s非零を動かすとmはC*全体を走る。有限平方根の選択は複素数上で常に存在する。
従って§1の全generic source IIの実正常形は、既NEW25上のこの全族とactual PGL同値である。

## 4. 全zero補集合と元sourceへの復帰

fixed f=u+zwの全sourceはA8で、元tensorへのmapは全parameter上regularかつ非零。
既NEW25は閉PGL不変なので、その逆像は閉集合。
§1–3のαRS非零dense openがこの逆像に入るため、逆像は全A8である。
α=0、R=0、S=0、全lower/central/zero parameterを単位formulaへ代入したのではない。
元の全source自体のirreducibilityとregular mapを使っている。

さらに全base affine normalizer、modular cのframe、unique HのPGL frameを戻せる。
既NEW25のPGL不変性により、全source II飽和像閉包Z_IIも同じ親へ入る。
Z_IIの既exact dimension24は、この25成分の真の境界という結論と整合する。
旧D98／tangent27はその親を識別した根拠ではなく、今回の全actual mapとnormalizerが根拠である。

本稿はsource I/IIIをこの境界へ外挿しない。両者は別の既約24閉包であり、全E帰属は未解決。
新算術なし。新成分性／局所scheme／総成分上界を主張しない。旧reportsとglobal／main／Gitは変更しない。
