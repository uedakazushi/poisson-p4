# constant modular source IIIの全十係数：B1成分へのactual極限

2026-11-08、PRIMARY `E_multiple_hyperplane`、担当 `/root/plane_collision_review`。
**純証明。新Jacobi/operator、点、minor、profile、jet又はGBは実行しない。**
全十係数源IIIを、既成分であるrank一long B1のactual PGL極限へ送る。
source IIIを新24成分として数えない。独立査読／root採録はこのproducerとは別である。

## 1. 全源とgeneric normalizer

`CONSTANT_MODULAR_DEGREE_TWO_WHOLE_REDUCTION.md` の必要十分表示を使う。
座標は(t,u,z,w)、f=u+z²、Jacの順序は(u,z,w)。全源は

```
Ptu=f−2zB, Ptz=B, Ptw=(2βz−γ)w+K2(u,z),
Puz=0, Puw=−2z, Pzw=1,
B=βu+γz+δw+e.
```

β,γ,δ,eとK2の全六係数を任意に取るA10である。
divP=−∂t、Yf=f、divY=1、PfP=fが全域で成立する。
δ≠0で次のactual affine変更とtensor scalarを取る：

```
T=δt, U=u, Z=z, W=βu+γz+δw+e, tensor scalar 1/δ.
```

全5×5斉次行列式はδ²、逆はt=T/δ、w=(W−βU−γZ−e)/δ。
得られる源は再び(T,U,Z,W)を(t,u,z,w)と書いて

```
Ptu=u+z²−2zw, Ptz=w, Ptw=Ktilde(u,z),
Puz=0, Puw=−2z, Pzw=1,
Ktilde=β(u+z²)+(γ−2βz)(βu+γz+e)+δK2.
```

これは全六括弧への直接代入であり、特にWのt括弧はY(B)を計算したもの。
K2→Ktildeはδ倍と固定二次式の加算なので全六係数を保ち、
逆はK2=[Ktilde−β(u+z²)−(γ−2βz)(βu+γz+e)]/δ。
以降はK=Ktildeを任意二次式として保持する。
δ=0は最後に同じA10の閉性で戻し、別の未知fiberへ一般点定理を外挿しない。

## 2. 既B1成分の正しい全source

旧座標を(X,Y,U,Z)とする。rank一longのg=0、
planar場−U∂U−Z∂Z/2、H=U+Z²のB1全源は

```
PXY=−F(U,Z)+Y/2−X²/4,
PXU=−U, PXZ=−Z/2,
PYU=−XU+2Z, PYZ=−XZ/2−1, PUZ=0.
```

Fは全二次式六係数である。divは∂X+X∂Y、Pfは−(U+Z²)。
これは原六括弧(19)へd=−3/2を代入した全sourceそのもの。
`RANK_ONE_FIVE_QUADRATIC_SUPPORT_SEPARATE_REVIEW.md` は、この全B1 sourceの
PGL飽和閉包が25次元全Poisson成分であることを独立認証している。
一つのwitnessから全Fの所属を推測していない。

## 3. 全Kを保った明示的PGL曲線

ε≠0で a=1/(8ε³)、b=1/ε とし、上の旧B1へ

```
Fε(U,Z)=K(a²U−2abZ−b², aZ+b)/(8ε³)
```

を代入する。任意KについてFεはU,Zの二次以下の式で、係数はLaurent多項式。
B1の全F六空間に属するため、係数を特殊な部分族へ制限していない。
次のactual affine変更を行う。tensor scaleは1、Euler補正は零である。

```
z=aZ+b,
u=a²U−2abZ−b²,
w=−2εX+z/2,
t=−X−4ε²Y+u/8.
```

行列式は−8ε³a³=−1/(64ε⁶)。両方向の逆を明示すると

```
X=(z−2w)/(4ε),
Y=(−t−X+u/8)/(4ε²),
Z=(z−b)/a,
U=(u+2bz−b²)/a².
```

全式はε≠0で可逆なPGL5変更であり、非線形DarbouxをPGLと呼んでいない。
H=U+Z²は正確にf/a²へ移る。

## 4. 六括弧の全係数恒等式

変換後は全て次のε多項式となる：

```
Ptu=u+z²−2zw+εu(z−2w),
Ptz=w+ε(z²/2−zw),
Ptw=K(u,z)+ε[t+(3z²−8zw+4w²−3u)/8],
Puz=0,
Puw=−2z−2εu,
Pzw=1−εz.
```

手展開を追えるようにcross項と最後の項を示す。
D0=−u∂u−z∂z/2、H0=2z∂u−∂z と置くと、essential座標の変更により
old XのcrossはD0−bH0/2、old YのcrossはX(D0−bH0/2)+aH0となる。
したがって新tと新wのcrossは厳密に

```
Pt,ess=[−1+ε(2w−z)]D0+(z/2−w)H0,
Pw,ess=−2εD0+H0.
```

tからu/8を引いた部分のw括弧は
8ε³F−4ε³Y+2ε³X²+z/4−b/4+ε²Xz−ε²bX+2ε²a。
−b/4+2ε²a=0、Yの逆から出るεXと−ε²bXが相殺し、
{u/8,w}=−z/4−εu/4と−εu/8を合わせると表示したPtwを得る。
これで全六括弧／全150斉次係数の恒等式が原式から確定する。

ε≠0はactual B1なので全元Jacobiを満たす。
上の全係数がε多項式であるから、全350 projective Jacobiもε=0まで恒等零である。
任意のKでこの議論が成立し、未知のcompatible高次項を捨てていない。

signed PfはPtu Pzw−Ptz Puwであり、全ε一次・二次係数が相殺して
正確にf=u+z²である。別方向にも
old Pf=−H、det=−8ε³a³、H=f/a²、8ε³a=1から同じ式を得る。
従ってactual五次D=h³(hu+z²)は全曲線で非零かつ一定。
ε=0は§1の全K正常形そのもので、rank四を失わない。

## 5. 全source帰属と零補集合

任意Kでε≠0の各点はB1の閉成分に入り、曲線全係数はε=0で正則。
よって§1の正常形の全六係数は同じB1閉成分に入る。
可逆normalizerを戻せば、元A10のδ≠0の全像が入る。
元A10は既約で、δ≠0は稠密開、元Pの全係数は十変数の多項式であり、
Pf=fは全域で非零。同じ閉成分の逆像は閉なのでδ=0の全β,γ,e,Kも含む。
δを分母とする正規形へ直接零を代入したのではない。

これによりconstant modularのdegf≤2で残ったsource IIIは全点B1へ帰属する。
field cの任意非零値、κ≠0、H/frameを戻す部分は既必要十分入口のactual
affine/PGL/非零tensor尺度に従い、B1のPGL不変な閉性で保たれる。
source IIのNEW25、source IのM11への別actual極限と合わせて全三源を処理出来るが、
このノートは他sourceの成分性を独立に再認証しない。
新成分増分零、smooth性・局所一意性・全軌道分類は非主張。
