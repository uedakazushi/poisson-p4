# unimodular二次pencil：primitive regular源の全延長はT16

2026-11-08、PRIMARY `E_multiple_hyperplane`。producer純証明、独立/root査読前。
Gの既非優越を変更せず、全C6/L24を含むE membershipだけを新たに証明する。新算術無し。

## 1. 正確な範囲と既全kernel

閉じた二形式α=J+β+γ、Pf(α)は非零constant、J非退化、
γ=df∧dg、f/gは独立な二次斉次式、そのpencilに非退化memberを持つとする。
まずgcd(γij)=1の部分を扱う。gcd非constantとγ零は削除せず、末尾の別全sourceへ戻す。

`UNIMODULAR_REGULAR_QUADRATIC_PENCIL_SUPPORT.md` §§1–4の全Segre分類と
独立 `CLASSIFICATION_UNIMODULAR_QUADRATIC_PENCIL_REVIEW_20261107.json` を入力とする。
四相異固有値は全低次延長が空。
残る十型では全24一次係数のclosed/wedge kernelが正確に八次元で、全βは

```
β=df∧v+u∧dg,  u/v constant one-forms.
α=Ω+dF∧dG, F=f+u.x, G=g+v.x, Ω=J−u∧v.             (1)
```

constant Pf条件はΩ∧df∧dg=0、Ω∧(df∧v+u∧dg)=0。
Ωが零又はrank二なら、非退化fのU=kerΩ上の二つの独立偏微分から
g−tfとv−tuがUを消し、Ω∧u∧v=0、Pf J=0となる既補題が全点を排除する。
従って全解でΩも非退化。
残八型はΩの全kernel上でPfΩ恒等零なので空であり、残すのはSS二blockとnil[2,2]の二型だけである。
低次kernelをgeneric pointから外挿したのではなく、既に保存された全kernelと全Ω sourceを用いている。

## 2. 全linear項は同時actual translationで消える

必要なら関数対(F,G)をdet一のGL2で取り直し、Hess fを非退化とする。
これは同じdF∧dGを表す操作で、未知のtensorを変えない。
(1)のquadratic/linear wedge条件は、Ωに関するPoisson括弧{F,G}Ωがconstantであることを意味する。
唯一のx0=−(Hess f)^{-1}uへactual affine平行移動してFのlinear項を消す。
移動後、{f,g+v'.x}Ωのlinear項が零なので

```
(Ω^{-1} Hess f).x(v'.x)=0.
```

ΩもHess fも可逆だからv'=0。従って両linear項が同時に消え、constant括弧も零である。
全u/vはこの一つのactual translationで保持される。
反転するのは存在を仮定したnondegenerate pencil memberとΩだけで、
Hess f零detのpencilをこのchartへ強制しない。
以後α=Ω+df∧dg、{f,g}Ω=0としてよい。

## 3. SS二blockはproduct source

旧全Segre表のSS残型は

```
f=(x1²+x2²)/2, g=(x3²+x4²)/2,
Ω=a dx1∧dx2+b dx3∧dx4, a*b≠0.
```

各plane内の可逆線形変換でf/gをxy/zwの非零倍に出来る。
従って逆bivectorは、constant二plane tensor J1/J2と
diagonal fields D1=x∂x−y∂y、D2=z∂z−w∂wにより

```
P=A J1+B J2+k D1∧D2
```

の形である。例えばΩ=a dx∧dy+b dz∧dw、γ=λd(xy)∧d(zw)で
matrix inverse規約なら(A,B,k)=(−1/a,−1/b,λ/(ab))。
固定volumeによるdualとの違いは非零constant scaleだけであり、membershipを変えない。

`UNIMODULAR_LINEAR_HAMILTONIAN_WHOLE_MEMBERSHIP.md` §3の全五period/edge bindingにより、
これは真completed T16のp=q=0部分そのもの。したがって全SS残型はT16に入る。
HQQ全pointwise theoremの完成を前提としない。

## 4. nil[2,2]を全sp4 polynomial族でCartanへ送る

旧全Ω/β kernelで

```
f=x1*x2+x3*x4, g=(x2²+x4²)/2,
Ω=t(dx1∧dx4+dx2∧dx3)+s dx2∧dx4, t≠0
```

が全源である。q=(x1,x3)、p=(x2,x4)、J2=[[0,1],[-1,0]]、K=−J2/tと置く。
bare skew matrixの逆を使いA=Ω^{-1} Hess f、B=Ω^{-1} Hess gと定義すると、直接block演算で

```
A=[[K,−(s/t)K],[0,K]], B=[[0,K],[0,0]], K²=−I/t²,
A³+A/t²=(2s/t³)B.                                (2)
```

s非零ならp(z)=t³(z³+z/t²)/(2s)としてB=p(A)。
今度はAを一点に固定せず、**全T∈sp(Ω)**上に

```
F_T=x^T ΩT x/2, G_T=x^T Ωp(T)x/2,
α_T=Ω+dF_T∧dG_T                               (3)
```

を定義する。ΩTとΩp(T)は対称（pはodd）で、Tとp(T)は可換。
従って{F_T,G_T}Ω=0、α_T²=Ω²≠0、dα_T=0が全Tで成立する。
四次元のskew inverseはPfで割ったcomplementary成分なので、
P_T=α_T^{-1}は元変数の次数二以下で、Tに関して多項式のactual Poisson familyである。
最高部だけを動かして未知の低次fiberを同一視していない。

sp(Ω)はaffine十次元の既約空間で、そのregular semisimple開は非空稠密。
その開ではsymplectic eigenframeでTの固有値を±a,±b、ab(a²−b²)≠0と取れる。
F_T/G_Tは二つのproduct invariant xy/zwの線形結合となり、
dF_T∧dG_Tのcross係数は非零constant倍ab(b²−a²)。
従って全generic Tは§3のT16へ入り、閉性からT=AもT16に入る。
これはRichardson定理又はcommuting varietyの既約性を外部入力として仮定する議論ではなく、
明示odd-polynomial family (3)の既約parameter空間とdense semisimple開を使う。

s零の全補集合は、固定t≠0、f/gを保ってΩ_sを全sで動かす実族へ戻す。
Ω_s+df∧dgは全sでPf=t²≠0かつ閉で、その逆もquadratic polynomial族。
全s非零がT16に入るので、s零も同じ閉T16に入る。
§2の元全u/vは固定actual translationの逆で全て復帰可能。

## 5. 全境界の正しいrouting

以上で、primitive gcd一・nondegenerate memberを持つ二次pencil最高部の
**全unimodular低次延長**はT16又は空となる。

* 四相異固有値と全八空Segre型は旧全30係数空性を使う。
* γ零は全零最高部の別証明へ戻す。
* rank一memberの平方等によりγに共通factorが生じる場合はfactor最高部の全入口へ戻す。
* 全member退化のpencilはこの証明のchart外。共通kernel無しは既全源空、
  共通kernel有りはfull三変数pullback延長の義務へ保持する。

gcd境界やsingular pencilをgeneric Qの平滑化だけで覆ったとは主張しない。
新成分零、全scheme/nonreduced厚みの主張無し。全A=c=0/Eの完成はまだ別課題である。
