# 非被約R・m三generic21・primitive degrees(1,0)の全first connection

2026-11-07、主義務G。producer純幾何。operator/G閉鎖はまだ未証明。
constant-root全sourceとの補集合を扱い、degree(0,1)/(0,2)、generic triple、m二は残す。

## 1. 全factorizationと接続

actual Dは被約既約quintic、rank二plane S沿い位数三、zero cubic fの支持は非被約。
G3のdouble normal rowがprimitive degree一、simple rowがdegree零とする。
UFDによりactual constant normal GL2で、nonzero amplitudeを自由に残し

`G3=z(Az+Cw)^2`

へ入る。A,Cは独立base linear forms。coefficient gcdは一。
normal root determinant Cのzero支持はRに含まれ、従ってCはRのline成分である。

zを保存する全conormal connectionはlower triangular β=[[p,0],[r,q]]。
moving double rowの固有connectionはq+dlogCなので、G3のparallel性は

`q=dlog(f/C²), r=−d(A/C)+(A/C)(q−p)`。

flatnessはp,q closedと同値である。rのtwisted closed式は上の式から直接従う。
pをlogarithmicとは仮定せず、global fβのregular homogeneous polynomial条件を課す。
flat connectionだけからglobalK存在は推論せず、後続batchで全三Kを独立に残して解く。

## 2. 全base/root配置

R=2L+Mではf=a²b、C=a又はb。
C=aに対するp0=(A=C=0)がL∩Mかそれ以外かでA=b又はc。
C=bでも同様にA=a又はc。
base stabilizerとnormal BorelでAのC倍を吸収し、nonzero diagonalで係数を戻せる。
この四配置は全A,Cを覆い、任意のAの一点を選んで全supportを推論しない。
R=3Lではf=a³、C=a、全A独立Cは同じbase affine orbitなのでA=bで尽くされる。

## 3. 二重line型のglobal closed p

全closed connectionは

`p=(dlogf)/3+αdlog(a/b)+γd(b/a)+δd(c/a)`。

a非零chartのclosed Laurent/polynomial formとa∞/b零のpole orderから全三方向を得る。
高位poleはexactで、primitiveはP1/a modulo constants。
C=aならq=dlogb。
f*rに出るa^{-1}daの係数はnonzero factorを除きA*b*(γb+δc)なので、global性はγ=δ=0を強制する。
両A=b,cで同じ結論。従ってこの二配置は全α∈Cのwhole linesである。

C=bならq=2dloga−dlogb。
f*rのb^{-1}db係数はnonzero factorを除きA*a²*(α−1/3)で、他のexact項はこれを打ち消せない。
従ってα=1/3が必要で、γ,δ全Cを残す。
この時rは

```
A=a: r=γdlog(a/b)−δdc/b+δc da/(ab),
A=c: r=−dc/b+c da/(ab)−γ(c/b)d(b/a)−δ(c/b)d(c/a).
```

どちらもf*rはglobal quadratic homogeneous one-formになる。
非log exact二方向を消していない。

## 4. 三重line型

C=a,A=b、q=dloga。全p=dloga+d(P2/a²)、P2 homogeneous quadratic modulo constants*a²。
f*rに残る高位項は2A*P2 da/aなので、gcd(A,a)=1からa|P2が必要。
従って全pは

`p=dloga+γd(b/a)+δd(c/a)`

で尽くされる。全γ,δを残す。

## 5. actual root-row stabilizerの全逆正常形

C=b,A=cのcanonical rowc z+bwについて、a'=ta,b'=rb,c'=sc+ubを考える。
normal z'=kz,w'=ℓw+mzでrℓ=sk、m=−uℓ/sとすれば、rowはnonzero scalar倍のc'z'+b'w'へ戻る。
constant simple zも保存する。pのexact parametersは

`δ'=tδ/s, γ'=t(γ−δu/s)/r`。

従ってδ非零ではδ'=1,γ'=0、δ零/γ非零ではγ'=1、両零は別源となる。
C=b,A=aではc'=sc+ub+vaを使い、normal diagonalsでa z+bwを保つと同じparameter作用を得る。
R=3L,C=a,A=bでもc shearとa/b scales、rowを保つnormal Borelで同じ三strataを得る。
全inverseはminus shearsとnonzero diagonal inverse。
πS amplitudeはnonzero tensor scaleで戻せる。

従って全first connectionsは

- R2L、C=a,A=b/c：二whole α lines。
- R2L、C=b,A=a/c：各δ1γ0、δ0γ1、δγ0の三fixed。
- R3L、C=a,A=b：同じ三fixed。

合計十一源で尽くされる。

## 6. raw source・finite operatorとdimension upper

βの各entryからhorizontal polynomial one-form
`ωij=3f βij−δij df` を作る。
ωij=−x_base×Vijならtracefree linear Vはcurl(ωij)/3。
raw mixed bracketは `{xk,n_i}=−Σj Vij,k*n_j`。
base Pfaffian liftはgradient fとする。
normal quadratic Kは初めから全k0z²+k1zw+k2w²を保持し、global actualD2零を全係数で解く。
P0全350 Jacobiも確認し、inconsistent全源は空、壁にだけ存在する源があればその全補集合を復活する。

残全源にはfullP1三十/fullP2九と、actualD3が全line `C*z(Az+Cw)^2` に入る必要条件を課す。
lineの非零amplitudeを固定せず、全completekernel像が零ならm三源を空とする。
最弱wholeprojective上界はB+k1+k2+5。

fixedplane source Bは、Rorderedlines四/πS amplitude一、constant simple root一、
moving row dataを足して得る。moving rowのdet指定はprojectiveP5に二条件なので三dim、
simple rootと合わせ四dim。intersection点がL∩Mならさらに一減する。
従ってC=a,A=b/cのwhole αはB=9/10、固定α壁は一減。
C=b,A=aの三strataはB=10/9/8、A=cなら11/10/9。
R3Lではline二＋πS一＋root data四＋exactparameter2/1/0でB=9/8/7。
parameter source upperを使い、特殊scalar first connectionに一般orbit dimensionを外挿しない。

一回のbounded二QQ(α)＋九QQ batchで全K/P0/30/9/actualD3lineと全parameter例外を扱う。
未達源は全supportを残し、Ore/actual divisorの構造へ移る。
本finite化だけをG完了とはしない。
