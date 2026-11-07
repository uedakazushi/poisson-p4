# m2 split normal roots：root-line class obstructionの別作者読解

2026-11-07、独立 reviewer `/root/quartic_review`。**PASS（純証明）**。
`PLANE_M2_SPLIT_ROOT_MATRIX_ENTRANCE.md` 及び
`PLANE_M2_ROOT_LINE_CLASS_OBSTRUCTION.md` の全量化・局所極・global bundle条件を読解した。
新しい線形代数、Jacobi、minor、GB、jet計算は無い。

対象は rank2 restriction plane S=P2、非零 **reduced** zero cubic R=(f=0)、
被約既約actual quintic D、normal位数m2、generic squarefree normal quadratic G2。
nonreducedR、double normal quadratic、m3への拡張はしない。
既済の Hamiltonian/associated-graded parallel補題、conormal connectionの全simple-pole源を入力とする。

## 全 root-matrix入口

twisted Sym2 の discriminant は degree零 invariantなので、parallel性から
`Δ(G2)=c f²`、c≠0。K=C(S)で二つのrootがsplitする。
UFD C[x,y,z] の Gauss補題を base/normal二gradingで適用し、全 common因子をhへ戻すと

```
G2=h ℓ m,    (ℓ,m)^t=M(u,v)^t,    h detM=f
```

を得る。各rowはprimitive homogeneous polynomial vector。
scalar cの平方根は複素非零constantに吸収できる。
normal二rootの相異性からdetM≠0。
row degree d1≤d2、deg h=3−d1−d2より全flagsは
(0,0),(0,1),(0,2),(0,3),(1,1),(1,2)だけ。
root rowのcodimension二basepointを除去していない。

## θの極を先に評価する必要

rootがsplitして相異なるので、接続の infinitesimal parallel transportは各root lineを保存する。
有限交換は infinitesimal connection に追加の off-diagonal を与えない。
したがってrational one-forms θiが存在し

```
dM+Mβ=diag(θ1,θ2)M.
```

β(E)=Id、homogeneous row degree diから `ιE θi=di+1`。
flatness `dβ−β∧β=0` をこの actual rational gaugeで変換すれば dθi=0。
これは単にtrace又はprojectiveflat性だけからdiagonal flat性を推測したものではない。

各 irreducible divisor T の generic DVRで、primitive rowは少なくとも一成分がunit。
対応する scalar row equationからθiを解けば、TがR外ならregular、R上なら高々simple pole。
**M inverseを代入しただけではこの結論は出ない**。
両rowが点で同時に零になるcodimension二 locusも、このdivisorial boundを破らない。

closed rational1formのsimple-pole局所表示
`(a dt+Σ b_j ds_j)/t+regular` では、dθのt^-2項がb_j|T=0を強制する。
残るresidue a|T=0はclosednessからconstant。
各componentのconstant `ri dlog fi` を引くと全divisorial poleが消える。
ambient affine coneはsmooth/normalなのでcodimension二に新しいrational coefficient poleは無く、残りはpolynomial。
全formはEuler weight零（係数degree−1）なので残るpolynomial1formは零。
従って θiは各R-componentに沿うconstant-residue log形式の和である。
β全体がlogであるとは結論しない。

## detM上の二重極の全係数

f reduced、h detM=fなので、任意のp|detMはdetMを一回だけ割り、hを割らない。
M modpはrank1。primitive条件で二rowはいずれもgeneric点で非零である。
residueをr1,r2とし、

```
β=M^-1 diag(θ1,θ2)M−M^-1dM
```

と置く。後半は高々simple pole。
前半のpotential p^-2係数は `adj(M)diag(r1,r2)M dp`。
二rowをv,λv、λ≠0と書くと、scalar部分はadj(M)M=0、
残る係数は(r1−r2)の**非零rank1 multiple**である。
よってβのsimple-pole条件は必ず r1=r2 を要求する。
片方のrowがgenericに零となる偽補集合はprimitive条件が排除しており、rank1行列のargumentを無条件に使っていない。

従って η=θ1−θ2 はhの成分だけをpoleに持つclosed log形式。

```
η=Σ_(fi|h) si dlogfi,      Σ si degfi=d1−d2.
```

hがconstantならη=0だがEuler contractionはd1−d2。
**flags(0,3)及び(1,2)は全connection/全Poisson源が空**となる。
chosen matrix、generic point族のrank、full30/9のcutには依存しない。
primitive root rowsの全basepointとR交点のcodimension二境界も保持する。

残る(0,1)ではdegree二h上のresidue和−1、(0,2)ではη=−2dlogH、
(1,1)ではη=0を得るだけである。これを全source dimension/E membershipに昇格しない。
nonreducedRではdetMのsimple零及びhとのcoprimalityを失うので、今回の空性証明の対象外。

## triangle (0,0)のactual曲線

入口の別結論も正当である。
f=xyz、constant independent rootsなら G2=κxyzuv、κ≠0。
全global normal-weight源 Π0+Π1+Π2 を normalの実GL5 rescalingで
`Πt=Π0+tΠ1+t²Π2` へ送る。
全coefficientがpolynomial、t≠0では元とPGL同値、t0のactualDは正確にκxyzuvでSNC log-symplectic。
このirreducible曲線全体を含む周囲componentを取れば、それはSNC点に会う。
既受理MPS Theorem1.4の**SNC locusに会う全component** scopeから全曲線、従って元源はMPS40へ入る。
一点が既知閉和に入るだけという逆境界推論やT01局所一意性は用いない。
conic+line/concurrent/nonreducedfのendpointをSNCと誤認しない。

本receiptは新component計数、残flags、normal double/m3、全rank2planeの完了を主張しない。
入力全文hashと独立読解段階を対応report/jobに固定する。
