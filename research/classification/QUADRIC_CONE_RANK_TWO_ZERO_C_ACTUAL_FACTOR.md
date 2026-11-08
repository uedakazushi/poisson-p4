# Rank二quadric cone・C零：全higher源に二次因子

2026-11-08、PRIMARY E、純証明。新M1／operator／parameter profile／jet／CAS計算は無い。
保存済みfull33源の一つのJacobi blockだけで、C零の全sourceを閉じる。

**命題。** `S=(h=q=0)`、`q=a*c−b²`のquadric cone上で一般restriction rankが二、
実Dが被約既約quinticであるrank四Poisson構造を考える。
既受理のactual A正規化でconstant three-variable part Cが零なら、その全sourceは空。
Fv零／非零、全sigmaとその全共鳴、全Y／T／P2、全vertexを含む。
constant C非零のsemisimple／nilpotent二sourceは本命題の補集合として保持する。

## 1. 元全33源と正当なOre入口

全sourceと全normalizer逆は `QUADRIC_CONE_RANK_TWO_SEPARATE_REVIEW.md` §§2–5を使う。
strict chart h=1、座標(a,b,c,v)において、元Jacobiのhighest-symbol帰結により
constant Hamiltonian pencil Dは零となり、全sourceは

`Pi=B(a,b,c)+(Y(a,b,c)+v*T(a,b,c)) wedge partial_v`

である。degB≤2、degY≤3、degT≤1。
この形を入口で仮定した旧失敗結果を用いず、D=0を証明した後のwhole源を用いる。
元fullprojective Jacobiは特に `Jac(B)=0` を含む。

C=0、Fv=0なら全actualPfaffianのh⁰係数が零でH因子となり、既に空。
Fv非零なら、既存SO(q)とactual scalar／normal h-scalingの全逆により
Fv=b又はFv=cの二軌道へ戻る。q-dual norm非零／零の両補集合であり、
restrictionの零divisor Rの種類を切っていない。

## 2. 保存済み全sourceから読むbase Bの完全式

入力 `CLASSIFICATION_QUADRIC_CONE_RANK_TWO_PENCIL_20261107.json` の
cases[2]（type=zero）、field `pencil_B0_point_and_all_directions` を用いる。
この配列のindex0はparticular、index1,…,33はfamily係数p1,…,p33。
元の全base係数に寄与するのは正確に

`A=p1, C0=p2, D0=p8, r=p7, s=p14, t=p20`

とFv三係数p28,p31,p33だけである。他の全係数はY、Tに残し、零としない。
同reportのFv matrixと、`TRIANGULAR_BOUND` reportのC零二recordにより
Fv=bはp31=1,p28=p33=0、Fv=cはp33=1,p28=p31=0。
ここでC0,D0は自由constant、正規化のconstant bivector Cとは別の記号である。

従ってFv=bの**全B**は

```
B_ab=A−r*b+s*a,
B_ac=C0−r*c+t*a−2q,
B_bc=D0−s*c+t*b.
```

Fv=cの**全B**は

```
B_ab=A−r*b+s*a+2q,
B_ac=C0−r*c+t*a,
B_bc=D0−s*c+t*b.
```

全six coefficientsを保持したliteral式である。
sigma、残weight一／二／三係数はBに現れないため、以下はそれらを一切制限しない。

## 3. Fv=b：Jacobiの五係数から全baseを決める

三変数b-vectorを `bvec=(B_bc,−B_ac,B_ab)` と定める。
Jacobiは `bvec dot curl(bvec)=0` と同値である。
上の第一式を代入して、半分のscalarは正確に

```
−t*a*b−r*b*c+2s*b²−D0*a+A*c+s*C0−r*D0−t*A.
```

ab,bc,b²の全係数から t=r=s=0、続くa,cから D0=A=0。
constantの残式は自動で零。C0にunit条件は無い。
従って全Poisson源で

`B=(C0−2q) partial_a wedge partial_c`。

これはgeneric sigmaだけの議論でも、別のB normalizerを追加した議論でもない。
保存済み全元source上の必要Jacobiの直接帰結である。

## 4. Fv=c：全係数比較

第二式で `−bvec dot curl(bvec)/2` は正確に

```
r*c²+t*a*c−2s*b*c+2D0*b−C0*c+r*D0−s*C0+t*A.
```

c²,ac,bcからr=t=s=0、b,cからD0=C0=0。
Aを反転したり零としたりせず、全Poisson源で

`B=(A+2q) partial_a wedge partial_b`。

ここでも全sigmaと全higher係数をそのまま残す。

## 5. chart因子から元actual五次のglobal因子へ

全Ore式の4×4Pfaffianは、Fv=bでは非零constantと符号を除いて

`(C0−2q)*(Y_b+v*T_b)`、

Fv=cでは

`(A+2q)*(Y_c+v*T_c)`。

他のY／T項がこの積を壊すことは無い。全sourceの唯一のbase bracketを用いた
元actualPfaffian恒等式だからである。
第一因子のdegreeは正確に二、第二は高々三。
従ってdegree五の元Euler-Pfaffianをhomogenizeすると、第一因子はそれぞれ

`C0*h²−2q`、`A*h²+2q`

という非定数quadric因子になる。残cofactorが零ならactualD=0であり、やはり対象rank四でない。
どのC0,Aを取ってもqのcoefficientは非零なので、単なるchart上のunitへ退化しない。
原projective tensorへEuler gaugeを戻してもactualE wedge Pi²は変わらない。
全既約quintic枝が空であることが従う。

以上はC零の全元源と、Fv零／両非零軌道の全逆を合わせたpointwise命題である。
constant C非零二source、他の導手曲面又は可約actualDの親成分は本稿で解決したとはしない。
