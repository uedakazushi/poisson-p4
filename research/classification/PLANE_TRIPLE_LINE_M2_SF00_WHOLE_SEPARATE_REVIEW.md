# 三重直線・m二 SF00：全 bad-vector 対による独立全点被覆確認

2026-11-08、PRIMARY G/E、別作者 `/root/quartic_review`。
`PLANE_TRIPLE_LINE_M2_SF00_WHOLE_COVERAGE.md`
(`cf75dd21b793ba9ebc6f9b6842a0b42d133cda450bd048d647892922584cdef6`) を全読解し独立再導出した。
**全 SF00 源の pointwise essential-three 八23成分への被覆を受理する。**
新 CAS、operator/source/module、minor、jet 又は GB は実行していない。
長い canonical affine/Ore 必要還元は正しい背景として保持し、本結論には以下の弱い証明だけを用いる。

## 1. 全源の仮定と二つの proper-factor 条件

対象は rank 四 Poisson、actualD 被約既約 quintic、restriction rank 二 plane、f=a³、
normal 位数二 SF、両 primitive roots が constant の全源。
既独立確認済みの全 first-symbol cover は四 sem、一 translation、全(t,d)∈C×C* nil である。
全 raw36/Euler6 P1 と全9 P2を戻す。
constant root u の bad v² vector C、v の bad u² vector C' は Euler 六で不変。
C=0 なら全 u-row が u 倍、C'=0 なら全 v-row が v 倍なので actualF はその proper factor を持つ。
**P2 は全九の全てが normal column を持たない**ため、inhomogeneous particular でもこの結論は保たれる。
従って指定既約源には両 C,C'!=0 が必要。

projective weight 一の元全 vector 式は producer (U)/(V) と一致する。
二根交換は X/Y を交換して κ を −κ にするので、対称式の符号も正しい。
curl だけでなく元式の残 base coefficients を用いる。
f が非被約でも literal polynomial identity は同じであり、reduced/log 仮定を追加していない。

## 2. 四 sem 源の全空性

X=diag(0,1,−1)、Y=±X、κ=±3/2 又は±5/2。
κ!=0 と curl の a-row より Ca=0、従って C(f)=0。
元(U)の ab/ac 成分は a(κI+DY)C=0、これと curl から DX C=−2κ C。
DX eigenvalues は0,±1、−2κは±3又は±5なので C=0。
実際 DY の eigenvalues と κ の比較だけでも同じ矛盾を得る。
全 higher source の u factor が従い、四 sem 全源は空。

## 3. 全二parameter nil 源の全空性

X=N=a∂b+b∂c、Y=tN+dN²、κ=0、d!=0。
N²=a∂c、N³=0。
t!=1/4 では ((4t−1)I+4dN) が可逆なので第一 curl の kernel は kerN=C∂c。
C=γ∂c を元(U)へ戻すと DY C=C(f)=0、残係数は (2t−1)aγ=0。
γ!=0 の必要条件は t=1/2。
同様に t!=4 では第二 curl の kernel は C∂c、元(V)から t=2 が必要。
両条件は同時に満たされない。

全除外 complement も検査した。
t=1/4 では第二 curl はなお kernel∂c、元(V)の係数7/4!=0なので C'=0。
t=4 では第一 curl と元(U)の係数7!=0から C=0。
従って **全 t、全 d!=0** の nil 源は空。
d=1 slice、generic t、translation の regular 性、追加高次方程式は用いていない。
canonical affine source の A零/nonzero 必要分岐より前に、元 whole bad-vector 対だけで全源を排除できる。

## 4. translation 源の全 higher は essential 三

P0 の X=a∂b=N1、Y=a∂c=N2、κ=0。
full P の row divergence は Ru=Rv=0、Rbase=X+Y+normal-linear vector。
canonical divergence-zero strict lift Π=P+E∧R/5 は Πuv=0、actualF の strict Casimir lift である。
K=C(u,v) 上で V1=Π#du/u、V2=Π#dv/v は base affine fields、normals を固定し commute。
linear parts の可逆組合せから

    V=(4V1+V2)/3=N1+p,
    W=(V1+4V2)/3=N2+q.

p,q は任意 K-rational degree 一 translations のまま保持する。
commuting condition は p_a∂c−q_a∂b=0、従って p_a=q_a=0。
Casimir の二式は producer 記載通りで、その bc determinant は K[a] で leading a² coefficient 一。
非零 polynomial なので F_b=F_c=0 を得る。
元 F の polynomial 性により u=0/v=0 も含む全域の恒等式である。
有理 translation を ambient PGL とみなしたものではない。
従って F∈C[a,u,v]、既約 homogeneous quintic の essential degree は三であり、
既 accepted essential-three 全八成分 theorem により **全点** が被覆される。

## 5. 結論と他枝

全 first-symbol 必要被覆と合わせ、四 sem 全空、whole nil 全空、一 translation 全八成分。
よって **R=3L・m二 SF00 の全 pointwise E 被覆**が証明された。
新成分は追加しない。
SF01/02/03/11/12、constant/moving double 等は本 receipt 単独の範囲ではなく、それぞれの別証明と合わせる。
R=2L+M、m三以上、他 surfaces 及び全分類の完成は依然として別義務。
