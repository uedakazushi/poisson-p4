# 全marked三次曲線代数のglobal normal form：復元した部分還元

2026-11-07、G、producer plane_collision_review。
VM交代前の未pushノートの数学内容を会話引継ぎから復元した**別identity**でありbyte同一ではない。
旧expected SHA256は`d605e543dad68d54c2533f43ab35905dff0ecbf010a4baf3aafcbe6018497fb3`、
旧manifest expected SHA256は`d9939ed0d5acfa781517be829fc2d42629dfacf70246bdd2bfcbbb357dd8ef36`。
新CAS無し、独立review前の部分証明。X全体の分類又は全G上界をこのノートだけで主張しない。

## 1. 全global algebraとframe

前ノート `RANK_TWO_QUINTIC_CANONICAL_TRIPLE_ALGEBRA_RECOVERED.md` の記法を用いる。
V=P1上のunital algebraはA=OV⊕OV(-1)^2。
H0(A(1))/H0(OV(1))=C²である。
nonzero rank一mark aがbaseのlinear sectionであれば、generic multiplication rankは三となって矛盾。
従ってそのquotient classは非零。
constant GL2のframe changeとbase linear translationにより(1,a,z)をglobal module basisに出来る。
全変換はblock triangularでconstant determinant非零、逆も同じ形を持つ。
rank一行列の恒等式M²=(tr M)Mはdense openからglobalへ延長し、
`a²=(tr a)a`を与える。base pointを反転していない。

## 2. Node、二markのquotient classesが独立

(1,a,b)がbasisとなる。rank一multiplicationのimageはgenericにa又はbのlineであり、
飽和性から全global積は

`a²=ell*a, b²=m*b, ab=0`,  ell,m∈H0(OV(1))

となる。これらの式は全ell,mでassociative。
base GL2、a,bの独立nonzero尺度、aとbの交換により全型は

| 型 | (ell,m) |
|---|---|
| B1 | (v,w) |
| B2 | (v,v) |
| B3 | (v,0) |
| B4 | (0,0) |

で尽くされる。v,wはP1の独立linear forms。
B1は三直線chain、B2はnoncoplanar concurrent triple、B3はdouble line+line、B4はsquarezero tripleに対応する。
全special fiberを含み、ell又はmの零点を除いていない。

## 3. Node、二markのquotient classesが従属

b=c a+h、c∈C*、hはbase O(1) section。
ab=0とa²=(tr a)aからh=−c tr a。
tr aがgenericに非零なら、multiplication by bは
`c(M_a−(tr a)I)`でrank二となりrank一条件に反する。
従ってtr a=h=0、b=c a、a²=0である。
これはT上のrestrictionの等式である。
X上ではs2−c s1がspan(u0,u1,u2)の非零sectionであり得るので、二つのΓを同一視していない。

## 4. 一つのnilpotent markの全global型

これはcuspとnode従属型の両方を含む。
basis(1,a,z)、a²=0、multiplication aがrank一よりaz=ell a。
一般式z²=P2+Q1 a+R1 zのassociativityは
`P2=ell²−R1 ell`を与える。
z→z−ellで

`a²=az=0, z²=kz+qa`, k,q∈H0(OV(1))

へ変換される。残るframe change
`a'=alpha a, z'=beta z+gamma a`、alpha beta非零、gamma constantでは

`k'=beta k, q'=(beta² q−beta gamma k)/alpha`。

base GL2も用いると全型は

| 型 | (k,q) |
|---|---|
| N1 | (v,w) |
| N2 | (v,0) |
| N3 | (0,v) |
| N4 | (0,0) |

となる。独立、k非零でq比例、k零でq非零、両零という全補集合である。
q比例型はgammaで消せる。必要なpivotはconstantのみ。
N1のv=0 fiberはcurvilinear triple、N2/N3のv=0 fiberはsquarezero、N4は全fiber squarezero。

## 5. 正確なscope

nodeはB1–B4とN1–N4の八marked global algebra、cuspはN1–N4の四型に還元された。
全OV-algebra isomorphismと逆変換を与え、generic field型だけでglobal familyを捨てていない。
しかしこれはTのalgebraとmarkのnormal formであり、
canonical三fold Xへの全extension、normal方向の係数、ambient Poisson descentの全normal formではない。
そのextension parameterを零と仮定していない。
この段階の部分結果を全G≤24又はE被覆と読み替えない。
