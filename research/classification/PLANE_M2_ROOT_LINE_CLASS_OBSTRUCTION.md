# m2 split normal roots：二つのdegree flagsは全空

2026-11-07、主義務G。新しい数学的証明。計算jobなし。
rank二restriction plane S=P2、zero cubic R=(f=0)はreducedかつ非零、
実Dは既約被約quintic、S沿い位数二、leading normal quadraticはgeneric squarefreeとする。
本ノートは root-matrix degree flags (0,3),(1,2) の全接続・全Poisson sourceが空であることを示す。
nonreducedR、normal double quadratic、m3、他支持へは単極の議論を流用しない。

## 1. 全root-matrixとnormal connection

`PLANE_M2_SPLIT_ROOT_MATRIX_ENTRANCE.md` のUFD証明により

\[
 G_2=h\ell m,\quad
 \binom\ell m=M\binom uv,\quad h\det M=f.
\]

各root rowはprimitive homogeneous polynomial vectorで、base degreesをd1≤d2とする。
f reducedなのでhとdetMは互いに素なR-component productsである。
conormal Hamiltonian connectionをnormal covector matrix βとして表す。
β(E)=I2、β−(dlogf/3)Iはhorizontal、各entryはR沿い高々単極である。
これはglobal H0Ω1(3)の完全源であり、log性を仮定していない。
flatness conventionはdβ−β∧β=0。

二rootはC(S)上でsplitして相異なるので、connectionは各root lineを保存する。
したがってrational one-forms θ1,θ2が存在して

\[
 dM+M\beta=\Theta M,\quad
 \Theta=\operatorname{diag}(\theta_1,\theta_2),\quad
 \beta=M^{-1}\Theta M-M^{-1}dM. \tag{1}
\]

flatnessからdθ1=dθ2=0、Euler contractionからθi(E)=di+1となる。

## 2. θiの単極はprimitive rowから証明する

あるirreducible divisor T のgeneric pointを取る。
primitive root rowの二成分はTで同時に割れないので、少なくとも一成分はlocal DVRでunitである。
そのunit成分を用いて、対応するrow equation

\[
 d(\text{row}_i)+\text{row}_i\beta=\theta_i\,\text{row}_i
\]

からθiを解く。
βがR外でregular、R上高々単極なので、θiも同じpole boundを持つ。
**MβM−1の形式的代入だけでこの単極性を推論しない。**
M−1自体がdetMで極を持つため、その推論では二重極を排除出来ない。

closed rational one-formに単極しかない場合、smooth generic pointで
θ=(a dt+Σb_j ds_j)/t+regularと書くと、dθの二重極からb_j|t=0となる。
残るresidue a|t=0はclosednessからconstant。
各componentでそのconstant logarithmic formを引けばdivisorial poleは全て消える。
codimension二にはrational one-formの新しい極が無く、同次degree−1のpolynomial one-formも無いので、
θiはR各成分に沿うconstant-residue logarithmic one-formsの和である。

## 3. detM沿いのresidue差は零

pをdetMのirreducible factorとする。
f reduced、h detM=fなのでpはdetMを一回だけ割り、hを割らない。
pのgeneric pointではMはrank一、primitive条件により両rowともnonzeroである。
θ1,θ2のp-residuesをr1,r2とする。
式(1)の第一項のpotential二重極は

\[
 \frac{\operatorname{adj}(M)\operatorname{diag}(r_1,r_2)M\,dp}
      {p\det M}.
\]

第二項M−1dMは高々単極であり、この二重極を消せない。
rank-oneかつ両row非零の2×2matrixに対し
adj(M)diag(r1,r2)M=0 mod p はr1=r2と同値。
実際、二rowをv,λv、λ≠0と書けば、scalar-diagonal部分はadj(M)M=0、
残る部分はr1−r2のnonzero rank-one multipleとなる。
したがってβの単極性は必ずr1=r2を要求する。

ゆえにη=θ1−θ2は**hの成分にしか極を持たない**closed logarithmic formであり、

\[
 \eta=\sum_{f_i\mid h}s_i\,d\log f_i,
 \qquad\sum_{f_i\mid h}s_i\deg f_i=d_1-d_2. \tag{2}
\]

ここでnonlogなβを捨てていない。
nonlog部分M−1dMを全て保持した上で、対角root-line connectionだけのclosednessを使った。

## 4. 全空の二flagsと残る三flags

h constantなら式(2)はη=0を強制し、一方Euler contractionはd1−d2。
degree flags(0,3),(1,2)はd1≠d2なので矛盾する。
従って **この二wholeflagsは全Poisson sourceが空**。
chosen matrix、generic minors、全30/9のdimension計算は不要。
root rowsのbasepoint/codimension二境界を含めた議論である。

残るflagsでは式(2)により次を要求する。

| flag | h | 全η |
|---|---|---|
| (0,1) | degree二R-component product | Σsi dlogfi、Σsi degfi=−1 |
| (0,2) | 一つのline H | −2dlogH |
| (1,1) | 一つのline H | 0、従ってθ1=θ2 |

さらに(0,0)はG2=κfuvで、triangleの場合は既actualSNC曲線により全source MPS。
conic+line/concurrentの(0,0)、上の残三flags、normal double quadratic、
nonreducedRとm3は未完了。
ηのfinite residue reductionだけをG上界又はE帰属と呼ばない。
本ノートの全数学はroot手読解を経たが、別作者独立receiptは別に記録する。
