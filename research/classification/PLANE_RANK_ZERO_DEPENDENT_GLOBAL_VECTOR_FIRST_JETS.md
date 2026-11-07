# rank零plane・dependent trace/common-fieldの全Jordan source

2026-11-07、主義務G。producer手還元。runner/current jobsのcheckpoint freeze中に保存する。
本noteはactual whole normalizations、base dimensions及び有限operator課題を定めるもので、
まだ全higher sourceのG非優越を証明しない。

## 1. 正確な入口

被約actual Dのrank零planeには新modular補題から全multiplicityで35/30/9 sourceを使える。
traceless co-Higgs symbolはglobal Vとconstant nonzero sl2 matrixの積とし、V≠0。
今回の全条件はtrace fieldが`T=tV`、t∈Cである。
pointwise proportionalでもlinear independentなcommon-image二planeは別noteで閉じた。
scalar/traceless零sourceとは混同しない。

normalssでは全weight零Jacobiが

`Π0=V∧((t+1)z∂z+(t−1)w∂w)+kzw∂z∧∂w`

を与える。
normalnilでは[T,V]=−bVとdependent性からb=0なので

`Π0=V∧(t Euler_normal+z∂w)+az²∂z∧∂w`。

いずれもimageはspan(V,∂z,∂w)に入りΠ0²=0。
actual rank四を持てばplaneに沿うDの位数は三以上。
既約quinticの位数五はDがbinary normal quinticとなりC上で因数分解するので、本Gでは三又は四である。
Π0がrank二だから全higherもrank二とは仮定しない。

以下Eijは係数行列、fieldはxj∂xi。
trace-free liftを使い、baseのconjugation、normal GL2、非零tensor scalarのみを使用する。

## 2. Vの全五Jordan型

nonzero traceless三×三matrixは次の全五型で尽くされる。

| V型 | whole representative | 全parameter | affine V-locus dimension |
|---|---|---|---:|
| regular semisimple | Rλ=diag(0,1,λ)_tf | λ∈C、λ≠0,1 | 8 |
| semisimple double eigenvalue | H=diag(1,1,−2) | none | 5 |
| mixed Jordan2＋distinct eigenvalue | H+N、N=E01 | none | 7 |
| regular nilpotent | J=E01+E12 | none | 6 |
| rank一nilpotent | N=E01 | none | 4 |

最初はnonzero eigenvalue differenceでscaleし、三座標のpermutationで任意のordered pairを選ぶ。
λ零/一はdouble-eigenvalue型へ別に保持する。
mixed型はnonzero double eigenvalueでscaleし、Hを保つdiagonal base changeでnonzero N係数を1へ直す。
nil型は全Jordan chain basisを使う。
各逆は使用した非零scalarとGL3の逆であり、arbitrary eigenvalue ratioを標本値へ固定していない。
V-locus dimensionsは順にconjugacy orbit六＋eigenvalue二、orbit四＋scale一、orbit六＋scale一、nil orbit六、nil orbit四。

## 3. normalssの全七parameter families

normalssではtensor scalar σがT,V,kを同時にscaleするのでtは不変。
semisimple/mixed型では§2のV-normalizationにそのscaleを使い、kは全Cのまま保持する。
nil Vではactual base gradingがVを任意のλ倍へ送る。
Jならdiag(λ²,λ,1)、Nならdiag(λ,1,1)が実現する。
σとλ=1/σを組み合わせてV及びT=tVを固定しkだけをscale出来るため、k零/非零は0/1の二whole casesへ直る。

| V | T | K | 全parameters | affine symbol base upper B |
|---|---|---|---|---:|
| Rλ | tV | kzw | λ≠0,1、t,k全C | 12 |
| H | tV | kzw | t,k全C | 9 |
| H+N | tV | kzw | t,k全C | 11 |
| J | tV | zw | t全C | 10 |
| J | tV | 0 | t全C | 9 |
| N | tV | zw | t全C | 8 |
| N | tV | 0 | t全C | 7 |

BはV-locus dimensionにt一、normal semisimple line二、K一を加え、K零では一つ下げた値。
全tを含み、特にt=0及びnormal factorが退化するt=±1も削らない。

## 4. normalnilではt零/非零を分ける

normal diagonal changeでz∂wとaを同じλnでscale出来る。
tensor scalar σを合わせるとT→σT、V→ηV、a→ηa、η=σλnであり、σ,ηは独立。
t≠0ならσ=η/t、λn=tでTとVを等しく出来る。
t=0は別に保持する。したがってT=εV、ε=0又は1の全二caseとなる。
元の任意の非零tへは上の可逆normal diagonalとtensorの逆で戻る。

semisimple/mixed Vでは§2のnormalization後もaは全Cを保持する。
nil Vではbase gradingを加えηを自由に動かしVを戻せるため、a零/非零をδ=0/1へ直せる。
T=εVも同時に保て、すべてactual変換である。

| V | T | K | 全parameters | affine symbol base upper B |
|---|---|---|---|---:|
| Rλ | V | az² | λ≠0,1、a全C | 11 |
| Rλ | 0 | az² | λ≠0,1、a全C | 10 |
| H | V | az² | a全C | 8 |
| H | 0 | az² | a全C | 7 |
| H+N | V | az² | a全C | 10 |
| H+N | 0 | az² | a全C | 9 |
| J | V | z² | fixed whole orbit | 9 |
| J | V | 0 | fixed whole orbit | 8 |
| J | 0 | z² | fixed whole orbit | 8 |
| J | 0 | 0 | fixed whole orbit | 7 |
| N | V | z² | fixed whole orbit | 7 |
| N | V | 0 | fixed whole orbit | 6 |
| N | 0 | z² | fixed whole orbit | 6 |
| N | 0 | 0 | fixed whole orbit | 5 |

一般t非零のBはV-locus dimensionにt一、normalnil line一、K一を加えた値。
T零では一つ、K零ではさらに一つ下げる。
normalizationでtを1へ直したことを、sourceのt parameterを捨てたことと混同しない。
失ったtはactual group orbitの方向としてBに含まれている。

## 5. 有限whole-source課題と未解決部分

合わせて二十一families：固定八、一parameter八、二parameter四、三parameter一。
full higher sourceは常にΠ1全三十、Π2全九。
fixed symbolのfull linear kernel upperの和をdとすると、ambient image upperはB+d+5。
各parameterのspecializationではactual orbit upper又は直接source countでbase dropを証明し、全minor補集合を戻す必要がある。

最初の有限batchは固定八＋一parameter八を全sourceで扱える。
generic six-shear kernelとD2injectivityが成立すれば、actual有限shearで全higherを消してrank二以下となる。
nonclosing sourceではcomplete original Jacobi/Pfaffianのfactor geometryを優先し、無限wall又は新jetのfarmへ進めない。
本noteのnormal formsだけでG非優越を宣言しない。

traceless零scalar源、Π0零のm四、common-image全四、normalnil b非零のaffine型、
以前閉じたm二全sourceは別の明示scopeとして保持する。
