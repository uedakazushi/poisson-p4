# concurrent m三 SF root geometry：独立確認した有限入口と残余

2026-11-08、PRIMARY G、別作者 `/root/quartic_review`。
producer `PLANE_CONCURRENT_M3_SQUAREFREE_ROOT_GEOMETRY.md`
(`2d5fb18de3ed1d33b28612a3f3e66c9ab06457e35f57884519293e1bf345a688`) を独立読解・再導出した。
**記載された部分的 whole geometry を受理する。**
finite constructions の全 candidate 列挙・global β/K 検査・全 higher source G 上界は未完了。
新 CAS、operator/source/minor 計算、jet 又は GB は実行していない。

## 1. 全 covers と c-degree の二分

対象は reduced concurrent f=ab(a−b)、rank二 plane、actualD 被約既約 quintic、normal 位数三 SF。
U=(P1\{0,1,∞})×A1 の有限 étale root covers は marked base curve cover の pullback。
π1=F2 の permutation image として 1/C2/C3/S3 を全て保持する。
C2 は二 punctures、C3 は二 punctures 又は全三 punctures に branch。
S3 の三 inertia は parity と surjectivity より (2,2,3) のみである。
Riemann–Hurwitz は Galois curve の genus 零を与え、三 marked branch values を合わせれば
標準 D3 model r→ζr,r→1/r を用いてよい。

splitting curve の generic field K 上で root pairs を primitive K[c] pairs に取る。
全 c での SF parallel 性から pair determinants と content は K[c] units。
非constant polynomial なら generic base 上の collision curve が U に入り矛盾する。
degree sum は deg_c G3<=2。
positive degree が一つだけなら他二 independent constant pairs との determinants がその pair を constant にする。
従って全 roots c-independent 又は (0,1,1) だけ。
transitive C3/S3 では三 degrees が同じなので後者は無い。
後者の unique fixed root が global primitive degree 一以上なら、residual coefficient degree は一以下なのに
c-degree 二が必要となり矛盾する。従って実 constant normal root u へ移せる。
K 上の argument を global UFD/Pic 零の仮定に置換していない。

c-independent primitive connected incidence の class は (e,3)、e<=2。
root curve normalization から normal P1 への degree は e。
全 content は base divisor R-supported、degree 2−e。
curve の rational parametrization は root maps の分類に使い、ambient PGL 同値には用いない。
最後の normal postcomposition だけが実 constant normal GL2 である。

## 2. constant root と C2 moving root の有限 construction

constant root では G=u(Au²+2Buv+Cv²)、C と T=AC−B² は非零 binary forms、R-supported。
C の全六 degree二 divisors、T の全十五 degree四 divisorsを保持する。
ordinary triangular connection の記載式は、parallel tensor と determinant twist から得られる。
C が二 distinct lines なら B²≡−T mod C は二点の値を有限に固定し、B mod C が有限となる。
C=L² なら ord_L T=0 は value/first jet の有限 square roots、ord一は不可能、ord>=2 は
B=L(sM+tL)。t は actual constant unipotent、s は零/非零を分けて diagonal scale で扱う。
T/C の nonzero amplitudes は p/q の dlog に現れず、r の共通 amplitude は normal diagonal で除ける。
全 leading amplitude は同一 connection の parallel space に戻す。
これは connection candidates の有限化であり、未検査の simple-pole/global-K candidate を実在とは呼ばない。

C2 fixed root degree 一の residual incidence は (1,2)、moving root map は degree 一。
branch cover r²=a/b 上で moving roots を r,−r、fixed root を Möbius ρ(a/b) とできる。
ρ(r²)=r の homogeneous cubic divisor は {0,∞,1,−1} 上に支持され、全20 multisets で尽くされる。
係数 (C,−A,D,−B) が ρ を scalar 以外一意に復元する。
AD−BC!=0 と global source 条件を戻し、もう一方の collision は r→−r で保持する。

## 3. C3 の全 degree二 maps と elliptic branch の全空性

二 branch C3 は P1_r、r³=a/b。
degree 一 map は normal PGL2 後 r、degree二 map は quotient involution τ のみで決まる。
collision equation τ(r)=ζr は

    ζC r²−A(1+ζ)r−B=0.

全 degree二 divisor は {0,∞,1,ζ,ζ²} 上の15 multisets。
この coefficient map は injective であり、detτ!=0 を課して全候補を得る。
σ² との collision はこの divisor の deck translate で、別の無記録 family は無い。

全三 branch C3 は smooth elliptic cubic w³=ab(a−b)。ramification point を origin に取れば σ=[ζ]。
normal map は degree二、Pic² line bundle F の complete pencil で involution τ_s(x)=s−x。
τ_s(x)=σx の唯一解は (1+ζ)x=s、1+ζ=−ζ² は unit。
U 上 collision が無いことは、この解が三 fixed points の一つであることを強制する。
H=3O、H−F の degree一 point は p=−s なので、p も同じ fixed subgroup の一つ。
従って全 maps は三 ramification points からの projection pencils。
norm/content 計算は G=a²u³+b(a−b)v³ を与える。
diagonal ordinary connection の一 character が (1/3)dlogf なので、Q=b(a−b) は fullD3 の像 C(f) に入る必要がある。
C(f)=kQ を mod b と mod(a−b) で読むと C(b)=C(a−b)=0、従って C(a)=0、k=0。
target coefficient 非零に反し、**全 elliptic C3 branch は空**。
他の rational covers をこの空性へ含めていない。

## 4. S3 degree二 maps の有限 construction

root field q=r+r^-1、t=(q³−3q)/2、normal map ψ(q) の degree は一又は二。
degree 一の primitive cubic b(u³−3uv²)−2av³ と全 R-supported line content を保持する。
degree 二では ψ の quotient involution τ に対し τ(q(r))=q(ζr) の polynomial を独立に展開すると

    Cζr4+Aζ²r3−(C+B)r²+Aζr+Cζ².

全 degree四 divisor は {0,∞}∪μ6 上の330 multisetsで尽くされる。
C は r4、A は r3、B は r2 coefficient から復元でき、coefficient map は injective。
既 collision q(r)=q(ζr) も同じ branch-preimage set 内にある。
従って all candidates を含む finite construction として正しい。
330個が全て実 connection/実P0 である、又は既に算術列挙された、とは結論しない。

## 5. 二つの c-dependent profiles と全 nonlog 部分

C/T c-independent、deg_c B=1 より C|B_c²。
primitive binary degreeから C=L²、B_c は L 倍となり、global A の polynomial 性から B0 も L 倍。
実 c-shear/scale とその逆で B=Lc。
A=c²+T/L² の polynomial 性は T=L²M2 を強制する。
全 r=(c/(2L))dlogM2−dc/L を戻すと、M2 に L が入れば L 上の非零 double pole となる。
R reduced の simple-pole source に反するので、M2 は他二linesの平方又は積だけ。
二 profiles は whole necessary cover であり、off-diagonal を最初から log としていない。

記載 raw X/Y/Z は αX=dlogf−3p、αY=dlogf−3q、r=−αZ/3 と一致する。
square profile では det(E,X,Z)+Zf=6f、nonsplit では零。
従って D2=0 から K20=2/0、K11=K02=0 が正しい。
flat horizontal rank二 connection とこの一意 K は U 上 integrable bivector を与え、polynomial 性で全域へ延びる。

z=cu+Lv,w=u の frame で rational parallel monomials を列挙すると、
square は b²w³,bzw²,z²w の三dim、nonsplit は b(a−b)w³,z²w の二dim。
残 z³/b や square-root coefficients は polynomial/rational でなく除かれる。
全 leading amplitudes を一つの scalar 点へ固定していない。

first-source safe counts も確認した。
c-independent β は全 raw fields ∂c 方向、K=0 なので P0=∂c∧W0、全 c-shears は P0 を変えない。
R choices5+anchor1+normal conjugacy<=3 で B<=9。
c-dependent square の c-shift along b は connection を変えず、along L は normal unipotent と重複。
nonsplit では L 外の一 shift を追加し B<=10。
c-scale は normal diagonal と重複し、同じ変化を二重に数えていない。

## 6. canonical repeated02 補足と受理限界

§10 の追加必要条件も確認した。
canonical divergence-zero lift の barX/barY/barZ と元 bad C=λ∂c を使う strict coefficient の a,b projection では
全他項と Euler correction が消え、(DZ C)_(a,b)=0。
Z=3∇L×∇A2、L=a では Z_b=−3A_c、従って A_cc=0。
canonical化を global stabilizer の存在仮定と混同していない。
これは repeated02 の追加必要条件だけであり、全 repeated source の G 閉鎖ではない。

受理したのは全 cover/c-degree dichotomy、有限 algebraic candidate constructions、elliptic C3 全空性、
二 c-dependent full first profiles/parallel spaces、安全な first-source counts 及び上記必要条件である。
**残 finite candidates の実列挙、global β/K と full350、higher30/9 compatibility、actualD3 realizability、
concurrent 全G/Eは未完了。** known concurrent m三成分の存在とも両立する部分証明として保存する。
