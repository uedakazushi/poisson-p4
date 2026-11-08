# 三重直線・m二 constant double 全八 first symbols の独立確認

2026-11-08、PRIMARY G、別作者 `/root/quartic_review`。
root 純証明 `PLANE_TRIPLE_LINE_M2_CONSTANT_DOUBLE_EIGHT_SYMBOLS.md` を独立に全読解・再導出した。
**全源の二 nil 固定＋六 semisimple 固定 P0 への必要被覆を受理する。**
八つの higher sources の存在、G 次元、E membership は未結論。
新 CAS、operator/minor、source profile、jet 又は GB は実行していない。

## 1. 全 irregular pool と full bad-vector 式

対象は actualD 被約既約、rank 四 Poisson、rank 二 plane、f=a³、normal 位数二、
double normal quadratic の primitive root degree 零の全枝。
constant normal 右 GL2 後 G2=γa³u²、γ!=0。
parallel leading relation は二 ordinary diagonal characters を同じ θ にし、off entry も closed。
全 closedness と tracefree field 式から Xa=Za=0。
各 X/Z は全五次元 divergence-zero affine Hamiltonian pool を保つ。
irregular exact potentials を logarithmic 部分へ限定していない。

a-row 零なので E∧Jac=0 は strict Jac=0 を含意し、全残 Jacobi は [X,Z]=0。
affine area form の Cartan identity から det(Xbc,Zbc)=d0a²、actual γ=3k−d0。
det(E,X,Z) をこの段階で零とはしない。

全 raw36 P1/Euler6 と P2 九の bad u-row v²C に対する式は

    E∧(DX C)−X∧C+P_(C(f))=0.

curl は 3DX C=0。
元式へ戻すと ab/ac coefficients が C_a X_b=C_a X_c=0 を強制するので、X!=0 なら C_a=0。
従って X∧C=0、X の全 image は一定 line C*C に含まれる。
C!=0 と DX C=0 より rank DX=1、(DX)²=0。
curl kernel だけから field normalform を結論したのではない。
X=0 では同式は C_a=0 を課すだけで、X!=0 の結論を外挿しない。

## 2. k-unit・nilpotent Z の全源排除

既 independently accepted equal-diagonal proof cc1c を次の正確な仮定で使える。
k!=0、X=Y、Z nilpotent、actual leading γa³u²、γ!=0。
Ru=0、Rv=−2ku と canonical strict lift より ratio Hamiltonian は

    Mnormals=(3k/5)(u,v),
    Mbase=(Z−2kI/5)base+B(u,v).

kI−Z が可逆なので全 normal-degree-one translation が消え、
base degree5−m の operator は Z|Sym^(5−m)+k(m−2)I。
nilpotence より m!=2 が消える。
translation は normal degree 一なので original leading は変わらず、全 F=γu²(a+Ta)³。
fraction field C(u,v) 上の positive-base-degree factors から Gauss lemma で既約性矛盾。
ここで k!=0 は明示仮定であり、γ!=0 から k!=0 と推測していない。
全 d0 と全 raw36/9 higher を保持する。

## 3. X!=0 の全 commuting pair と actual inverse

Xa=0、rank-one square-zero X を a-line 保存 base GL3 で X=(ra+sb)∂c に戻す。
一般 Za=0、trace Z=0 の式へ [X,Z]=0 を代入する。
s!=0 なら Z_b=0、Z=(Cb+qa)∂c。従って d0=0、γ=3k!=0、nilZ で前節の全空性へ入る。
s=0 なら r!=0、全 Z=pa∂b+(Cb+qa)∂c。
k!=0 は同じ全空性。k=0 では d0=−rp、γ=rp!=0 なので p!=0。
c'=c+ell*b が q を q+ell*p へ動かすため、全 q を消せる。

old b=βb、old c=χc、tensor σ=βχ、relative normal scale e とすると

    r'=βr,    p'=eχp,    C'=eβ²C.

β=1/r を取り、C!=0 では e=1/(β²C)、χ=β²C/p、C=0 では e=1、χ=1/p。
全 units と逆が存在し、C=0 を別閉補集合として保持する。
従って X=Y=a∂c、Z=a∂b+εb∂c、k=0、ε=0/1、γ=1 の全二固定模型を得る。

## 4. X=0 の全 semisimple spectral walls

γ=3k!=0。
Z の bc block が eigenvalue 零のみなら、a-column を含む全三次元 Z も nilpotent で前節により空。
非zero eigenvalues λ,−λ の場合は実 base conjugation と a-column shear により
Z=λ(b∂b−c∂c) に戻る。relative normal scale は k と Z を同時に変えるため λ/k が残り、k=1 にできる。

translation failure は λ/k=±1。
それ以外では全 normal-degree m の operator は Sym^(5−m)(Z)+k(m−2)I。
m=3 では base spectrum 0,±λ,±2λ より λ/k=±1,±1/2、m=4 では ±2、m=5 は無し。
translation failure もこの全六壁に入る。
壁外では m=3/4/5 が零となり、全 F=γu²(a+Ta)³ で矛盾。
m=2 の Sym³(Z) に他の zero-weight monomials があっても、元の m二 leading は指定 γa³u² なので
未知の新 summand を加えられない。degree 一 translation はこの leading を変えない。

従って X=Y=0、k=1、λ∈{±1/2,±1,±2} の全六固定源が必要残余である。
λ=0、k=0、nilpotent と各 normalizer 境界は上の明示別枝で扱った。
全二 nil＋六 semisimple は必要被覆であり、その full higher30/9 は未検査。
他 root degrees、SF、f=a²b、m三以上及び全分類を本八表から結論しない。
