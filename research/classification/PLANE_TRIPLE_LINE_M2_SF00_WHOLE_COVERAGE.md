# R=3L・m二SF00：元bad-vector対とcanonical Casimirによる全点被覆

2026-11-08、root純手証明、独立確認前。新CAS/operator/module/jet/GB無し。
actualD被約既約quintic、rank二restriction plane、f=a³、normal位数二、generic SF、
primitive根degrees(0,0)の全sourceを扱う。
全first-symbol入口は `PLANE_TRIPLE_LINE_M2_SQUAREFREE_CONSTANT_ROOT_SYMBOLS.md` と
同別作者reviewで受理済み。四sem固定、一translation固定、whole(t,d非零)nilが全源を必要被覆する。

## 1. 両constant rootの元全vector式

raw P0はP_f、mixed uX/vY、normal kappa*u*vで、X,Yはtracefree linear、Xa=Ya=0。
全P1はraw36/Euler六、全P2九を許す。
u-rootのbad v² u-rowをC、v-rootのbad u² v-rowをC'とする。
両方ともEuler六gaugeで変わらず、零ならそれぞれu/vがfullactualDのproper因子となる。
P2九はnormal columnを持たないので、inhomogeneous particularでもこの必要性は変わらない。
ゆえに全対象でC,C'とも非零が必要である。

元projective quad(i,j,u,v)のv³又はu³係数を取り、strict Jacobi自体を零とは置かずに

```
E∧(kappa*C+DY*C)−(2Y−X)∧C+P_(C(f))=0,       (U)
E∧(-kappa*C'+DX*C')−(2X−Y)∧C'+P_(C'(f))=0.   (V)
```

を得る。これは既独立SF02証明で使った全vector式とその二根交換である。
curlは

```
(2kappa*I+4DY−DX)C=0,
(-2kappa*I+4DX−DY)C'=0.
```

を与える。以下ではcurlだけでなく必ず元(U)/(V)にも戻る。

## 2. 四semisimple全源は空

全canonical modelsはX=H=b∂b−c∂c、Y=±H、
kappa=±3/2（同符号fields）又は±5/2（反対fields）。
kappa非零だからcurlのa-rowからCa零、従ってC(f)零。
(U)のab/ac成分はa*(kappa*I+DY)C零なので、DY*C=−kappa*C。
これをcurlへ代入するとDX*C=−2kappa*C。
しかしDXの固有値は0,1,-1、−2kappaは±3又は±5である。
従ってC零、既約actualDと矛盾する。
translationsやnonsemisimple境界を先に落としたのではなく、全first-symbol還元後の四全fiberの空性である。

## 3. whole二parameter nil源は全t,d非零で空

N=a∂b+b∂c、N²=a∂cとしてX=N、Y=tN+dN²、kappa零、t任意、d非零。
curl対は

```
((4t−1)N+4dN²)C=0,
((4−t)N−dN²)C'=0.
```

t≠1/4なら(4t−1)I+4dNは可逆でNと可換、従ってC∈kerN=C*∂c。
その場合DY*C=0、C(f)=0で、元(U)のbc成分は

`(2t−1)*a*C_c=0`。

C非零が可能なのはt=1/2だけ。
対称にt≠4ならC'∈kerN=C*∂c、元(V)から

`(2−t)*a*C'_c=0`。

C'非零が可能なのはt=2だけである。
tが1/4,4のどちらでもない領域では両必要値が一致せず矛盾。
t=1/4では第二のcurlが依然kernel∂c、(V)の2−t=7/4が非零なのでC'零。
t=4では第一のcurlがkernel∂c、(U)の2t−1=7が非零なのでC零。
従って**全t、全d非零のnil sourceが空**。
d一のslice、generic t、canonical translationのregularity又は高次局所環は用いていない。
先に得たcanonical affine/Ore/A分岐の必要補題は正しいが、この短い元係数障害だけで残余全体を排除できる。

## 4. 一つのtranslation源は全essential三

残るP0はX=a∂b、Y=a∂c、kappa零。
fullrawP/P2九でrow divergence Ru=Rv=0。
canonical divergencezero strict liftPi=P+E∧R/5はPi_uv=0、actualFのCasimir liftである。
K=C(u,v)上で

`V1=Pi#du/u, V2=Pi#dv/v`

はnormalsを固定し、全base上affine fieldsである。
Pi_uv=0とstrict Poissonから[V1,V2]=0。
linear partsは(4N1−N2)/5、(4N2−N1)/5、N1=a∂b、N2=a∂c。
全higher translationsをp,qとして

```
V=(4V1+V2)/3=N1+p,
W=(V1+4V2)/3=N2+q.
```

と出来る。p,qはKの任意rational homogeneous degree一vectorsでよく、polynomialとは仮定しない。
commutationはp_a∂c−q_a∂b零、従ってp_a=q_a=0。
actualFはV,Wに殺されるから

```
(a+p_b)F_b+p_c F_c=0,
q_b F_b+(a+q_c)F_c=0.
```

係数determinantはK[a]でleading a²係数一、非零。
従ってF_b=F_c=0。元Fはpolynomialなのでこれは全域の恒等式で、actualDはa,u,vだけに依存する。
非零既約quinticはessential二以下ではないため、本質三変数の既八23成分へ全点入る。
rational translationsをglobal PGLへ変換したという議論ではなく、元polynomialの二偏微分を直接零にしている。

## 5. 全source結論と残余

独立受理済みの全SF00 first-symbol coverと合わせ、四sem/nilwhole全空、一translation全essential三。
従って**R=3L・m二SF00の全点は既八23成分へ被覆される**。
独立読解前はproducer純証明の状態を保つ。
SF01五源の全空性、SF02/03/11/12既空、normaldouble既八23被覆を別に合わせれば全R3L/m二の全点被覆となる。
R=a²b、m三又は全G/E・上界は本稿だけでは閉じない。
