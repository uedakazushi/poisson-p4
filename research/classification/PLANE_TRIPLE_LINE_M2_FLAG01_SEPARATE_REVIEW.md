# R=3L・m二SF01：五つの固定first symbolへの独立全源還元

2026-11-08、PRIMARY G、`/root/log_manuscript`。
`PLANE_TRIPLE_LINE_M2_SQUAREFREE_FLAG01_FIVE_SYMBOLS.md`
（SHA256 `b15ee914ffeb875519d66cedb703743e8ca7c0d8f5fcb44f1e0e75eead3a55ef`）を
別作者として全読解し、以下の係数を独立に手導出した。**必要還元としてPASS**。
新point、operator、matrix、CAS、jet又はGBは実行していない。
残五symbolの全higher sourceの存在、像次元、成分性及びE帰属は未結論である。

## 1. 全irregular connectionから全P0へ

actualDを被約既約quintic、restriction planeのzero cubicをf=a³、normal位数を二、
normal quadraticをgeneric squarefree、primitive root degreesを(0,1)とする。
既独立受理 `PLANE_NONREDUCED_M2_SPLIT_ROOT_ENTRANCE.md` のpole budgetを用いる。
root rowのprimitive unitからcharacter poleを先に評価しているので、
root basepoint、非被約Rの高位極、exact部分をlog-only sourceに落とさない。

全root sourceはconstant normal Borel及びactual base座標により

```
G2=gamma*a²*u*(b*u+a*v), gamma!=0,
M=[[1,0],[b,a]],
eta=theta1-theta2=-dlog a+d(L/a), L=r*b+s*c,
theta=theta1=dlog a+d(P2/a²).
```

ここでP2 modulo Ca²は全五exact係数を持つ。
`beta=M^-1 diag(theta1,theta2)M-M^-1dM` を直接掛けると

```
beta11=theta,
beta22=theta-d(L/a), beta12=0,
beta21=-d(b/a)-(b/a)*d(L/a).
```

全P0はbase P_f、mixed `uX` と `uZ+vY`、normal `k*u²+kappa*u*v`
である。targetのv²項零からnormal v²係数零を得るので、先にKをsliceしていない。
tracefree mixed gaugeでclosed diagonalからXa=Ya=0となる。
`iE iXvol/f=dlog f-3beta11` と他の三entryを戻すと

```
Y-X=(0,-3s*a,3r*a),
Z=(s*a,-2s*b,s*c+3a+3r*b).
```

後者のcontractionは厳密に
`(-3b(a+L),3a²+3rab,3sab)` であり、`-3f beta21` と一致する。
tracefree条件はcontraction kernelのEuler multipleも除く。

全未比較Xを `(0,A*b+B*c+p*a,C*b-A*c+q*a)` と置く。
actualDのu²、uv係数は、座標順(a,b,c,u,v)のEuler-signed Pfaffianから

```
u²: 3f*k-det(E,X,Z)-Z(f),
uv: 3f*kappa-det(E,X,Y)+X(f)-Y(f).
```

ここではX(f)=Y(f)=0、Z(f)=3s*fであり、直接のdet展開は

```
det(E,X,Y)=3a²*(r*Xb+s*Xc),
det(E,X,Z)=3a*((a+r*b)*Xb+s*b*Xc).
```

全u²係数のa³,a²b,a²c,ab²,abcと全uv係数をtargetへ比較すると、順に

```
B=0, k=p+s, A=-kappa,
s*kappa=0, r*kappa=s*C,
gamma=3*(kappa-r*p-s*q)!=0
```

を得る。したがってproducerの必要P0表示は全係数比較と一致する。
P1三十／raw三十六modulo Euler六、P2九をここで消していない。

## 2. kappa非零の全sourceは空

kappa非零ならr=s=0、従ってY=Xである。
全P1のu-columnのu非可除項をv²Cbadと書く。
Euler六方向はこのCbadを変えず、P2はnormal columnを持たない。
Cbad零ならfull u-rowはuで割れるためactualDにu因子が生じる。
よって対象ではCbad非零が必要である。

既独立確認済みの元projective Jacobiのcurl恒等式
`(2kappa I+4DY-DX)Cbad=0` は、今 `Y=X` なので
`(2kappa I+3DX)Cbad=0` となる。
DXはC,p,qを保ったtriangular matrixで固有値0,-kappa,kappaを持つ。
従ってこの行列の固有値は2kappa,-kappa,5kappaで全て非零となり、矛盾する。
ここでnoncanonical raw Jacobi自体をstrictに零と仮定していない。

## 3. 全higherを保つratio-Casimir operator

以後kappa零。canonical liftの正当性を
`PLANE_M3_CANONICAL_TRANSVERSE_OFFDIAGONAL_OBSTRUCTION.md` §2のBV恒等式で再確認した。
quadratic lift Pに対し `Pbar=P-Etotal wedge div(P)/5` はdivergence零で、
元projective Jacobiからstrict Jacobi零が従う。
rank四の場合、五変数volumeでdualしたPbar²のclosed homogeneous one-formは
actualD Fの非零constant倍のdFであり、Pbar sharp dF=0である。

全P0に対してdiv(P0)のbase成分は-(X+Y)、normal成分は(0,2ku)である。
従ってcanonical symbolは

```
Xbar=(4X-Y)/5, Ybar=(4Y-X)/5,
Zbar=Z-(2k/5)Ebase, Kbar=(3k/5)u².
```

zeta=v/uに対するHamiltonian fieldは正確に

```
base: A(zeta)-(2k/5)Ebase+B(u,v),
normal: (3k/5)(u partial_u+v partial_v),
A=Z-zeta*(X-Y).
```

Bは全P1に由来するnormal degree一の全rational columnである。
P2はnormal columnを持たないのでこの表示から落ちるが、Fにはその全寄与を残す。
Aはzetaを保つ。field C(u,v)上のbase translationでBを消す条件は
`det(k I-A)!=0` であり、translationをactual PGLとは呼ばない。
translationはnormal degree一なので、Fのnormal degree二のleading部分を変えない。
normal degree m、base degree d=5-mの全項に対するoperatorは

```
Sym^d(A)+k*(m-2) I.
```

scalarは `(-2k/5)(5-m)+(3k/5)m=k(m-2)` から直接得る。
したがって3/5と2/5を別のnormal degree conventionで取り違えていない。

全m=3,4,5のoperatorが可逆なら、translated Fは元leading
`gamma*a²*u*(b*u+a*v)` のみである。
これはC(u,v)[a,b,c]内でpositive base degreeの非単元因子に分解する。
元の既約FはC[u,v]上primitiveでありGauss lemmaによりC(u,v)上も既約なので矛盾する。
この議論はdense rational座標での必要式であり、非regular PGL同一視やformal jetの代用ではない。

## 4. s零の全complementと五つのspectral壁

kappa=s=0ならgamma=-3rp非零からr,p非零、k=p非零。
Aは `3*((1+r*zeta)*a+r*b) partial_c` で、C,qを保ったままnilpotentである。
translation operator kI-A及び全m=3,4,5のoperatorは可逆となる。
よってこの全sourceは上のGauss矛盾で空である。

s非零ならC=0。Aはa,b,c順のtriangular matrixで固有値(s,-2s,s)を持つ。
Jordan blockの消滅を仮定せず、degree dでの全固有値はs(d-3j)、0<=j<=dである。
operatorの非可逆値を全て戻すと

| normal degree m | base degree d | k/sの全候補 |
|---|---:|---|
| 3 | 2 | -2,1,4 |
| 4 | 1 | -1/2,1 |
| 5 | 0 | 0 |

translation failureの1,-2もこの和に含まれる。
従って必要値は正確に `k/s in {-2,1,4,-1/2,0}`、
すなわち `p/s in {-3,0,3,-3/2,-1}` である。
k零及び全resonanceを無断で反転していない。

## 5. actual normalizerと全higherの戻し

old c=c+t*bのactual base shearはr'=r+s*t、q'=q-t*pを与える。
t=-r/sでr'=0、q'=(rp+sq)/s非零となり、逆はminus shearである。
root row b*u+a*vもこの操作で保たれる。

old base=(alpha*a,beta*b,chi*c)、old normal=(mu*u,nu*v)、
nu=beta*mu/alpha、tensor scale=beta*chi/alpha²とすると、base P_fの係数3を保ち

```
p'=chi*p/alpha, s'=chi*s/alpha,
q'=beta*q/alpha, r'=beta*r/alpha.
```

Zaの係数、Zcの3a係数、normal kの変換もこの同じ式と一致する。
chi=alpha/s、beta=alpha/q'でs'=q'=1、r'=0へ戻す。
反転するのはこの分岐で非零と証明したs,q'、選択したcoordinate unitsのみである。
全P1/P2は同じinvertible actual変換で同時に運ばれる。

従って必要残余は全て

```
p in {-3,0,3,-3/2,-1},
X=(0,p*a,a), Y=(0,(p-3)*a,a),
Z=(a,-2b,c+3a), K=(p+1)u²,
G2=-3a²*u*(b*u+a*v)
```

の五固定P0の全higher sourceへ入る。これは必要入口の有限化であり、
五sourceの非空性、全G上界、全E帰属又は新成分数は結論しない。
SF00、既空SF02、normal double、f=a²b、m三及び全分類は別義務として保持する。
