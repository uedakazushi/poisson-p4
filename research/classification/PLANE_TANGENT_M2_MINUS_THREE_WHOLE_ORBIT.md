# tangent m二・t=−3：残二parameter全源も一実軌道へ

2026-11-08、`/root/quartic_review`。純構造証明、新しいM1／ambient／profile計算無し。
`PLANE_M2_TANGENT_MOVING_WHOLE_SOURCE_REDUCTION.md` §3の全(r,s)を扱う。
s!=0は一actual PGL軌道、s0はその閉包に入り、Eの残義務を一canonical tensorへ縮める。
既知成分名・成分性・軌道次元の下界は本稿では主張しない。

## 1. Fullsourceとcanonical Hamiltonian

座標(a,b,c,u,v)で全源は

    Pab=3(ac−b²), Pac=a(ru+sv),
    Pau=−3(bu+cv), Pav=−3(au+bv),
    Pbc=asu+crv, Pbu=v(ru+sv),
    Pbv=−3(bu+cv), Pcu=Pcv=Puv=0.

全(r,s)で元350式は恒等零、全9P2は元weight二系から零と一意に決まる。
行発散R_i=Σ_j∂jPijのnormal成分はRu=0,Rv=6u。
canonical divergence-zero lift Π=P+(E∧R)/5はEuler ambiguityに依らない。
ρ=v/u、M=Π#dρとするとdense field上で

    M=P#dρ+(6/5)E,
    Ma=−9a/5+3cρ²,
    Mb=−9b/5−3cρ−u(rρ²+sρ³),
    Mc=6c/5, Mu=6u/5, Mv=6v/5.                         (1)

これは全higher coefficientsを戻したexact expressionで、unknownを零にしたものではない。

## 2. 必要な全actual coordinate map

λを任意として次の六定数（source r,s,λのpolynomial）を置く。

    f=−2λs/3, e=−(2λr+4λ²s)/3,
    A=λ²e, B=−2λe+λ²f,
    D=−e−λf+2λr/3+λ²s,
    C=(rλ²+sλ³)/3−λD.

actual PGL変換を

    a'=a+2λb+λ²c+A u+Bv,
    b'=b+λc+C u+Dv,
    c'=c+e u+fv,
    u'=u, v'=v−λu                                      (2)

とする。det=1。逆はu=u',v=v'+λu'から順にc,b,aをback-substituteする明示三角逆。
λ又はrを反転する必要は無い。
bare Borel部分は元P0をliteralに保ち、残六shiftはpositive-weight actual shearなので、
全Poisson sourceとassociated-graded P0を保つ。

ρ'=ρ−λとすると、(1)を(2)でpushしたMa',Mb'を直接整理して

    Ma'=−9a'/5+3c'ρ'²,
    Mb'=−9b'/5−3c'ρ'−u'[(r+5λs)ρ'²+sρ'³]

を得る。残成分も(1)と同形である。
六定数は順にMa'のρ³,ρ²,ρ,定数項、Mb'のρ,定数項を消した全解。

## 3. Mだけの一致からtensorを推測しないための完全性

旧全K1 dimension八は、実六shearと表示r,sの二方向で完全に張られる。
この全同型及び全9P2のrank9は前のfinite receiptに保存されている。
canonical Hamiltonian Mへのweight一写像はこの全八次元上でinjectiveである。
実際、六shiftを(A,B,C,D,e,f)と書くと、Maの変化は

    3u[A+Bρ−eρ²−fρ³]。

従ってρ²/ρ³からe,f、定数/ρからA,Bを一意に読む。
続くMbのρ²/ρ³からresidual r,s、定数/ρからC,Dを読む。
これはfullkernelの全generatorを含む写像であり、特定Hamiltonian一致を全tensor一致と仮定していない。

(2)の像を既定Euler quotientへ戻すと同じwhole sourceに入る。
Mの一致と上記injectivityでP1がP1(r+5λs,s)に一意に決まり、
全9P2は元rank9 affine系で一意なのでP2=0となる。
したがって正確にprojective tensorとして

    gλ_*P(r,s)=P(r+5λs,s).                              (3)

ρ chartをPGL変換と誤認していない。(2)自身は全P4上のpolynomial線形可逆mapである。
dense u!=0上で得た(3)は正則tensor sectionの等式だからu=0にも延びる。
u0上の点又はsource componentを捨てたものではない。

## 4. 全unit枝・零枝・canonical既約性

s!=0ではλ=−r/(5s)でr'=0。
共通normal scaling(u,v)→(su,sv)はP0を保ち、residual(r,s)→(r/s,1)。
したがって全unit枝はP*=P(0,1)の一actual軌道へ入る。
全inverseは上記三角逆とs^-1で、source unitはsだけ。

s=0ではactual Dは

    −[r a u²+3c(bu+cv)] [3(bu+cv)+r v²]

と分解する。さらに全(r,s)はaffine polynomial Poisson族、s!=0は稠密なので、
全s0も同じ軌道閉包に入る。canonical点の既知帰属が分かれば全源へ閉性で戻せる。

actual Dの既約性は先に確認出来る。
K=C(c,u,v)上で

    A0=a−cρ², B0=b+cρ+u(rρ²+sρ³)/3

とすると、実Dの全式は

    F=3u²[su A0²−(ru+sv)A0B0−3cB0²].                  (4)

canonical r0,s1ではquadratic discriminantはv²+12cu。
これはrank三irreducible quadricで、UFD valuationが一なのでK内のsquareではない。
従って(4)はK[a,b]で既約。
元のpolynomial Fはa²係数3u³とa係数の−v⁴項からcontent1で、
Gauss lemmaによりC[a,b,c,u,v]上でも絶対既約。actual F!=0なのでgeneric rank四。

## 5. 第二planeを無断に仮定しない

−9/2で使えた(c=u0)は、ここではP_ac=avが残りPoisson planeではない。
同じconcurrent-m三定理をそのまま適用しない。
canonical u1 chartでA0=a−cv²,B0=b+cv+v³/3とすると

    F=3(A0²−vA0B0−3cB0²).

このchartのSing FはA0=B0=0。
そこでは全affine bracketsが零で、rankzero surfaceである。
その射影閉包は

    [a:b:c:u:v]=[3CV²:−3CUV−V³:3CU²:3U³:3U²V]

とparametrizeされる。これはknown-parent認識への構造的入口であり、
parametrization又はdivisor形だけからT40等への帰属を結論しない。
次に行うなら、この**一whole orbit**のactual known-family認識を最弱目標とする。
