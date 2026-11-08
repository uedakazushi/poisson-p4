# 残るnode/cusp-flex二源：cubic一般leafのadditive residue障害

2026-11-08、PRIMARY E、`/root/quartic_sources`。全磁気係数を保持した純手証明。新CAS、point、operator、rank/profile、GB又は局所変形jetは無い。target上のno-algebraic-leafを独立に証明する。whole T37への適用は、凍結した共通proper leaf定理の受理を別依存として明示する。

## 1. 全源とintrinsicな関数体計算

`TWO_LINES_PLANE_CUBIC_THREE_MAGNETIC_QUOTIENTS.md` と独立査読により、元六磁気係数は全域actual shearで `R=Ih²+Jhy` へ移る。I又はJを反転せず全零補集合を保持する。h=1上で

```
P=U wedge partial_z+V wedge partial_w+(-d*w+R)*partial_z wedge partial_w.
```

nodeとcuspのC,U,V,dは `TWO_LINES_M3_WHOLE_T36_QUADRATIC_PENCIL_EXCLUSION.md` §1に印刷した全式であり、両Cは既約。各cubicのnormalization parameter tauを使い、一般開で

```
U=a(tau)*partial_tau, V=b(tau)*partial_tau,
[U,V]=d(tau)*V.
```

と書く。このnormalization／有理座標はintrinsic leafの計算にだけ用い、actual PGL mapやcomponent所属の証明とは呼ばない。

cubic上のvolume d tau wedge dz wedge dwに対するkernel一形式は

```
alpha=(-d*w+R)*d tau-b*dz+a*dw.
```

k=b/a、m=d/aと置くとJacobi関係は k'=mk。u=w-kz、mu'/mu=-m、s=mu*uを使えば

```
(mu/a)*alpha=ds+(mu*R/a)*d tau.                    (1)
```

z係数はk'-mkによって正確に相殺する。aやmuの零・極を避ける一般開だけで議論するが、最後のrational differentialはP1全体のresidueを持つ。一般点のalgebraic leafを否定するにはこの開で十分である。

## 2. Node-flex/tangent-through-flexの全I,J

```
C=y²-x(x-1)²,
x=tau², y=tau(tau²-1),
a=tau(tau²-1), b=-tau*a, d=tau²-1,
k=-tau, mu=1/tau, s=w/tau+z.
```

dx/d tau=2tauを使えばU_x/(2tau)=a,V_x/(2tau)=bであり、yの微分でも同じ式を得る。R=I+Jyだから

```
mu*R/a=J/tau+I/[tau²(tau²-1)]
       =J/tau-I/tau²+I/[2(tau-1)]-I/[2(tau+1)].
```

S=s+I/tauとすれば(1)は

```
Omega_node=dS+J*dlog(tau)
            +(I/2)*dlog((tau-1)/(tau+1)).           (2)
```

residueはtau=0でJ、tau=1でI/2、tau=-1で-I/2、無限遠で-J。従って(I,J)≠(0,0)なら少なくとも一つが非零。特にJ=1,I=lambdaの**全finite lambda**でtau=0のresidueが一である。

## 3. Cusp-flex/secant-through-flexの全I,J

```
C=y²-(x+1)³,
x=tau²-1, y=tau³,
a=tau(tau²-1), b=tau*a, d=tau²-1,
k=tau, mu=1/tau, s=w/tau-z.
```

ここでも全U/V係数とJacobi関係に一致する。R=I+J tau³なので

```
mu*R/a=J*tau/(tau²-1)+I/[tau²(tau²-1)].
```

S=s+I/tauとすれば

```
Omega_cusp=dS+((J+I)/2)*dlog(tau-1)
             +((J-I)/2)*dlog(tau+1).              (3)
```

residueは(J+I)/2,(J-I)/2、無限遠で-J。全て零となるのはI=J=0だけ。J=1の全lambdaでは有限二residueの和が一で、lambda=±1でも片方は一である。

## 4. Nonzero additive residueは有限代数coverでも除けない

有理一形式omega(tau)d tauがP1上で非零residueを持つとし、滑らかな一般三次元開でdistributionが

```
partial_z, D=partial_tau-omega(tau)*partial_S
```

により張られるとする。不変irreducible代数面があれば、C(tau)[S,z]上のheight一primeをirreducible polynomial Fで生成出来る。partial_z invarianceとz次数からFはzに依らない。tauが定数というvertical面はD tau=1に反するので、FはSについて正次数を持つ。

従って面の関数体でSはC(tau)上代数的であり、有限拡大L/C(tau)にDが延び、DS=-omegaとなる。char零でtraceとderivationは可換なので

```
d Tr_(L/C(tau))(S)=-[L:C(tau)]*omega(tau)*d tau.
```

左辺はrational functionの微分であり、P1の全点でresidue零。一方右辺には非零residueがある。矛盾。この議論はramified coverや任意の複雑なalgebraic leafを排除し、rational first integralだけを試した陰性ではない。計算用開から除いたtauの点のresidueも、関数体のrational primitiveの必要条件として有効である。

## 5. 全sourceでの結論と保持する補集合

(2)(3)により、両sourceの(I,J)≠(0,0)部分では、cubic成分の一般rank二symplectic leafはその三次元cubic成分にZariski denseである。その一般点を含むinvariant代数面は次数を問わず存在しない。

I=J=0ではOmega=dSでalgebraic first integralが実在する。この零磁気cokernelを陰性に含めない。Jを反転せず、J=0,I≠0もresidueで直接扱った。元六磁気はactual shear正逆で戻せるので、この非零条件は二つの線形cokernel invariantの全域条件である。

共通定理 `T37_WHOLE_CLOSED_LOW_DEGREE_LEAF_SURFACE_THEOREM.md` が別作者により受理されれば、whole T37のrank四境界のactualD各componentの一般smooth rank二点はdegree≤5のinvariant surfaceに覆われる。この定理と本稿を合成すると、**両sourceの全(I,J)≠(0,0)点はwhole T37外**、従って両whole23 sourceの親としてT37は不可能となる。

これは共通定理のproper Hilbert／marked一般点／flat ideal延長を必要な依存として使う。単にT37に一枚だけleafが存在するという弱い主張では足りない。special limitがdouble H等の別componentへ逃げないこと、embedded点だけをtarget一般rank二点へ残せないことも、その共通定理の一部である。

このproducer保存時には共通定理の別作者receiptを待っており、global新成分数・相互相違・全分類完了は先取りしない。旧二次/三次局所零結果、固定frameの失敗、normal-cone全u因子証明は不変のまま、本稿は全leafの別global幾何による判定である。
