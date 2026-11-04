# R-primary二重DVRと残余tailの間の有限gluing：長さ16/19の上界

2026-11-04。canonical separationで得た

```
S_F=0:_{D0}F^infinity, B_R=D0/S_F=QQ[[U,epsilon]]/(epsilon^2)
```

を用いる。F=v^2-nu*z1^2、U=rho112は実際の同R0 scalarである。
J_R=image(T0→B_R)、B0=image(D0→D4)と書く。
本稿は純module論で、新しい係数展開・source jet・点・prime・GBを使わない。

## 1. F²はwhole D4をannihilateする

actual W1 tailの被約環はNで、nilradicalはsquare-zeroである。
FはNへのprojectionのkernelに入るので、D4内のFはそのnilradicalに入る。
従って**F² D4=0**。これはnilのgeneric長さだけでなく、既受理の全tail完成環を使う。

D0→D4のkernelがT0だからF²D0⊂T0であり、B_Rへ射影すると

```
F² B_R subset J_R.
```

## 2. 有限長上界16

B_RのQQ[[U]]-basisを(1,epsilon)とし、F=f(U)+g(U)epsilonと書く。
既存normの計算からord_U f=4。epsilon²=0なので、F²の乗法matrixは

```
[[f², 0], [2fg, f²]]
```

であり、determinantはf^4、U-orderは16である。
DVR上の有限free moduleのdeterminant length公式から

```
length(B_R/F² B_R)=16,
length(B_R/J_R) <= 16.
```

## 3. Contact4から得るのは包含のみ

D4→Nへの写像を通すとT0は零になる。Nのzero-degree subalgebraはN tailへ単射に入るので、
T0はactual D0→N0のkernelにも含まれる。
Rは全次数でtorus-fixedであるため、受理済みの全N∩Rの長さ4は、そのdegree0制限でも同じである。
positive-weight方程式のR制限は恒等零なのでcontact idealを変えない。
したがって

```
image(J_R -> B_R_red=QQ[[U]]) subset (U^4).
```

この像が(U^4)と等しいとはまだ主張しない。

J_Rは有限colength idealなので、DVRのHermite形により

```
J_R=(U^m+b(U)epsilon, U^n epsilon),   b modulo U^n.
```

epsilon倍で閉じるためn≤m。上のcontact包含からm≥4、F²の包含からm≤8。
colengthはm+nであり、全体として

```
4 <= m <= 8,  0 <= n <= m,  m+n <= 16.
```

b、m、n、元R0のepsilonへの作用は未計算である。

## 4. 残余F-torsionがtail内で欠く部分は長さ19以下

D0/S_FからD4/image(S_F)への誘導写像のkernelはJ_Rである。
さらにD4/B0=Zで、既受理のexact同定はZ=omega_(QQ[[r,s]]/(r,s)^2)、長さ3である。
したがって正確な列

```
0 -> B_R/J_R -> D4/image(S_F) -> Z -> 0
```

があり、後者は**有限長19以下**となる。
これはimage(S_F)がtail内でfinite-colength moduleであることを与える。

しかしker(S_F→D4)=T0 intersect S_Fはまだ決まらない。
S_F自身やこのkernelが有限長であるとは限らず、そのhidden socleの有無も未決である。
有限cokernel boundを有限kernelや元Dのdepthへ取り違えない。

Code `code/ambient_B18_T0_finite_R_gluing_bound_20261104_v63.py`。
Report `reports/AMBIENT_B18_T0_FINITE_R_GLUING_BOUND_20261104_v63/report.json`。
