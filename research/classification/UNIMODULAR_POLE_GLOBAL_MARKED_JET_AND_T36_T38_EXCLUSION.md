# h5残軌道のglobal親排除：proper marked-hyperplane jet条件とT36/T38

2026-11-08、PRIMARY `E_multiple_hyperplane`。純証明、独立root読解前。
新operator／point／rank profile／局所変形jetの延長はない。
同一点全D101／M47とwhole source一軌道還元は既immutable証明書を入力とする。
本稿の「二次jet」は既tensorの固定幾何的不変量であり、変形次数を増やす計算ではない。

## 1. marked Poisson hyperplaneの全quadratic lift

projective bivector πのPoisson hyperplane ℓ=0をmarkし、局所regular PGL frameでℓ=hとする。
任意のraw homogeneous quadratic lift Pについて、H上のnormal成分は
E4∧(P^{hi}|h0)=0を満たす。従ってUFDで

```
P^{hi}|h0 = q(x) x_i
```

となるlinear qがある。L^h=−qを取り、L^i=(P^{hi}+x_i L^h)/hとすれば全L^iはlinearで、
P−E5∧Lのh行は零である。このh-Casimir quadratic liftは一意である。
実際h行零の二liftの差E5∧Aにはh A^i−x_i A^h=0が必要で、A^h=ch、A^i=cx_iとなり差は零。
これは局所regular frameでlinear係数操作だから、全marked-Poisson-hyperplane incidenceの相対vector bundleとして構成出来る。

h一chartで一意のquadratic affine tensorをC+L+Qと分ける。
T_Q(v)=dQ(v)の24×4行列にLを加えた24×5行列に対する

```
rank[T_Q | L] <=4                                      (1)
```

はclosed条件。linear frameとtensor scalarは可逆変換、translationはL列にT_Q(a)を加えるだけなので、
(1)はaffine frameに依存しない。rankT_Q低下も許す弱いdeterminantal条件であり、
Lが像に含まれるという非closedになり得る条件を使っていない。

actualD非零rank四openで、ℓ|Dかつℓ-Poissonを満たすmarked incidenceは
Poisson parameter×dualP4の相対closed部分で、ℓを忘れる射はproper。
(1)はこのincidence内でclosedにglueするので、そのforgetful像K_evenはrank四open内でclosed。
したがって全source C+QのPGL飽和閉包のrank四部分はK_evenに入る。
sourceごとのmarkedℓは動いてよく、極限で他factorと合流してもpropernessがmarkを保持する。
sourceのgenericD=h5を仮定していない点が、zero-log限定の前稿からの追加である。

P*のD=h5ではmarkの唯一候補はh。前稿の明示minorによりそこでrank[T_Q|L]=5。
従って**P*は全C+Q斉次化sourceのPGL飽和閉包から外れる。**

## 2. 適用する追加既知成分の正確なsource

`code/additional_rank4_families_20261103_v62.py`の式と旧各component定義を用いる。
次の成分は、その全定義sourceがh行零のh²J+Qであるから§1が適用する。

|成分|全source入力|
|---|---|
|R4_TRIPLE_HYPERPLANE_QUADRIC|legacy/v05/triple_hyperplane.tex|
|R4_H_QUARTIC_ESSENTIAL5|legacy/v06/nilpotent22.tex|
|R4_ZERO_LOG_5H|legacy/v06/zero_log_component.tex|
|R4_MINIMAL_M|legacy/v07/minimal_thickness.tex|
|R4_A0|legacy/v07/nil22_v7.tex|
|R4_DPM|legacy/v07/dplus_fullp4.tex|

fullDPMの二つの共役sign、全Sp4/GL/scaleと各source境界も同じclosed incidenceに入る。
別sourceが同成分へdenseに入るという推測は使わず、componentの定義sourceを使っている。

R4_CMINUS_MAGはgeneric γh zのlinear磁気項を持つので、L0とは呼ばない。
しかし全sourceの最高Qはyに依存しないためrankT_Q≤3、任意Lを追加してもrank≤4。
従ってその全source／PGL飽和閉包もK_evenに入り、P*を含まない。

R4_CPLUS_EXTはこの議論で扱わない。旧全familyのgeneric projective tangentは25であり、
その全成分上ではrankD≤100というclosed条件が成立する。
P*のrankD101により同成分も親ではない。
他のMPS sourceの(1)判定は別pure inventoryへ分け、本稿は全MPSを排除しない。

## 3. T38のwhole closed zero-surface条件

`GA9_LINEAR_W_WHOLE_T38_AND_A5_EXCLUSIONS.md` §§2–3がtrue completed T38の全generic源から
rank-zero twisted-cubic cone surfaceを明示し、universal projective zero incidenceのpropernessで
全T38へdimZero≥2を延長した。全period／edge零補集合もこの閉性に含まれる。
前稿のP* zero台はconic＋二直線で正確に一次元だから、**P*はwhole T38の外**。
T36はgeneric zero台が曲線であり、このsurface議論を適用しない。

## 4. P*のorder二以上zeroは一点だけ

全zero台は
C={h=y=0,2xw+z²=0},ell1={h=x=y=0},ell2={h=x=z=0}。
共通点p=[h:x:y:z:w]=[0:0:0:0:1]を除くと、各既存tensorのchart firstjetに以下の非零行がある。

* Cのx一chart：π_hz=h（線形係数一）。
* ell1のz一chart：π_xw=x+高次。
* ell2のy一chart：π_zw=x+h+高次。

よってp以外にはorder二zeroは無い。pのw一chartの全0/1jetは零。
局所順(h,x,y,z)のquadratic leading tensorは全六式で

```
Q=(0,h²,−hx,2hx,−2x²,−h²−hz).                         (2)
```

これは元fullP*をw一へdehomogenizeして二次項を取った固定式である。
新local deformation ansatzではない。

## 5. pのcubic one-formはquadric-pencil cone外

`GA9_LINEAR_W_WHOLE_T36_QUADRATIC_PENCIL_EXCLUSION.md` §§1–3のwhole定理を使う。
T36の全closed成分は「あるorder二zero qで、α_Q=i_(E4∧Q)vol4が
quadric-pencil imageのclosed affine cone C22（零を含む）に入る」を満たす。
quadric-pencil map Gr(2,Sym²V*)→P(H0Omega1(4)), <F,G>→[GdF−FdG]は全域regular/properなので、
共通factorを持つpencilも含む。この全scopeを狭めない。

(2)から全cubic one-formは

```
α = x(−h²+hz+2xy) dh + h(h²−xy) dx −hx² dy −h²x dz.
```

特に非零。quadratic jetの方向射影空間のx一chartでは

```
α=(−h²+hz+2y)dh−hdy−h²dz.
```

s=y+hzとすればα=(2s−h²)dh−hds。
このfoliationのrational first integralがあれば、まずkernel方向∂z−h∂yを消去して
C(h,s)に属する。hunitでv=s/h²とすれば、残るkernel derivativeは

```
D=h∂h−∂v  on C(h,v).
```

この有理関数体でDの定数体はCである。実際T=A/B、A,B∈C[h,v] coprimeについてDT0なら
A|DA,B|DB、両cofactorは同じpolynomial f。
Dはtotal polynomial degreeを上げないからfはconstant。
A=Σ h^k a_k(v)とするとDA=fAはa'_k=(k−f)a_k。
非零polynomial a_kがあるのはf=kかつa_k constantの場合だけ。
従ってAとBはいずれも同じh冪のscalar倍で、coprime条件からTはconstant。

よってαには非定数rational first integralが無い。
特に独立quadric ratio F/Gをfirst integralに持つC22の非零元ではない。
pは唯一order二zeroなので、P*はT36のglobal closed conditionを満たさず、
**P*はwhole T36の外**。

## 6. 残す義務とscope

P*のwhole actual orbitの既知親はまだ同定していない。
全h5還元をGのnon-dominanceと取り違えず、この一軌道を未処理として残す。
T19、T37又は他のMPS sourceについて本稿は未証明の包含／排除を足さない。
既ambientのgap二を根拠に別点／higher jet／minor系列を再開しない。
