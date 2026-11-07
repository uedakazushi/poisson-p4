# rank一nil共通場のK零残余：suspensionと有限Q orbit

2026-11-07、主義務G。producer手証明。全恒等式とfixed-Q linear ranksは次の有限jobで別に保存する。
旧dependent six-source linear-only reportのfalseと大きい上界は変更しない。

## 1. a非依存liftは全kernelから実際に構成出来る

残るsymbolはV=N=b∂a、K=0である。
normalssはT=tNの全t-line、normalnilはT=N又はT=0。
`CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_DEPENDENT_20261107.json`の全D1 kernelsは、
base-base bc成分にa係数を持たず、a依存項が`a∂a∧L(y)`だけである。
ここでy=(b,c,z,w)、Lはnormal-linearなbase vector field。
full Π2九係数もa非依存なので、Euler_total∧Lを引くと同じprojective tensorを

`Π'=Q(y)+∂a∧W(y)`                                      (1)

というa非依存quadratic liftで表せる。
これはdivergence又はglobal stabilizerの仮定ではなく、full kernelのliteral Euler gaugeである。
Π2のbc三係数f2(z,w)も**全てQに保持**する。
残るab/ac六係数及びcommon Π1十二係数がWの全十八変数となる。

全projective Jacobiをaの次数で比較すると

`[Q,Q]=0`, `E_y∧[Q,W]=0`。                              (2)

逆も成立する。従ってQを固定すれば、全original JacobiはWへのinhomogeneous linear systemになる。
Wのpositive-weight raw二十係数のうち`W→W+f(z,w)E_y`というEuler二冗長性を固定し、
W_zのz²係数、W_wのw²係数を零として十八係数にする。
このgaugeは(2)とactual Dを変えない。

## 2. constant factorのsuspension lemma

Qが`Q=v∧Z`、vがyのnonzero constant field、`L_vQ=0`とする。
PfQ=0なのでactual Dはa非依存で、符号を除き`E_y∧W∧Q`の係数である。
そのv微分は`E_y∧[v,W]∧Q`。
(2)をvでwedgeすると、`[Q,W]=v∧[Z,W]+[v,W]∧Z`の前項が消え、ちょうどこの微分が零となる。
従ってDはaとvの二独立constant directionsに依存せず、本質高々三変数。
Gの既約被約quinticなら、本質二以下はbinary factorizationに反するので本質三である。
既存 `ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md`及び別作者reviewの全source上界二十三を適用出来る。
このlemmaはfullPoisson fiberの次元を小さいと仮定するものではない。

## 3. normalssの通常kernelとnormalnil trace一はconstant factorへ入る

normalssでp=t+1,q=t−1とする。
通常のfull十六kernel（t=±1,±3以外。±1/3ではfull rank十六を別に確認済み）の非共通四係数を
B,E∈C²と書く。Qのnormal anchorは

`X_B=tzw∂z+q w²∂w`, `X_E=−p z²∂z−tzw∂w`。

Q^{bc}のbase-linear部分は`((q/2)wB−(p/2)zE)∧(b∂b+c∂c)`、
残るf2は全binary quadratic。
strict Q JacobiはΔ=det(B,E)について`(3t+1)Δ=(3t−1)Δ=0`を与えるので、全tでΔ=0。
B,Eをcommon constant vに沿わせると、全Qはvをfactorに持ち、L_vQ=0。
全f2を含め§2へ入る。

normalnil trace一では非共通四係数D,E∈C²で、anchorsは

`X_D=2z²∂z+z(z+2w)∂w`, `X_E=z(z+2w)∂z+2w(z+w)∂w`。

base-linear partは`(zD+(z+w)E)∧(b∂b+c∂c)`。
strict Jacobiはdet(D,E)=0を強制するため、同様に全sourceが§2へ入る。
rank二の大きいlinear sourceを数えるだけの上界三十三を、新しいsource geometryで置き換えた。

## 4. normalnil trace零の全七cohomology

非共通係数はD,B,E∈C²及びfである。
normal anchorsは

`X_D=z²∂w`,

`X_B=−2zw∂z−3w²∂w`,

`X_E=z²∂z+2zw∂w`, `X_F=z²∂z+zw∂w`。

b,cのanchorを`X_b=D_bX_D+B_bX_B+E_bX_E`、
`X_c=D_cX_D+B_cX_B+E_cX_E+fX_F`とし、base-linear partは
`(−2wB+zE)∧(b∂b+c∂c)`。
全f2を保持したstrict Jacobiは正確に

`2det(D,B)+f E_b=0`、`det(B,E)+f B_b=0`、`det(D,E)−f D_b=0`。  (3)

f=0なら三vectorsはcollinear、全first column零でもFを含め全vectorsが∂c-lineなので§2へ入る。
f非零かつfirst column非零の場合、Plücker恒等式と(3)から

`4D_bB_b+E_b²=0`。

actual normal shear w→w+h zは

`(D,B,E,f)→(D−hE−h²B,B,E+2hB,f)`。

B_b非零ならh=−E_b/(2B_b)でD_b=E_b=0へ直せる。
base shear c→c+s bで残るB_cを消せる。
B_b零ならconic式からE_b零、D_b非零であり、c-shearでD_cを消せる。
base a,bの同時scaleとc-scaleでfirst nonzero entryとfを1へ直し、全sourceは次の二anchor orbitへ入る。

| 型 | D | B | E | f | Q^{bc} |
|---|---|---|---|---|---|
| A | (1,0) | 0 | (0,1) | 1 | −zb+f2 |
| B | 0 | (1,0) | (0,−1) | 1 | zb−2wc+f2 |

Aのanchorsは`X_b=z²∂w`, `X_c=2z²∂z+3zw∂w`。
Bのanchorsは`X_b=−2zw∂z−3w²∂w`, `X_c=−zw∂w`。
上の全変換はP0のweight零symbolを保つactual changesで、逆も明記したnonzero scalars又はminus shear。
first columnを選んだ二つの標本ではなく、nilpotent conicの全非零支持を覆う。

## 5. normalssの全critical periods

normal permutation z↔wとbase a→−aはactualにt↔−tを実現する。
従ってt=1,3だけを扱えば全四critical valuesを覆う。

t=3では通常B,Eに加え、Qのc-anchorへ`f w²∂z`が加わる。
strict Jacobiは`det(B,E)=0`, `f E_b=0`。
f零又はB_b零のbranchは§2へ入る。
残る全sourceはE=0、B_b非零、f非零であり、actual Borel base changesとscalesでB=(1,0),f=1へ直る。
canonical anchorsとbase bracketは

`X_b=3zw∂z+2w²∂w`, `X_c=w²∂z`, `Q^{bc}=wc+f2`。        (S3)

t=1ではextra c-anchorは`f(zw∂z+w²∂w)`。
strict Jacobiは`det(B,E)=f E_b`, `f B_b=0`。
constant-factor branch以外はB=(0,−f)、E_b非零で、同じactual changesでE=(1,0),f=1へ直せる。
canonical formは

`X_b=−2z²∂z−zw∂w`, `X_c=w²∂w`, `Q^{bc}=−zc+f2`。       (S1)

全f零、first-column零、intersectionをconstant-factor branchに保持し、open nonzero rowsだけから全supportを推論していない。

## 6. 全central extension f2も0/1へ直る

actual shear b→b+u(z,w)、c→c+v(z,w)（u,v全linear）により

`f2→f2−h_b u−h_c v+X_b(v)−X_c(u)`。

ここでQ^{bc}=h_b b+h_c c+f2。
各canonical anchorで二次式の二係数を任意に消せる。

| anchor | coboundaryの全image | 残すf2 |
|---|---|---|
| A | span(z²,zw) | γw² |
| B | span(z²,w²) | γzw |
| S3 | span(zw,w²) | γz² |
| S1 | span(z²,w²) | γzw |

例えばAではimageは`(v_w−u_z)z²−2u_wzw`、Bでは`−u_z z²−v_w w²`。
S3では`2v_z zw+(v_w−u_z)w²`、S1では`−v_z z²−u_w w²`。
これらのshearsはP0をhigher common項だけ変え、全W sourceを保つ。
必要ならW→W+f(z,w)E_yで十八係数gaugeへ戻す。

actual diagonal`(a,b,c)→(λ^−2 a,λ^−1 b,λ^−1 c)`とtensor scalar λを組み合わせると、
P0及び全canonical anchorを保ったままγをγ/λへ送る。
従ってγ零又は非零の二whole orbitで尽くされ、γ=0,1だけを調べればよい。

## 7. 八fixed Qに対する全linear moduleと次元target

Wの全十八係数は

`W_b=z(l0 b+l1 c)+w(m0 b+m1 c)+r0z²+r1zw+r2w²`、

`W_c=z(l2 b+l3 c)+w(m2 b+m3 c)+r3z²+r4zw+r5w²`、

`W_z=W0_z+B_a zw+C_a w²`、`W_w=W0_w+D_a z²+E_a zw`。

W0はA/Bでbz∂w、S3でb(4z∂z+2w∂w)、S1で2bz∂z。
全Jacobiは(2)の140×18 affine linear systemであり、八fixed Qについてその全kernel、particular solution、
full original350との同値を確認する有限jobが十分。

Qの全orbit dimensionsは、anchor部分＋central coboundary二＋nonzero central amplitude一を戻し

| Q型 | δ=0 | δ=1 | fixed-symbol base B |
|---|---:|---:|---:|
| A | 5 | 6 | 5 |
| B | 6 | 7 | 5 |
| S3 | 5 | 6 | 6 |
| S1 | 5 | 6 | 6 |

である。各orbit上の全W fiberはactual normalizerでisomorphicなので、fixed-Q affine fiber dimensionをkとすると
full higher dimensionはorbit dimension+k以下。
plane六、tensor projectivization一を戻したimage upperはB+orbit+k+5。
従って必要なk upperはAで9/8、Bで8/7、S3/S1で8/7。
Wにactual shearの方向も全て含めているので、それを再び加算又は減算しない。

このrank targetを満たすか、actualDに固定factorがあるかを有限jobで確認して初めてG閉鎖とする。
足りない場合にはそのwhole affine familyとirreducibilityを保存し、新componentの可能性も保持する。
全scalar source、Π0零m四又は全Gの完了にはまだ拡張しない。
