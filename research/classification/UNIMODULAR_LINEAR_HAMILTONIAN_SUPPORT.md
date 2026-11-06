# Unimodular最高部に一次Casimirがある場合の全一次像排除

2026-11-07。G。全二次 affine Poisson 構造 P のうち `div P=0`、
Pf(P) が非零定数で、最高二次部 Q に非零の一次 Casimir がある全点を扱う。
**結論は、primitive 一次評価の階数が高々三となり、W16を優越しないこと。**
特殊 cubic の正則列や一般固有値を仮定せず、全低次延長を保持する。

## 1. 零点のない affine Hamiltonian への正確な還元

Q♯dℓ=0、ℓ≠0 を一次式とし、`V=P♯dℓ=Ax+c` と置く。
Pf(P) が非零定数だから P は全 affine 空間で可逆であり、V は零点を持たない。
また ℓ(V)=0、div V=0、[V,P]=0。
最後は V が P の Hamiltonian 場であることから従う。

従って `c∉im A`、`im A+Cc⊂kerℓ` であり、rank A≤2、tr A=0。
これは最高部だけの正規化ではなく、元の全 P から得られる必要条件である。

A=0 なら c≠0、[c,P]=0 から L(c)=0。
A が非冪零なら rankとtraceの条件により固有値は(a,−a,0,0)、a≠0で半単純である。
平行移動で定数部の im A 成分を除き、`V=S+c0`、c0≠0、[S,c0]=0 とできる。
次数二以下の二ベクトル空間において ad S は半単純、ad c0 は冪零で可換。
`(ad S+ad c0)P=0` から両作用が P を消す。
特に [c0,P]=0、従って L(c0)=0。
いずれも primitive L の評価が非単射なので、真の閉部分集合に入る。

以下では A を非零冪零とする。
P の全体 scalar 倍と ℓ の逆 scalar 倍は V を変えず、評価の階数も保つので利用できる。

## 2. 全六つの nilpotent 対

座標を (x,y,z,w) とする。cをim Aで割り、Jordan中心化群で整理すると、
V と ℓ の全対は次の六型となる。

| A型・定数類 | V | ℓ |
|---|---|---|
| (3,1)、長い鎖top | y∂x+z∂y+∂z | w |
| (3,1)、短い鎖top | y∂x+z∂y+∂w | z |
| (2,2)、非零定数類 | y∂x+∂y+w∂z | w |
| (2,1,1)、長い鎖top | y∂x+∂y | w |
| (2,1,1)、短い鎖top・鎖のcovector | y∂x+∂z | y |
| (2,1,1)、短い鎖top・他のcovector | y∂x+∂z | w |

最初の三型では ker Aᵀ と c の直交条件により ℓ の直線が一意に定まる。
(2,1,1)の長いtopでは、残る二つの長さ一の鎖をGL2で動かしてℓ=wにできる。
短いtopでは許されるℓはspan(dy,dw)にあり、im Aᵀ の直線dyか、その補集合かの二軌道。
後者は `w↦w+t y` とscalarでdwへ取れる。
c∈im A は V に零点が生じるため最初から存在せず、除算で落とした枝ではない。

## 3. 全45係数の線形問題

P♯dℓ=V は六括弧のうち三つを完全に固定する。
残る三括弧には二次以下の各十五係数、計45係数を全て残す。
div P=0 と Pf(P) の正次数係数零は、これらの未知数について affine 線形である。
後者が線形なのは、Pfaffianの各積が必ず固定括弧を一つ含むためである。

四次元で div P=0 のとき、BV恒等式は Jacobi と d Pf(P)=0 の同値を与える。
従ってこの線形系を解くことは全 Jacobi の支持を解くことそのものである。
下記の源でも全四 Jacobi 多項式を独立に直接展開し、恒等零を確認した。
各源の自由パラメータは全て複素数であり、表示したPfaffianの非零条件だけを課す。

### (3,1)、長いtop

固定括弧は `(Pxw,Pyw,Pzw)=(−y,−z,−1)`。
全残括弧は四変数源

\[
 Pxy=n_0+n_1z+n_2wz-n_3y+n_3z^2,\quad
 Pxz=n_1+n_2w+n_3z,\quad Pyz=n_3.
\]

Pf=−n0、rank四条件はn0≠0。全一次項はxに依存しない。

### (3,1)、短いtopは空

固定括弧は `(Pxz,Pyz,Pzw)=(−y,−z,1)`。
自由括弧をA=Pxy,B=Pxw,C=Pywと置くと
Pf=A+yC−zB=k、div P=0から

\[
 V(C)=1,\qquad V(B)=C,
 \quad V=y\partial_x+z\partial_y+\partial_w
\]

が従う。第一式の定数項は coeff_w C=1。
第二式のw係数は `2 coeff_(w²) B=1`。
しかしPfのzw²係数は `−coeff_(w²) B` であり、定数Pfの条件に反する。
従って次数二以下の全解は存在しない。

### (2,2)、非零定数類

固定括弧は `(Pxw,Pyw,Pzw)=(−y,−1,−w)`。
全六変数源は

\[
\begin{aligned}
 Pxy&=n_1+n_2w-n_5wy+n_5z,\\
 Pxz&=n_0+n_1w+n_2w^2+n_3y+n_4wy+n_5wz-y^2,\\
 Pyz&=n_3+n_4w+n_5w^2-y.
\end{aligned}
\]

Pf=n0、rank四はn0≠0。やはり全一次項はxに依存しない。

### (2,1,1)、長いtop

固定括弧は `(Pxw,Pyw,Pzw)=(−y,−1,0)`。
全十変数源は

\[
\begin{aligned}
 Pxy&=n_0+n_1z+n_2w+n_3z^2+n_4wz+n_5w^2+n_8x-n_8y^2/2,\\
 Pxz&=n_6+n_7y+n_8yz+n_9wy,\\
 Pyz&=n_7+n_8z+n_9w.
\end{aligned}
\]

Pf=n6、rank四はn6≠0。
一次係数行列でPxzの行はPxwの行に比例し、PywとPzwの行は零。
従って全評価階数は高々三。

### (2,1,1)、短いtop、ℓ=y

固定括弧は `(Pxy,Pyz,Pyw)=(−y,1,0)`。
全十変数源は

\[
\begin{aligned}
 Pxz&=n_0+n_1y+n_2w+n_3y^2+n_4wy+n_5w^2+n_9x-n_9yz+z,\\
 Pxw&=n_6+n_7y+n_8y^2+n_9wy,\\
 Pzw&=n_7+n_8y+n_9w.
\end{aligned}
\]

Pf=n6、rank四はn6≠0。
一次係数行列でPxwの行はPxyの行に比例し、PyzとPywの行は零なので階数は高々三。

### (2,1,1)、短いtop、ℓ=w

固定括弧は `(Pxw,Pyw,Pzw)=(−y,0,−1)`。
全十変数源は

\[
\begin{aligned}
 Pxy&=n_0-n_7y-n_8y^2-n_9wy,\\
 Pxz&=n_1+n_2y+n_3w+n_4y^2+n_5wy+n_6w^2+n_7z-n_8x+2n_8yz+n_9wz,\\
 Pyz&=n_7+n_8y+n_9w.
\end{aligned}
\]

Pf=−n0、rank四はn0≠0。
一次係数行列でPxyの行はPxwの行に比例し、PywとPzwの行は零なので階数は高々三。

## 4. 全primitive一次像への帰結と残る最高部

ここでdiv L=0なので、固定Darboux frame後のLはそのままprimitive一次項である。
全六型、半単純型、定数型で線形評価が非単射となり、この条件はaffine基点変更後も同じ全型表で成立する。
`REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` 第14節にはinjectiveなW16の元があるため、
全像は真のdeterminantal閉集合に入り、一般βを優越しない。

最高二形式が `γ=dℓ∧df3` の場合、対応するQはℓを線形Casimirとして持つ。
従って本結論は全 rational(1,3) 最高部に適用でき、特異三次式、正則列の失敗、低次項の全境界を含む。
最高Qに線形Casimirがない場合、特に二つの二次式によるJacobiansや残る非安定型は別に残る。

公開再現は `code/classification_unimodular_hamiltonian.py`、
全源と全45係数への代入は `CLASSIFICATION_UNIMODULAR_AFFINE_HAMILTONIAN_20261107.json`。
契約 `jobs/unimodular_affine_hamiltonian_sources.json`、全源の実測2.310秒。
独立検証 `CLASSIFICATION_UNIMODULAR_AFFINE_HAMILTONIAN_REVIEW_20261107.json` は
別の Fraction 実装で全45係数行列を再構成した。
六型の係数行列階数は順に41,39,39,35,35,35で、第二型の拡大階数40から空性が従う。
全五源の自由度・全 Jacobi/Pfaffian/divergence/Hamiltonian 恒等式・全15minor零を確認した（0.068秒）。
統合担当は全六対の被覆、semisimple/constant補集合、短鎖空性の手証明も独立に読解し、
本ノートの全線形Casimir枝のG非優越を受理した。
