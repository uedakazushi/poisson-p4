# モジュラー階数二冪零層の primitive 一次像は非優越

2026-11-07。G と E_multiple_hyperplane に共通する全支持の像上界。
一般 rank4 Poisson 構造で、重複超平面を無限遠へ移した affine modular 線形部が
rank2 nilpotent となる全点を扱う。全点の既知成分帰属は主張せず、
**primitive β の像次元が高々15**という最弱の結論を証明する。

## 1. 全 affine モジュラー型と既存の有限全源

線形 Jordan 型は(3,1)または(2,2)。定数部を image で割り、中心化群で整理すると五型となる。
以下の fixed source 次元は、全90係数・全 Jacobi を保持した全源の上界。

| 型と定数類 | fixed source 上界 | affine modular orbit 次元 | β像上界 |
|---|---:|---:|---:|
| (3,1), image 内 | 10 | 12 | 13 |
| (3,1), 長い鎖 top | 5 | 14 | 13 |
| (3,1), 短い鎖 top | 6 | 13 | 13 |
| (2,2), image 内 | 14（§2） | 10 | 15 |
| (2,2), image 外 | 7 | 12 | 13 |

全源の証拠は `MODULAR_NIL31_COMPONENT.md`（n7開と二つの補集合）、
`MODULAR_NIL31_LONG_SUPPORT.md`、`MODULAR_NIL31_SHORT_SUPPORT.md`、
`MODULAR_NIL22_TOP_SUPPORT.md`。
固定した線形部の GL4 orbit 次元は、中心化群次元6,8よりそれぞれ10,8。
affine 平行移動が加える image は二次元。
(3,1) の cokernel 上で長い top の orbit は二次元、短い top は一次元、
(2,2) の非零 cokernel orbit は二次元。これで表の orbit 次元を得る。

## 2. (2,2)、定数部が image 内の全 source は高々14次元

モジュラー場を `Z=p1∂q1+p2∂q2` とする。
全90係数の線形 reduction は、全 P を

\[
 P=F(q,p)\partial_{q_1}\wedge\partial_{q_2}
   +\sum_{i,j=1}^2M_{ij}(p)\partial_{q_i}\wedge\partial_{p_j}, \tag{1}
\]

へ制限する。M は次数高々二で p のみの式、F は q に関して affine-linear。
これは `CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json` 第四 chart の
全19変数表示であり、Ppp を追加仮定で零にしたのではない。
div P=Z は F の二つの q 係数を M から一意に定める。
従って M を定めた後の自由度は、p の任意の次数二以下の式 F0(p) の六個だけである。

`D=det M` は Pfaffian の非零定数倍。
最高 Q の rank≤2 により deg D≤3。
第一添字 sharp の modular 恒等式 `P♯dD=DZ` から

\[
 M\nabla D=-Dp,\qquad \operatorname{adj}(M)p=-\nabla D.       \tag{2}
\]

後式は D≠0 の開から多項式恒等式として全体へ延びる。
p=0 で ∇D=0 なので `D=D0+D2+D3`、添字は斉次次数。

\[
 H=\operatorname{Hess}(D_2)+\tfrac12\operatorname{Hess}(D_3)
\]

とすれば Hp=∇D。二変数の一次 Koszul syzygy により、(2)の全 M は

\[
 M=-\operatorname{adj}H+p\otimes V(p),\qquad \deg V\le1.      \tag{3}
\]

実際、`(adj M+H)p=0` の各行は `(p2,−p1)` の affine-linear 倍であり、
2×2 adjugate が線形であることから (3) を得る。
行列式を取ると必要条件

\[
 D=\det H-V\cdot\nabla D                                   \tag{4}
\]

となる。この条件の全解を分類する必要はなく、次元上界だけを用いる。
パラメータは D2 の3係数、D3 の4係数、V0 の2係数、V1 の4係数。
D0 は (4)の定数項 `D0=det Hess(D2)` で定まる。

### Hess D2 が rank2

次数一の式は `V0·∇D2=(det H)1` であり、V0 の二係数を一意に定める。
次数二の式における線形写像 `V1↦V1(D2)` は rank3。
従って cubic 条件を課す前でも M の自由度は

\[
 3+4+2+4-2-3=8
\]

以下である。

### Hess D2 が rank1

rank1 二次形式の全体は二次元。座標を取り D2=p1²/2 と書く。
`D3=a p1³+b p1²p2+c p1p2²+d p2³` とすれば、次数一の式は
`V01=c,d=0` を独立に課す。
次数二の写像 `V1↦V1(D2)` は rank2 なので、必要条件だけで

\[
 2+4+2+4-2-2=8
\]

以下となる。残る次数二・三条件を無視することは上界を弱めるだけである。

### Hess D2=0

D0=D2=0 なので rank4 には D3≠0 が必要。
最高式は `V1(D3)=−D3`。
非零 binary cubic の三根が相異なる場合、二重根の場合、三重根の場合の
stratum 次元はそれぞれ4,3,2、`T(D3)=0` となる線形場 T の空間次元はそれぞれ0,1,2。
これは代表 `p1p2(p1−p2),p1²p2,p1³` の係数比較で直接確認できる。
Euler/3 の負号が特解なので、V1 の解空間次元はこの stabilizer 次元に等しい。
V0 の二自由度を加えても、各 stratum の M 次元は高々6。

全ての場合に M の次元は8以下。F0 の六自由度を戻して fixed source は高々14次元。
この証明は一般点の Hessian が可逆という仮定を置かず、その補集合を全て扱った。

## 3. symplectic 値の正規化と三つの実 Taylor fiber

各 fixed source を全 affine 変換で飽和した空間の次元は、表の source 上界と modular orbit
次元の和以下。symplectic 基点を原点へ移した開では、GL4 が非退化定数 bivector 全体に
推移的に作用するため、定数値を J に固定すると六次元を引ける。

(3,1),image 内の全三源では Pxy に任意の二変数二次式 F0(z,w) を加えられる。
(2,2),image 内では (1)の F0(p) が任意である。
各 affine 基点 a で、その二変数に関する

\[
 (p_1-a_1)^2,\qquad(p_1-a_1)(p_2-a_2),\qquad(p_2-a_2)^2
\]

の三方向は、全 Poisson 条件と全 modular 型を保ち、定数項と一次項を全く変えない。
従って fixed-J locus の primitive β 写像の実 fiber は少なくとも三次元。
これは三つの接ベクトルだけでなく、全係数で存在する affine 三平面である。
群で移した点についても、元の座標に戻した基点を中心として同じ三平面を取れる。

これより image 内の二型では

\[
 10+12-6-3=13,\qquad 14+10-6-3=15.
\]

残る三型は Taylor fiber を使わなくても、それぞれ `5+14−6=13`、
`6+13−6=13`、`7+12−6=13`。
重複超平面を無限遠へ取る perspective は primitive β を保つことが
`REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` §3 で証明されているので、
超平面を動かす四次元を β 像の上界へ足してはいけない。

従って、モジュラー線形部が rank2 nilpotent である全 repeated-H 支持の
primitive β 像は高々15次元、W16 へ優越しない。

## 4. 独立確認と未解決部分

別担当が (2)の符号、二変数 Koszul 表示、det 恒等式、Hess rank1 の二つの独立条件、
binary cubic stabilizer 次元、三つの全 Taylor fiber と β 像の上界を独立に確認した。
本証明に Gröbner 計算や標本点の数値階数を使わない。
全八源などの実 Jacobi 検算は元の report に保持する。

この非優越は、残る rank2 nilpotent 各源が既知成分へ入るという結論ではない。
E の全支持被覆は該当 source ノートに残す。一般 β の未知優越成分を除くためには、
他の modular 型と被約非正規因子をさらに処理する必要がある。
