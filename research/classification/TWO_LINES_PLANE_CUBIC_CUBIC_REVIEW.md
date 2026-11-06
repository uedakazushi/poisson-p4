# 平面二直線・三次曲線の残三源：固定切片による独立監査

次の三つの既約実源の飽和閉包は、それぞれ射影次元23の全 Poisson 成分である。
全源と固定切片の算術を、producer を import しない標準ライブラリ Fraction 実装で
独立に再構成し、すべて受理した。三成分相互および既存成分との相違、平面配置の
全尽くしは `TWO_LINES_PLANE_CUBIC_SUPPORT.md` の別の幾何証明に依存する。

| 源 | 全 D 階数 | 実源接像 F | ker D | 真の法線 |
|---|---:|---:|---:|---:|
| node_two_flexes | 100 | 48 | 50 | 2 |
| node_flex_secant | 100 | 48 | 50 | 2 |
| cusp_tangent_cusp_third_secant | 101 | 48 | 49 | 1 |

## 実際に独立再構成したもの

全六つの magnetic パラメータを保持する源について、strict Jacobi、実 Pfaffian、
全六微分方向を再構成した。各元点で全350×150 Jacobi 微分 D、PGL・scale・
鉛直方向、全源接像 F を再計算した。鉛直階数24、群等の像47、実源像48である。
全実方向と公開法線を D が消すこと、両者の独立性、公開された D の非零 minor
を有理数で確認した。従って表の階数は下界・上界の双方が確定している。

source report は最初の二源について
`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json`、最後について
`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_THIRD_SECANT_20261107.json`。
法線は既存 `NORMALS` / `TWO_NORMALS` report、三次証明書は
`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_CUBIC_SUPPORT_20261107.json` にある。

## 固定切片と三次初期式

公開された r 本の標準座標方向を U、n=102−r 本の法線を R とする。
D|U の公開 minor が可逆で、F⊕R=ker D であるため
\[
 \mathbb C^{150}=F\oplus S,\qquad S=U\oplus R,\quad\dim S=102.
\]
これは実際の線形補空間である。群作用による直積分解は使わない。

選択された r 本の Jacobi 式の陰関数定理により、U 座標は法線変数 t の形式級数となる。
その二次項を `T_ij` とする。公開 T_ij がすべて U に属することと
\[
 (D T_{ij}+Q_{ij})|_{\text{selected rows}}=0
\]
を独立に確認した。可逆 minor により、これが固定切片の一意な二次項である。
一般には全350行で右辺が零になる必要はなく、残る行が二次初期式を与える。
その全係数空間も独立に再構成した。

各公開 λ について全150列で λD=0、全 raw normal quadratic について λQ_ij=0
を確認した。従って、同じ λ を実際の IFT 残式へ適用すると二次以下が消え、
三次初項は正確に `λ B(R,T₂)` となる。公開された全三次係数ベクトルと全 λ の
pairing を独立に再計算した。固定切片の一意 T₂ だけを使うため、任意の ker D
補正に対する不変性は要らない。generic_structure もこの形式論理を独立に受理した。

## 被約支持に必要な最小証明

- **node_two_flexes**：全二次空間は
  `t0²−16t1², t0t1−4t1²`。その支持は `t0=4t1` だけである。
  三次 functional [3] のこの直線上の値は
  `230605056000000*t1³`。従って支持は原点である。
- **node_flex_secant**：二次式は `q=t0²+(9/4400)t1²`。
  三次 functional [1] を q で割った剰余は
  `−174254976000/11*t0*t1²`。
  `t1²q−t0(t0t1²)=(9/4400)t1⁴` より根基は極大イデアルである。
- **cusp_tangent_cusp_third_secant**：二次空間は零だが、三次初項は
  `840558307758750*t0³` なので原点だけを支持する。

以上の初期式は真の IFT 残式のイデアルに由来するため、102次元線形切片上の
局所支持は零次元である。元空間は48本の線形式でこの切片へ切られるので、
Krull の次元不等式から元局所次元は48以下。一方、実既約源の微分は階数48を持つ。
その飽和閉包は従って全空間成分であり、鉛直24と射影化1を除いて次元23となる。
局所唯一成分性、稠密な一軌道、正確な非被約局所環は主張しない。

## 保存と再現

- 独立実装：`code/classification_two_lines_plane_cubic_cubic_review.py`。
- 再利用した独立算術：`code/classification_two_lines_plane_cubic_component_review.py`。
- 契約：`research/classification/jobs/two_lines_plane_cubic_cubic_review.json`。
- 独立証明書：`reports/CLASSIFICATION_TWO_LINES_PLANE_CUBIC_CUBIC_REVIEW_20261107.json`。

全入力 hash、独立実装 bytes、minor、pairing、支持の証明を報告内に保存した。
実行時間は約9.05秒。既存出力を保護した一時 checkout で契約を gate し、実装を実行する。
既存報告のある場所では上書きを拒否する。この検証は全源成分性の算術と次元論を
確認するものであり、台帳の相違証明や全 Poisson 構造の分類完成を代替しない。
