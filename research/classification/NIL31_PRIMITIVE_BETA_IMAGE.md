# 新しい nil31 成分の primitive β 像は次元高々13

2026-11-07。課題 G。一般五点の尽くしを前提にせず、新しい25次元成分が
全16次元の primitive 一次底を優越するかを判定する。
本ノートの非優越証明は全源の実 fiber を使い、標本の微分階数を上界と読み替えない。

## 1. 全源と、一次 jet を変えない三次元族

入力は `reports/CLASSIFICATION_RANK2_NIL_31_COMPONENT_20261107.json` の
全10変数源、局所化 n7≠0 である。
モジュラー場を rank2 nilpotent Jordan 型 (3,1)、定数項零へ正規化した全該当枝であり、
その PGL 飽和閉包が25次元全空間成分であることは別の全周囲計算と独立確認による。

源の全六括弧のうち、n0,n1,n2,n3,n4,n8 は P12 の項

\[
 6\bigl(n_0+n_4x_3+n_8x_4+n_1x_3^2+n_2x_3x_4+n_3x_4^2\bigr)
 \partial_1\wedge\partial_2
 \tag{1}
\]

だけに現れる。他の五括弧はこれらに依存しない。
従って、全ての二変数二次以下多項式 F(x3,x4) を (1) に加えても、
**全 Jacobi 条件と実 Pfaffian は保持される**。
これは第一階の接方向だけではなく、全 F に対する実 Poisson 族である。

任意の symplectic 基点 p を取る。その座標を `(p1,p2,p3,p4)` とすると

\[
 (x_3-p_3)^2\partial_{12},\quad
 (x_3-p_3)(x_4-p_4)\partial_{12},\quad
 (x_4-p_4)^2\partial_{12}
 \tag{2}
\]

の三方向を任意に加えた affine 三次元族は、P(p) と p における全一次微分を厳密に固定する。
三方向は多項式 bivector として一次独立であり、Pfaffian も不変なので rank4 の開から出ない。
源を affine 座標変換で動かした任意の点でも、その変換で p を引き戻して (2) を作り、
押し戻せば同じ三次元の実族を得る。

## 2. 全域の次元上界

**命題。** この成分の固定 symplectic 基点と定数値 J における primitive β 像の次元は高々13。
従ってこの新成分は W=A16 を優越しない。

**証明。** モジュラー場を固定した全源は10次元。
Jordan 型 (3,1) の中心化群は次元6、核の平行移動は次元2なので、
場 `Z=Nx` の Aff4 軌道は `20−6−2=12` 次元である。
したがって全 affine 飽和実支持の次元は高々22。

固定基点で定数値が非退化となる開において、評価写像は GL4 同変であり、
非退化な定数 bivector の六次元空間への GL4 作用は推移的である。
よって固定値 J の支持の次元は高々16。
第1節の三次元実族は評価値も一次 jet も固定するため、この固定 J の支持内で同じ β を持つ。
全ての実点で β の fiber は次元少なくとも3であり、像の次元は高々 `16−3=13` となる。

有限 perspective による primitive 正規化は固定 J の一次部だけから一意に決まり、
同じ一次 jet を持つ (2) の全族へ同じ可逆変換を適用する。
従って三次元族を primitive gauge へ戻す際にも次元は失われない。
任意の重複超平面を無限遠へ送る perspective も β を保つので、
`C0_S61_SLICE_SCALAR_AUDIT_20261106.md` 第1節の全橋を通じて PGL 飽和後の主張へ戻る。
像閉包と成分の境界については、この β 写像が定義される開での閉包の連続性を用いる。□

この証明は一般五点の正しさにも、下記の一点での微分計算にも依存しない。
全源・全低次項・基点の移動・primitive gauge を保持した fiber 下界である。
独立担当および統合担当が、実三次元 fiber、場の軌道12、評価の次元減少6と
perspective の適用範囲を読解確認した。

## 3. 一点の exact 微分による別の下界

全10変数を1とした既存の source witness と、affine 基点 p=(0,0,1,1) を使った。
整数化された源の定数 Pfaffian は −48 で、非退化である。
全10源方向、GL4 の16方向、平行移動4方向、scale1方向を取り、
定数六係数と一次24係数への微分を作る。

定数値を C=P(p)、中心座標を y=x−p とすると、perspective の四方向は

\[
 G_a(y)=-2(a^Ty)C+y\wedge(Ca),\qquad a\in\mathbb C^4.
\]

これを定数零として加えた30×35行列の階数は23、
定数評価の階数は6、perspective の階数は4であった。
従って primitive β 微分の階数は `23−6−4=13`。
定数評価の kernel を先に取り、その一次部を perspective 四方向で割る別の計算順でも同じ値を得た。
可変な Darboux 正規化の微分を無視した計算ではなく、評価 kernel が定数値を固定する接空間そのものである。

`code/classification_component_beta_differential.py` と
`jobs/nil31_primitive_beta_differential.json` による有限 QQ 計算は約0.43秒で完了した。
`reports/CLASSIFICATION_NIL31_PRIMITIVE_BETA_DIFFERENTIAL_20261107.json` は
全行列、非零 minor、全 kernel、入力 hash と producer bytes を保持する。
この数値の literal な独立再構成は未実施である。
第2節の非優越にこの下界の独立再計算を追加前提としない。
