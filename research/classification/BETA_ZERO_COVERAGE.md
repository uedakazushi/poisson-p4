# 零一次部分の全被約台を既知成分で覆う

2026-11-06 / v64-classification-first。課題E_zero。

## 主張と入力定理

元134式の零一次部分を $X_0=\{\beta=0\}\subset S61$ とする。
**既存のQ45被約22成分定理と、既知成分の同定を入力とすれば、
$X_{0,\mathrm{red}}$ は既知52個の階数四成分の閉包の和に含まれる。**
これは特殊底一点 $\beta=0$ の全台の被覆であり、全特殊底や全S61の尽くしではない。
nil指数・一般厚み・全軌道・交叉の計算は要らない。

入力Q45定理は `legacy/v07/quadratic_classification.tex`, Theorem Q22。
その内容は、$[J,Q]=[Q,Q]=0$ を満たす斉次二次テンソルの被約台が下の22成分の和になるというもの。
この定理の全場合分けを今回再証明したわけではない。
ただし以下の**22成分から既知の全Poisson成分への包含**は、名称の類似でなく、
実際の族・稠密軌道・多項式変形から接続する。
元134式を実際に制限した証拠は
`reports/CLASSIFICATION_E0_BETA_ZERO_20261106.json`（二次70行、零64行、45変数）。

SNCの先行定理はMPSのTheorem1.4：SNC log symplectic locusと交わる全Poisson成分は40個。
一次資料の適用範囲を確認した記録は `research/REFERENCES.md`、
`research/notes/MPS_CATALOGUE_20261103_v62_REVIEW.md`、
`legacy/v23/proofs/literature/SS_geometry_and_MPS_identification.txt` のTheorem1.4参照。
40個の既知族との同定は既存の図式対応と局所Torelliの証明を入力として保つ。

## 22件の包含表

| Q45成分 | 包含先 | 全族を覆う理由 |
|---|---|---|
| G1 | MPS40の和 | 一般A,D,EでSNCへの実PGL極限、残りは閉包 |
| G2 | MPS40の和 | 一般A,B,EでSNCへの実PGL極限、残りは閉包 |
| T1 | MPS40の和 | 非零QスカラーとSp4飽和の全体にSNC極限 |
| T2 | MPS40の和 | 同上 |
| T3 | MPS40の和 | 同上、別の正重みを使用 |
| T4 | MPS40の和 | 同上 |
| E1 | MPS40の和 | 同上 |
| E2 | MPS40の和 | 同上 |
| E3 | MPS40の和 | 同上 |
| E4 | MPS40の和 | 同上 |
| NB | R4_TRIPLE_HYPERPLANE_QUADRIC | Q45代表の二倍が既知の二次部分、Sp4軌道閉包 |
| M0 | R4_CA | 全9パラメータの三次ポテンシャルを明示 |
| M1 | R4_CB | 全9パラメータの三次ポテンシャルを明示 |
| M | R4_MINIMAL_M | 同じ代表と、10次元Sp4×Gm稠密軌道 |
| nil22開枝 | R4_H_QUARTIC_ESSENTIAL5 | 同じ代表と、10次元Sp4×Gm稠密軌道 |
| Cplus | R4_CPLUS_EXT | 10次元稠密軌道の代表が拡張族c=-1,beta=0に入る |
| Cminus | R4_CMINUS_MAG | 三磁気係数全体の座標変更、kappa非零開から閉包 |
| Dpm | R4_DPM | Q(sqrt3)上の同じ代表と、10次元Sp4×Gm稠密軌道 |
| U1 | MPS40の和 | 全磁気項のconfluenceからG1型、上のG1包含を適用 |
| A0 | R4_A0 | 同じ代表と、10次元Sp4×Gm稠密軌道 |
| AB | T01 | 全磁気項を含む、五独立一次式上のlogcanonical族への変形 |
| L0 | R4_ZERO_LOG_5H | 定義に使う二パラメータ族そのものが既知族 |

T1..T4は**Q45の名称**であり、MPS内部IDのT01..T40とは異なる。
特定MPS成分への一意な所属を決める必要はない。有限個の既知閉成分の和への包含で十分である。

## 閉包への移行の共通補題

$\Phi(Q)=[h^2J+Q]$ はQ45から全射影Poisson空間へのmorphismである。
非零な定数部分が固定されるので、$Q=0$ でも定義される。
$\Phi$ はSp4-equivariantで、$Q$ の非零定数倍はhの再尺度化と
テンソル全体の定数倍によって射影同値になる。
従って、閉じたPGL5不変集合Cに $\Phi(Q_*)$ が入れば

\[
 \Phi\bigl(\overline{(\mathrm{Sp}_4\times\mathbb G_m)Q_*}\bigr)\subset C.
\]

同様に、パラメータ空間の稠密開の像がCに入れば、その全像の閉包もCに入る。
有限個の既知成分の和も閉じている。この初等的な閉包操作を各行で使う。

## 非冪零10成分：実際のSNC極限

座標を(h,x,y,z,w)とし、係数 $x^\alpha\partial_i\wedge\partial_j$ の
pullback重みを $\alpha\cdot\omega-\omega_i-\omega_j$ とする。
Q45の正常形は旧Q22定理の式をそのまま用いる。
次の重みでは $h^2J+Q$ の全項が非負重みで、零重み項はlogcanonical項だけになる。

| 型 | omega | 極限Pfaffianをxyzwで割った係数kappa |
|---|---|---|
| G1 | (3,2,2,1,0) | AD |
| G2 | (2,0,0,1,1) | AB |
| T1 | (3,2,2,1,0) | 3 |
| T2 | (3,2,2,1,0) | 3/2 |
| T3 | (3,0,2,3,0) | 3/2 |
| T4 | (3,2,2,1,0) | -3 |
| E1 | (3,2,2,1,0) | 1 |
| E2 | (3,2,2,1,0) | 1/2 |
| E3 | (3,2,2,1,0) | 1 |
| E4 | (3,2,2,1,0) | 1/2 |

G1のAD非零、G2のAB非零はそれぞれの源の稠密開条件。
T/EのQをlambda倍すればkappaはlambda二乗倍され、非零lambdaで消えない。
各全項・重み・極限の六係数・Pfaffianの正確な恒等式を新しい証明書に保存した。
最初に提案した全10件共通の重みはG2のzw項で不適切だったので、G2/T3の二件を上の通り訂正した。
不適切な共通重みを証明として用いていない。

重要なのは、射影退化因子が**単なるrank4ではなく**

\[
 \tfrac12 E_5\wedge Q_{\log}\wedge Q_{\log}
       =\kappa\,hxyzw\,\partial_{hxyzw}
\]

になることである。これは五つの独立な座標超平面の被約SNC因子である。
非零epsilonでは重み族は元の $h^2J+Q$ の実際のPGL5軌道にある。
PGL5は連結なので各全空間既約成分を保つ。元の点を含む任意の成分はその軌道閉包、
従ってこのSNC極限を含み、MPSの40成分の一つである。
よって元の点は既知MPS40の和に入る。Sp4、非零スカラー、パラメータの閉包で10行全体を覆う。
SNC因子を持つ条件の一般の開性を仮定しておらず、特殊極限から元の因子もSNCだと推論していない。

## M0とM1：全9パラメータの実包含

Q22の最小冪零正常形の9パラメータをt1..t9とし、c=0,1/4に対して

\[
\begin{split}
 H_c={}&t_4x^2z+\tfrac12t_5xz^2+(t_6+c)xzw+\tfrac13t_7z^3
       +\tfrac12t_8z^2w+t_9zw^2\\
     &-t_1x^2w-\tfrac12t_2xw^2-\tfrac13t_3w^3
\end{split}
\]

と置く。c=0では $f=h^2+x^2/2$ として正確に $P_A(f,H_0)$、
c=1/4では $f=h^2+x^2/4$ として正確に $P_B(f,H_{1/4})$ に一致する。
両方の全六係数差の恒等的零を今回検証した。
源の9方向を標本に固定せず、元の全族をCA/CBに含めている。
既知CA/CBの定義は `legacy/v06/ambient_families.tex`, (PA),(PB), TheoremAB。

## 五つの稠密軌道による包含

M、nil22開枝、A0、Cplus、DpmのQ45成分はいずれも既存のQ22定理で10次元既約である。
今回、Sp4の十Lie方向とQのスカラー方向からなる全60×11係数行列を正確に作り、
次の10角非零minorを得た。

| 型 | 正確なminor |
|---|---|
| nil22開枝 | 16 |
| A0 | 243/128 |
| M | 6879707136 |
| Cplus | -8257536 |
| Dpm | (-9499+3641sqrt(3))/30233088 |

各代表は指定Q45成分内にあり、その群軌道もそこに含まれる。
次元10の下界と既存の成分次元10が一致するので、その軌道閉包はQ45成分全体である。
単一代表の所属を全族の所属に取り違える穴は、ここで閉じる。
Dpmの演算体はQQ(sqrt3)で、rを超越変数や浮動小数としていない。
最後のminorが非零なのはsqrt3が有理数でないことからも分かる。

nil22開枝、M、A0、Dpmの $h^2J+Q_*$ はそれぞれ
`code/additional_rank4_families_20261103_v62.py` の同じ名前の代表そのものである。
CplusではQ*=16倍した旧Cplus正常形であり、$h^2J+Q_*$ は同コードのc=-1,beta=0。
従って共通補題により、各行のQ45成分全体が指定の既知閉成分に入る。
必要な代表・全60×11係数・選択行列式を証明書に保存した。

## 残る五行の全族包含

**NB。** Q22のNB代表の二倍は、既知triple-hyperplane-quadric代表の二次部分に一致する。
Q22はNB成分をそのSp4軌道閉包と記述する。
非零倍をhの再尺度化で吸収する共通補題から、NBの全斉次化像が指定成分に入る。
代表の全六式はQ22と `legacy/v05/triple_hyperplane.tex`, (brackets) にある。

**Cminus。** 旧v8 `sections/cminus_component.tex` は三磁気係数a,b,cを残して包含を証明する。
$X=x+z,Z=x-z,Y=(y+w)/2,W=(y-w)/2$ のsymplectic変更とテンソル8倍、hの再尺度化で、
源は $P^{xy}=h^2+x^2,P^{yz}=-xz,P^{zw}=h^2+x^2-2z^2$ と
$P^{yw}=-xw-Ax^2-Bxz-Cz^2$ になる。
ここで $A=a+b+c,B=2a-2c,C=a-b+c$。
剪断 $y'=y+Az,w'=w+Ax$ の後、$\kappa=C+2A\ne0$ でfiberを再尺度化すると
既知CminusMag族の $\beta=B/\kappa,\gamma=0$ になる。
$\kappa\ne0$ は三係数空間の稠密開であり、その補集合も既知成分の閉性で覆われる。
この全三係数族のSp4飽和がQ45のCminus成分を稠密に覆うことをQ22から用いる。

**AB。** v8 `sections/ab_toric.tex` は、正規化したコアの**任意の二次磁気項**B(h,u,v)を含む族に、
二つのパラメータr,sを持つ実際の多項式変形を構成する。
一般rs非零では $h,L_+,L_-,M_+,M_-$ が行列式4klrsの独立な一次座標で、
括弧がlogcanonicalそのものになる。r=s=0が元の全磁気項族。
源の一般lambda非零コアはsymplecticな対角変更で同じコアへ正規化でき、
三つの磁気係数は可逆な二変数変更でBの二次項へ移る。
Bのu二乗/v二乗はfiber剪断で除去でき、v8の四係数族で覆われる。
従って全AB成分の斉次化像は標準トーリック成分T01に含まれる。
元の一つのp=(1,2,3)代表だけの包含を使っているのではない。

**U1。** v8 `sections/u1_t16.tex` のconfluenceは、一般epsilon非零で
$h^2J+N(s_1,0,0,s_2,c)$、すなわちG1型に入る実多項式族である。
特殊化の磁気商座標は $(1,-\eta/2)$。
全二次磁気項B(h,x,z)へのfiber剪断の像は4次元で、商座標は
$(b_{hz}-2b_{hx},b_{xz}-4b_{xx})$ の二つ。
この第一座標が非零の開は、共通fiber再尺度化とetaで尽くされる。
その開は全六磁気係数空間で稠密なので、閉じた既知MPS40の和は全Bを含む。
beta=0のB(x,z)はこの六変数族の特殊部分として入る。
beta=0だけを見ると第一商座標は恒等的に零なので、そこでの開性を誤って主張しない。
最後にQ22のU1正常形飽和の閉包へ移す。
旧表の名称T16をMPSの内部番号と無条件同一視せず、新しいG1包含だけで十分である。

**L0。** Q22の零モジュラー成分L0は、
`legacy/v06/zero_log_component.tex`, (zeroQ)のQ(a,b)のSp4飽和閉包として定義される。
その全 $h^2J+Q(a,b)$ は、既知R4_ZERO_LOG_5Hの定義に使う族そのもの。
従って全L0像の包含が直ちに従う。

## 結合と証拠の範囲

入力Q22定理によりX0の各閉点のQは、表の22閉成分の少なくとも一つに入る。
各行の全像を覆ったので、その点は既知52成分の閉包の和に入る。
これが被約台について求めた包含である。
逆向きの全52成分のbeta0ファイバーの分類は要求せず、成分を追加して数えてもいない。

新しい正確な恒等式とminorは
`reports/CLASSIFICATION_BETA_ZERO_FULL_COVERAGE_20261106.json`。
実装は `code/classification_beta_zero.py`。
残る実族包含の既存係数証拠は
`legacy/v08/data/reproduced/{ab_toric_certificate,u1_t16_certificate,u1_t16_independent_audit,cplus_component_certificate,independent_cminus_audit}.json`。
本文の閉包・密性の議論と、それぞれの証明書の射程を併用する。
Q22全場合分けとMPSの外部証明は既存の依存であり、今回の小計算がそれらを再証明したとはしない。

ジョブ契約は `jobs/beta_zero_crosswalk.json`。
E0が作業中にdoneへ変わった際、一回の継続計算で旧task-idのgateがinactiveを返した。
根の担当者によるbeta0継続の明示割当は存在し、数学入力は変えていない。
現在は専用のprimary E_zeroへ再束縛し、以後の計算前にそのgateを通す。
