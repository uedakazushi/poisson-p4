# normal CI22・rank二全支持の排除：独立確認

2026-11-07。producer `plane_collision_review`、reviewer `root`。
**全normal CI22・全零divisor Rの不存在証明をPASSとして受理する。**
non-normal CI22又はprojected scrollへnormal Hartogsを拡張しない。

## 1. global sourceと全normal pencilの幾何

normal Gorenstein surfaceの滑らかな部分のrank二bivectorは、adjunction/Hartogsで
OS(1)へ延びる。Koszul ACMから全五sectionである。
JU=ι(dq2)ι(dq1)Uを用い、Eとの外積がell_U(x) times全tangent three-vectorとなることを、
contraction公式dq_i(E)=2q_iから別に計算した。したがって全restrictionをJUでlift出来る。
残差のEuler kernelをH1(OS)=0でglobal WSへliftし、linearnormalityでambient linear Wへ。
全sourceはJU+q1A+q2Bであり、S上だけのgaugeをglobalにあると仮定していない。
同じliftでnormal contractionsはIに入り、局所qi/Xj²の補正はI²となる。
従ってstrict{q1,q2}∈I²という全conormal方程式を使える。

regular pencilのskew centralizer/normalizerの既受理全Jordan表を読み直した。
dimN≤1は全sourceに共通quadric因子を与える。
dimN≥2でnormalな残りはsem221、mixed[(21),(11)]、nil41の三型。
geom kernel三以上はrank二quadricを含みintegral性に反し、長さ二以上の二Jordan blockの
common eigenkernelはq1-isotropicな二空間なのでsingular line、normalでない。
single-eigenvalue nil partitionでも同じ判定を全partitionに適用した。
三/二/一固有値のpencilPGL2正常化に連続modulusは残らない。

all-member singular pencilでは、M2はplaneを含みnon-integral、M1+regular2はsingular line。
common radical二以上もbase零次元からreducible。normalとして残るのはcommon radical一の
smooth ellipticquartic coneだけであり、四つの相異lambdaを全て保持する。

## 2. 三有限源を元350式へ戻した独立算術

`classification_normal_ci22_rank_two_review.py` はproducer/CAS/legacyをimportしない。
外積contractionから全五JUを別構成し、全q-gradient Casimir、全JU線形和のJacobi/Pf零を確認。
constant skew全十係数とq1/q2二補助係数の15×12 conormal行列を全て再構成し、
rank十、完全augmented kernel二、表示R/T二generatorの全張りを確認した。
各全九parameter源の全150係数、全350Jacobi、全126actualPfを独立Fractionで再計算した。

sem221のac,bc,a*u0,a*u1,d*u2,d*u3,d(a+b)、
mixedのac,bc,a*u0,a*u2,d*u3,d*u4,d(a+b)−2b*u2、
nil41のc²と2ad−5bcを、全original350行への各exact functionalとして確認。
det(ad−bc)非零なら前二型はc0/a,d非零、b=−a、JUはそれぞれU4又はU1のみへ。
この全三parameter縮約のJacobi零とactual固定x4/x1因子も直接確認した。
nil41では被約点上c0/ad0、det非零openに点は無い。
det零の全補集合は初めから共通quadric因子で処理し、捨てていない。

実行前gate後0.080秒で全算術PASS。
保存行列が疎coefficient辞書であるため最初のdecoderを訂正したが、
数学の入力・scope・producer原bytesは一切変更していない。

## 3. ellipticconeの全modulusを別に検算

q1=Σxi²、q2=Σlambda_i xi²（base四座標）、vertex vとする。
全conormal源はJU+∂v∧(q1V+q2W)。JUのvertex-free coefficient cが零なら全sourceが∂vを共有しrank二以下。
c非零ではJU全体のstrict compatibilityにより、元Jacobiは
Ebase∧(q1 Lie_V Jbase+q2 Lie_W Jbase)=0を必要とする。
missing dのcomponentをxi-axis（i≠d）で読むと、Lie derivativeで生き残る項は
Jbaseのcomplementary product xi*xdのxd微分だけであり、係数は非零scalarを除き
`(lambda_d−lambda_i)(V_d+lambda_i W_d)`。
相異なる二つのiでV_d=W_d=0、全dでV=W0となる。
残るJUはrank二以下で、全相異lambda・全vertex係数でrank四は無い。

以上は全Rのreducedness・integrality、normal-coneのnode、vertical不存在を仮定しない。
数学的に得た因子/空性だけを受理し、全E帰属、全schemeの厚み又は総数上界を主張しない。
