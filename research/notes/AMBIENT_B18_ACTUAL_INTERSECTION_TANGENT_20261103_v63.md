# B18原点のactual相対接空間：CA=CBの39平面

2026-11-03。個別CA/CBの全formal型を相対的に同一視せず、同じ物理150座標でcotangent embeddingを構成する。
対象は指定P0を持つopen B18 orbitとそのG移動である。

## 元のcommon R29とμ(I)=10

以前のactual quadratic-image reductionは、**同じ29 target関数**をbase Rとする有限free rank2 normalizerをCA/CB双方について与える。
各actual imageはR⊕I t、t²=nu、t(0)=0、I⊂m_Rである。
新しいabstract formal modelが構成するR座標を、ここへ同じ座標として持ち込まない。

maximal idealの二乗を取ると、odd部はm_R I t、even部はm_R²（nu I²⊂m_R³）である。
従ってactual cotangentはm_R/m_R²とI/m_R Iの直和。
個別actual ringの既受理のintrinsic edim39とdim R=29から、**元のRに対してもμ(I)=10**となる。
実seven linear classes、WU/WVのtwo classes、actual U³のone classは、既存の全source監査でI/mI内に独立に構成されている。
よってこの7+2+1は全basisであり、未計算の追加cotangent方向はない。

## even部をsource kernelの二乗から正しく復元する

全GL作用と九source fiberから28個の独立なeven一次方向Jを、rawテンソルから独立に再構成した。
ramification kernel pに沿う**純粋な実source curve**を全有限5×5行列でpushforwardし、全150の一次係数零と二次係数Wを確認する。
Rはactual target関数なので、このcurveでR(p)=R0+O(p²)。trace-zero t(p)=O(p)かつI⊂m_Rだから、odd tIはO(p³)である。
従ってWはgenuineな**even Rの一次embedding**だけを与える。
任意のsource二次jetを新しい接方向に加えた議論ではない。

共通29 target rowsをCとするとM=(J|W)[C,:]は単元であり、

    E=(J|W) M⁻¹,  E[C,:]=Id29

がeven cotangent embeddingである。
再構成した29次minorはCAが9216/25、CBが−576/25。

## odd部の実10 basis

全actual odd first-R-jetからseven actual generator rows

    18,19,28,29,60,67,72

を選び、全150係数のcombination行列Lを再計算する。
zero-full-firstjetの全quadratic coefficientからrows46,95の実generatorを取り、全150係数Qを正規化する。
actual monic critical U-axisの全cubicからLとQを差し引くと、唯一のrow63が残り、係数はCAが12/125、CBが−4/243。
それを正規化した列をC3とする。

B=(L|Q|C3)は**I/mIの全10 classesの物理embedding**である。
七つのactual monic elimination、全Neumann補正、mIを殺すcritical評価は既存の独立監査を明記して再使用する。
単にsource p²やp³のspanをcotangentと呼び替えない。

## actual CA/CB/Dの接空間

T_CA=(E_CA|B_CA)、T_CB=(E_CB|B_CB)を同じphysical150 rowsで作る。
29 common rowsと上の10odd rows、計39個の**実target座標**で正規化すると、両150×39行列は全係数で同じになる。
従って全CA、全CBの接空間はともに39次元であり、和も39、交わりも同じ39平面である。
actual intersection Dはideal sumなので

    T_D=T_CA∩T_CB=T_CA=T_CB, dim T_D=39。

q106=1のscalar chartを取ると、射影のembedding dimensionは38となる。

## Krull次元、multiplicity、depth上限

全reduced intersectionはN27∪R19で、P0は両方の成分に入る。
coneでは唯一のtopdim成分N28と、incomparable minimal component R20がある。
従ってactual Dのcone dimは28、射影dim27。
Nに沿うactual Dのgeneric lengthは1、P0でのactual reduced Nのmultiplicityは2と既に独立に受理されている。
multiplicity加法性にはtopdim成分だけが寄与するので、**e(D,P0)=2**。
lower-dimensional Rやembedded成分をこの和へ足さない。

minimal prime R20もassociated primeなので、depth-localization不等式からdepth D≤20（射影≤19）。
これは正確なdepthやfull primary decompositionの計算ではない。

結果はcone dim28/edim39/e2/depth≤20、射影27/38/e2/depth≤19。
**T_Nが37次元であることだけから、D_redの接空間やnilradicalのlinear visibilityは推論しない。**
N∪Rの接空間はN一成分の接空間とは異なり得る。
generic B18のactual relative full ideal、full D原点環、正確なdepth・nilpotence・primary typeは未決として残す。

再現:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_actual_intersection_tangent_20261103_v63.py --output <new-output-directory>

全even/odd embedding、正規化した共通39平面、厳密minorとsource Wを保存する。120秒・2.2GB、higher jet/GBなし。
