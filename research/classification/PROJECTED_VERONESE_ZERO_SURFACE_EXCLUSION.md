# 滑らかなprojected Veroneseをrank零曲面に持つ枝の全排除

2026-11-07。主課題 `E_reduced_nonnormal`。
複素P4の滑らかなprojected Veronese曲面S上で零となるPoisson二ベクトルは零だけである。
特にSはrank4 Poisson構造のrank零曲面にはならず、既約非正規五次因子のこの具体的な枝を排除する。
他の非ACM曲面、特異射影、特異・非被約rank零scheme、既約非正規五次全体は排除していない。
本証明の独立作者reviewは別途行う。新成分を数える結果ではない。

## 1. 対象の曲面とPGL軌道

Sを

    [a:b:c] ↦ [a²−b² : b²−c² : ab : ac : bc]

の像とする。これはv2(P2)⊂P5をrank3の対称行列に対応する中心から射影したもの。
Veroneseのsecant varietyはrank≤2対称行列のdeterminant hypersurfaceなので、
中心はsecantとtangentの外にあり、射影は閉埋め込みである。
同じ型の任意の滑らかな射影Veroneseは、O_P2(2)の六切断から五次元線形系を選ぶことで得られる。
非secant中心はrank3対称行列全体であり、GL3の合同変換で一軌道となる。
従って全てPGL5を除いて上のSに一致する。
これは「全rank零曲面がVeroneseである」という分類仮定ではない。

このSはO_S(1)=O_P2(2)を持つが、ambient一次切断は五個、h0(O_S(1))は六個なので
h1(I_S(1))=1であり非ACMである。ordinary-double認識定理の自動適用は用いない。

## 2. 全150係数、全280式と完全商五次元

五斉次座標をx0,…,x4、二次同次二ベクトルをΠとする。
全150係数は、十個の添字対i<jと五変数二次単項式十五個の直積で並べる。
Euler列の外積とH1(T_P4)=0により

    H0(Λ²T_P4) = {全150係数Π}/{Euler∧全linear場},
    dim(Euler∧全linear場)=24.

各点で射影二ベクトルが零であることはx∧Π=0と同値。
上の五個の二次式をq0,…,q4と書き、全十個の三添字に対し

    q_i Π_jk(q) − q_j Π_ik(q) + q_k Π_ij(q) = 0          (1)

を課す。各式はa,b,cの六次同次式であり、その全28係数、合計280式を使う。
Sの像が閉埋め込みなので、pullbackでの零とS上でのsectionの零は同値である。
選んだnormal chart、標本点、因子条件、Jacobiの一部分で(1)を置き換えない。

QQ上の完全280×150行列Aについて

    rankA=121, dimkerA=29,
    全Euler24方向⊂kerA,
    dim H0(Λ²T_P4⊗I_S)=29−24=5.

121次minorの行列式は正確に1。
完全kernel29本とA・kernel=0、Euler24とcanonical五方向の和の階数29を保存したので、
rankの下界と上界が一致する。canonical化は

    Π ↦ Π − Euler∧div(Π)/5

で行い、全kernelのcanonical像の行空間RREFを五源の基底として固定した。
その全十括弧は `CLASSIFICATION_PROJECTED_VERONESE_ZERO_20261107.json`
の`complete_canonical_source`に、全280×150行列、全kernel、minorと共に保存してある。
パラメータをz0,…,z4とする。予想した三次元に制限したのではなく、完全空間は五次元である。

## 3. rank4はJacobiを課す前にも存在しない

全sourceに対する実anticanonical五次式は

    D = Σ_i (−1)^i x_i Pf(Πからi行・i列を除いた4×4行列)

である。これはEuler∧Π∧Πの係数で、Euler gaugeに依らない。
完全五変数zを保ったまま **D≡0** が成立する。
x0=1でEuler補正して得る実affine二ベクトルのPfaffianも恒等零となる。
従って、Jacobiをまだ課さない全線形空間の時点でprojective rank4は不可能である。

注意：4×4 Pfaffian五個自身は四次式であり、一般のhomogeneous liftでは
「x_iと共通五次式の積」にはならない。
最初のproducer中の共通比例を探す探索診断は正しい幾何恒等式ではなく、証明に使わない。
followupの `classification_projected_veronese_support.py` が上記の正しいEuler収縮を検証する。
元のproducer/certificateのimmutable bytesは変更していない。

## 4. 全Jacobi idealは原点支持：非零Poisson自体がない

さらに全350のprojective Jacobiを五源に代入する。
全てzの同次二次式であり、その350×15係数行列の階数は14。
その行空間は次の十四式の行空間とQQ上で一致する。

    r0=z0²−z4²,       r1=z0z1−z4²,    r2=z1²−4z4²,
    r3=z0z2,          r4=z1z2,        r5=z2²−z4²,
    r6=z0z3,          r7=z1z3,        r8=z2z3,
    r9=z3²−z4²,
    r10=z0z4,         r11=z1z4,       r12=z2z4,      r13=z3z4.       (2)

全350式から(2)を得る十四個の有理row functionalと、逆に全350式を(2)で表す係数を保存した。
従って単なるJacobiの部分集合でなく、**完全な固定SのJacobi ideal**が(2)である。
以下の多項式恒等式は、その根基をただちに定める。

    z0³ = z0*r0 + z4*r10,
    z1³ = z1*r2 + 4z4*r11,
    z2³ = z2*r5 + z4*r12,
    z3³ = z3*r9 + z4*r13,
    z4³ = z0*r10 − z4*r0.

全rが原点で零となることと合わせ、根基は(z0,z1,z2,z3,z4)。
従ってS上で零となる非零Poisson sectionは一つもない。
Euler方向だけのliftはprojective sectionとして零であるため、別のPoisson構造を作らない。
この支持排除はMora、形式積、flatness、有限jetの積分可能性、源のgeneric点の外挿を用いない。

## 5. 全軌道への結論と証明書

PGLはPoisson条件とsectionの零曲面を保つので、§1の唯一軌道性により全滑らかprojected Veroneseに
同じ排除を適用できる。これは一つの具体的な非ACM rank零曲面分岐の全排除である。
他の曲面の有限分類、非正規五次因子全体の排除又は全Poisson成分の尽くしを主張しない。

- 線形制限・完全五源・全Jacobi/Pf：`code/classification_projected_veronese_zero.py`
- 完全十四式・両方向functional・全cube恒等式・正しいD：`code/classification_projected_veronese_support.py`
- 証明書：`CLASSIFICATION_PROJECTED_VERONESE_ZERO_20261107.json`, `CLASSIFICATION_PROJECTED_VERONESE_SUPPORT_20261107.json`
- gate済み有限契約：`jobs/projected_veronese_zero.json`, `jobs/projected_veronese_support.json`

全線形源とJacobiの再構成4.58秒、支持・functional証明1.01秒。
独立作者はsurface substitution、全kernel/Euler商、全350式、十四式・根基、PGL軌道を別算術で確認する。
