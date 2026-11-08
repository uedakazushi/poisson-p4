# R2L m二constant-double最後の全p源：Ore全包含と唯一critical fiber

2026-11-08、PRIMARY G/E、producer `/root/quartic_review`。
rootの手係数導出を、既34源の完全核への別 exact substitution と全域幾何で確認する。
**新しい原operator、profile、minor、jet又はGBは計算していない。**
以下は全pのsourceを一つの必要critical fiberへ送る証明であり、そのcritical全支持の閉鎖は別担当の義務。

## 1. 元の全sourceとglobal係数関数

残sourceはp!=4上で

    {a,c}=−a², {b,c}=2ab,
    {c,u}=(p−3)au, {b,v}=au,
    {c,v}=bu+pav,
    other brackets=0,
    actualG2=(p−4)a³u².

全P1は固定忠実basis三十のkerJ1、全P2は原九係数を保持する。
P2をhomogeneous kernel三次元へ制限せず、affine particularを許す。
凍結34reportのgeneric及び全allowed wall証明書で、全p!=4のJ1はrank19、kerは十一dim。
p=−2では全u因子による空性が別にexact確認済み。
原行列/完全核/minorの別作者全再構成はrootが担当する。

全generic K1 basisで、c依存はPbcのc倍normal一次だけである。
B(u,v)=Pbc,c、Euler field B∂bを用いると

    P'=P+E∧(B∂b)=Q(a,b,u,v)+∂c∧W(a,b,u,v).

Wは元全sourceとの相関を持ったまま保持する。
ab成分のP2三係数は全arbitrary F2(u,v)、他の六P2係数はWへ入る。
全projective Jacobiのc係数/定数係数はstrict[Q,Q]=0及びE_y∧[Q,W]=0を与える。
元PがstrictPoissonであるという追加仮定は置かない。

T=raw Pav[u²]=Qav[u²]、m=raw Pcu[v²]を**global linear coefficient functionals**とする。
これらをgeneric basisの自由座標名へ依存させない。
全Qの形は

    Qab=a(Lu+Mv)+F2,
    Qav=T u², Qau=Quv=0,
    Qbu=Uu²+Vuv,
    Qbv=au+Ru²+Suv+W2 v².                         (1)

原Wとの混同を避け、最後のscalar coefficientをW2と書く。
原bad u-rowはm*v²∂cだけである。

## 2. 十一完全核へのexact確認

rootが手で得た全linear identitiesは

    V=−2(7p+2)m/[3(p−4)(p+2)],
    M=10m/[3(p−4)],
    W2=−(23p−2)m/[3(p−4)(p+2)],
    S−2U−L=U−3(p−4)T,
    2W2−2V−M=2V.                                 (2)

有限job `plane_double_line_m2_last_nil_ore_functionals.json` のgate後、
凍結reportの全十一generic basisを別sparse/Fraction実装で(1)/(2)へ代入した。
**0.219秒 / 97,192KiB、全PASS。**
Pab/Pacのc非依存、Euler除去後の全Q shape、元bad Ca/Cb零とm、及び全九P2の像も同時に確認した。
原J1/J2や新しいparameter点を再計算していない。
verifierはpublic JSON内のexact Python sourceとして保存し、凍結共有runnerを変更せず再生できる。

## 3. p=5等のgeneric-basis polesを落とさない全域論

J1(p)のentriesはregular polynomial in p。
全p!=4でrank19なので、K=kerJ1はB=A1\{4}上のrank十一vector bundleである。
Euler除去、m/T及び(1)の各係数は、このbundle上のregular linear mapsである。
shape条件と(2)の分母を(p−4)(p+2)で払った関係はpolynomial linear functionals。
generic full K上で零なので、locally free Kの局所frameと基底Bのdense generic pointを使えば全fiberで零。
従ってp=5を含む全chosen-minor/核basis denominator壁、二つのp²−10p+27 rootsも失わない。
generic basisに現れるp−5を反転して全sourceを捨てたのではない。

p=−2では(2)を割らず、既全u因子空性を用いる。
p=4は元actualG2零でnormal位数が変わる別義務。
よってp!=4,−2の全actualsourceに(1)/(2)が成り立つ。
この全域論は元34rank証明書の正しさへ依存し、それらの独立算術reviewとは区別する。

## 4. 二strictQ係数だけによる全sourceの三分

(1)からstrictcyclicを元定義で展開すると

    Jabu=T V u³,
    Jabv=T(S−2U−L)u³+T(2W2−2V−M)u²v.

F2の全係数はQuv=Qau=0及びb非依存から寄与しない。
(2)を入れると全複素解で

    TV=0, T(U−3(p−4)T)=0.                         (3)

T=0ならQは∂b∧V0、V0はb非依存。
module E_y∧[Q,W]=0へ∂bをwedgeしてactualF_b=0、PfQ=0からF_c=0。
全actualFは(a,u,v)だけに依存し、既 essential-three八成分へ全点被覆される。
W全係数とF2全係数を保持するconstant-factor補題の適用である。

T!=0かつp!=−2/7なら、(2)/(3)からm=0。
元全u-rowがu倍で、全九P2はnormal columnを持たないため、affine particularもactualFのu因子を除けない。
指定既約源は空。
従って真の必要残余は

    p=−2/7, T*m!=0                                (4)

の一fixed P0全fiberだけである。
T/mの零閉補集合は全pで処理済み。

## 5. critical全fiberの正確な入力とunit正規化

p0=−2/7はgeneric chosen-minor/complete-kernel denominator壁ではない。
従って既generic十一核をexactにspecializeすればcomplete criticalK1を得る。
新しいM1 profileを計算する必要は無い。
全P2九はそのまま戻す。
(2)/(3)の値は

    V=0, M=−7m/9, W2=−7m/18,
    U=−90T/7, S=L−180T/7,
    PfQ=−90T²u4/7.

この式だけでcritical全支持を空又は既成分へ認識したとは主張しない。
元full350とactualPfによる全source義務を残す。

normalizationはgeneric_exhaustionが提案し、本稿でも別に手確認した実PGLである。
push diag(A,B,C,U0,V0)、tensor scalar A/C とするとP0を保つ条件は
C=B²/A、V0=BU0/A。
全fiber上で

    m'=A³m/(B²U0), T'=A²T/(BU0).

A=1、B=m/T、C=m²/T²、U0=T²/m、V0=Tを取ればm'=T'=1。
全係数は(4)上のnonzero unitsで逆も明示的。
全残sourceを同時にtransformするので、一つのhigher点を選んだsliceではない。
T=0/m=0は前節へ戻し、失わない。

## 6. 正確な未了事項

これで全pの必要sourceは空、既essential三、又は(4)の**一つの全normalized fiber**へ送られる。
critical支持の最弱rowspace/actual-factor証明は別担当へ渡した。
同fiberの空性、G上界又はE membershipが独立に閉じるまで、R2L constantdouble全閉鎖を宣言しない。
他のR/m、HF八fixed源の個別E、全分類の完成は別義務である。
