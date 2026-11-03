# 正規五次因子の認識：必要方向、閉点同値、非被約基底の監査

2026年11月3日。v14 `sections/normal_quintics.tex` を必要方向から全scheme認識まで読み、
各次数と基底変更の論証を独立に確認した。以下の三つを区別する。

1. 一般階数4のprojective Poisson構造と**正規な実D**から、正準affine lift、Pfaffian–勾配、
   高さ3の最小Buchsbaum–Eisenbud分解を得る普遍的な必要方向。
2. 正規な固定五次Fの四つの閉点認識条件と、一意なprojective Poisson構造。
3. 実際の最大minor determinantal scheme上の、非被約基底を含む全scheme同型。

三つとも記載した仮定の下で有効であり、未解決の証明上のgapは見つからなかった。
保存manifestは `reports/NORMAL_QUINTIC_RECOGNITION_AUDIT_20261103_v62.json`、コードは
`code/normal_quintic_recognition_audit_20261103_v62.py`。
旧Fisher例の二つの独立計算は再実行せず、保存結果とsource hashを照合した。
これらの例を普遍命題の証拠として代用していない。

## 1. Projective Jacobiから正準affine Jacobiへ

S=C[x0,…,x4]、Euler E=Σxi∂i、divergence δを用いる。
任意の二次lift Prawから

    P=Praw−E∧δPraw/5

を取る。δ²=0、δ(E∧A)=5A−tr(A)EからδP=0。
垂直なE∧Aにδ=0を課すとAがEの倍となるため、そのbivectorは零。
従ってこれはprojective bivectorの一意な線形splitであり、任意のC基底変更と両立する。

このliftがaffine Poissonである点を、旧版の「基本対応」という引用だけで終わらせず確認した。
T=[P,P]/2は三次係数の3-vectorである。projective JacobiはE∧T=0、
BVのδ–Schouten恒等式とδP=0はδT=0を与える。
同次divergence恒等式

    δ(E∧T)=5T−E∧δT

からT=0となる。係数環を非被約C-algebraに替えても同じ多項式恒等式である。
従って追加のaffine Jacobi仮定は不要である。

α=Σbi(P)dxi=i_(P∧P/2)(dx0∧…∧dx4)とするとBVからdα=0。
αのEuler重みは5なのでCartanにより

    F=(Σxi bi)/5, dF=α, ∇F=b(P).

x0=1のprojective affine PfaffianはΣxi bi=5F。
従ってFはπ∧πが定める**実際の**退化切断を表し、被約支持から後で作った五次式ではない。
交代行列のPfaffian恒等式はP∇F=0を与える。

## 2. 正規性と高さ3の最小分解

実D=(F=0)が正規なら、Dは被約であり、projective空間で異なる超曲面成分が交わることから既約でもある。
三次元Dの特異集合の次元は高々1なので、斉次Jacobian ideal J_Fの高さは3以上。
他方、5×5交代行列のsubmaximal Pfaffian idealの高さは3以下。
従って高さは厳密に3である。

ここでBuchsbaum–Eisenbudを適用する環は**体上の正則多項式環S**であり、height=grade。
char0なので必要な2も可逆である。従って

    0→S(−10)→S(−6)⁵→P S(−4)⁵→(∇F)^T S→S/J_F→0

が完全。全非零成分が正次数なので最小であり、五偏微分はC上一次独立、低次数のsyzygyはない。
S/J_Fは高さ3のGorenstein環、深さ2で、J_Fは飽和。全Jacobian曲線Yは純一次元CMとなる。
Hilbert分子1−5t⁴+5t⁶−t¹⁰からdegY=40、pa(Y)=101を得る。
ここではC=Yredの既約性、滑らかさ、次数を仮定していない。
後の不活性成分・E6/E8帰着へ必要なのはこの方向である。

## 3. 閉点認識の逆方向と一意性

正規五次Fに対しK_F={P∈Alt5(S2):P∇F=0}と置く。
以下は同値である。

- Fが一般階数4のPoisson構造の実退化因子。
- htJ_F=3、S/J_FがGorenstein、五つの四次偏微分が一次独立。
- K_Fが有理関数体上でrank4の行列を含む。
- K_Fは直線で、その生成元のPfaffian vectorが非零。

第二条件からは高さ3 Gorenstein idealの自己双対な最小分解が
(1,5,5,1)、shift(0,4,c−4,c)となる。
Hilbert分子のt=1での一階零条件は−20+5(c−4)−c=0なのでc=10。
従ってBE行列の成分は二次で、そのPfaffian五本はJ_F,4の基底となる。
任意の基底変換をcongruenceで勾配基底へ直せる点も確認した。
p(CPC^T)=det(C)C^(−T)p(P)であり、欲しいM∈GL5に対して
d⁴=detM、C=d M^(−T)を選べばdet(C)C^(−T)=M。
この段階は複素閉点の議論なので第四根を取ってよい。

第三条件から直接戻る場合は、rank4のPの核が一次元なのでb(P)=h∇F。
正規Fの偏微分は共通因子を持たない。もし共通因子があればEuler式からFもその因子で割れ、
正規性の特異次元条件に反する。UFDでhの分母は定数となり、両辺次数4からhは非零定数。
従って同じ高さ3分解を使える。

P∇F=0を微分するとHessianの対称性から(δP)·∇F=0。
分解には一次係数のsyzygyがないのでδP=0。α_Pは∇Fの定数倍で閉じているためBVから[P,P]=0。
これで認識の逆方向が成立する。

一意性は単なる次元標本ではない。任意のR∈K_Fの各列は二次syzygyなので、完全性からR=PA、
Aは定数5×5行列。PA∇F=0とkerP=S∇F、次数比較からA∇F=c∇F。
五偏微分の一次独立性がA=cIを強制するためR=cPである。
従ってprojective Poisson構造は一意であり、K_Fを解けば復元できる。

## 4. 非被約基底でのnormal openの意味

全scheme命題で「正規な五次族」と言うときは、**幾何的fibresが正規である開集合**への写像を意味する。
非被約A上の全hypersurfaceそのものが正規環であるとは言っていない。
平坦な普遍projective hypersurfaceのgeometrically normal locusをU_norと取る。
この意味を保持すれば、Artin基底やnilpotentを含む基底も対象にできる。

M(F):Alt5(S2)→S6⁵は1050×150のcentrality行列。
**実際の全150次minor ideal**を零とし、ある149次minorが単元のchartを取り、
核生成元のPfaffian係数の少なくとも一つが単元となる開へ制限したschemeがZ_norである。
rank149という被約な点集合に置き換えてはいけない。

局所環A上で149次blockを消去するとMはdiag(I149,residual column)となる。
残る901成分は、それぞれ採用blockへ一行を加えた150次minorの単元倍。
従ってnilpotent成分も含めて残余列は**厳密に零**となる。
核は一つの係数を1に取れる自由直和因子A·Pで、P∇F=0はA上の真の恒等式となる。

## 5. Fibrewise BEから相対完全性へ：還元性を使わない補題

任意の局所環(A,m,κ)で、有限自由加群の有界complex C_n→…→C_0を取る。
κ上で正のhomological degreeのhomologyが零なら、このcomplexは正次数でsplit exact、H0は自由である。

証明は左端からの有限Gaussian消去である。
最左の非零C_n→C_(n−1)はκ上単射なので、全列rankのminorがAで単元となる。
恒等写像の可縮直和因子を取り除くと、complex条件から次写像はその成分で零。
この操作を繰り返すと次数0の自由加群だけが残る。
初めからcokernelの平坦性を仮定せず、Aのdomain性やreducednessも使わない。

b=b(P)のPfaffian恒等式により、相対BE列は完全性の前から真のcomplexである。
各斉次次数を固定すると、それは有限自由A加群の有界complex。
各residue fieldでは、幾何的正規性・非零Pfaffian・§3から高さ3となり、体上のBEが完全性を与える。
residue fieldが非代数閉なら代数閉包で確認し、忠実平坦性で降ろせる。
上記補題を**全ての斉次次数**へ適用すれば、相対BE列はA上で完全となる。
各graded cokernelもA-flatとなり、任意のbase change後も完全である。
非CMな一般Aへheight3のBE定理を直接適用しているわけではない。

## 6. Nilpotentを含む正確な復元とgluing

∇FはPの核の四次成分なので、相対完全性のdegree10部分から

    ∇F=c b(P), c∈A

となる。各residue fieldでc≠0なのでcは単元。
zero divisorを含むAの分数体で比例を論じる必要はない。
degree5のsyzygyには直前のS_A(−6)⁵が寄与しないため、§3と同じ微分からδP=0がA上で厳密に従う。
相対多項式微分ではcは定数なのでα_P=c^(−1)d_xFは閉じ、BVから[P,P]=0。
Euler式によりF_P=c^(−1)Fである。

逆にPoisson族では、正準liftは§1の線形splitで与えられる。
あるP係数が単元なのでM(F)P=0は150次minorをA上で厳密に零にする。
各residue fieldでrank149なので、149次単元minorのchartが族を覆う。
そこで復元する核直線は元のP直線に一致する。

Fの単元倍はkernel lineを変えず、Pの単元u倍はF_Pをu²倍にする。
従って両projective構成はglueし、相互に逆でbase changeと両立する。
§3の複素閉点の基底正規化と違い、相対復元ではcの平方根も第四根も取らない。
このため全scheme同型P_nor≅Z_norはnilpotentを落とさず成立する。

## 7. 今回の研究への適用と限界

新しい明示的πについて、元Jacobi、一般階数4、**実F**の正規性を別々に確認できれば、
§1–2と[正常Dの閉点分類監査](NORMAL_ACTUAL_EIGHT_SUPPORT_AUDIT_20261103_v62.md)を接続して、
πはPGL5と非零尺度を除きq5,1型と同定できる。
この同定は証明された存在型の認識であり、具体的なPGL5変換行列の算出とは別である。
例えば高さ3の実Jacobianがあれば、可約・非被約Fではその高さが高々2となること、
hypersurfaceのS2とR1から、正規性を判定できる。
支持の被約化の正規性だけではこの適用はできない。

本稿は相対認識の同型まで監査したが、任意基底上のPGL軌道分類、全Poisson schemeのreducedness、
全既約成分数、非正規または非被約の実Dの分類をそこから結論していない。
v14認識自体はv22の閉点分類と独立であり、determinantal schemeのnilpotent方向も保持する。
