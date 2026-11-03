# 全射影CA/CB交叉の被約成分：N27と新しい階数二のR19

2026年11月3日。ここではS61に制限せず、全射影Poisson空間の既知被約成分CA,CBの交叉を扱う。
全Poisson空間の新しい既約成分を加える結果ではない。

## 定理

全射影閉包の交叉の被約台は厳密に

    (CA∩CB)_red = N27 ∪ R19

である。右辺は相異なる二つの既約成分で、次元は27と19。
Nは既存の一般階数四の共通成分である。
Rは一般階数二の構造だけからなり、S61へ入らない。
従ってS61の全被約交叉が唯一Nだったことと両立する。
指定Q_(1)は階数四なので、Qを通る交叉成分は引き続き唯一Nである。

全交叉の大域的primary decompositionや、Rの一般点でのscheme厚みはまだ決めていない。
Qでの正確な218関係・主平方零nilradicalの局所定理は別に保持する。

## 1. Rの実際の八係数族

座標(h,x,y,z,w)でq0,…,q7を動かし、

    f=−4q4*x²,
    H=3h²q7z+2h q3 xz+h q6z²+q0 w x²+q1 x²z
                       +(q2/2)xz²+2q4 w xz+(q5/3)z³,
    P=f∂xy+∂y∧(−H_w∂z+H_z∂w)

とする。Rはこの族のPGL5飽和の射影閉包である。
これはCAのg=0境界の実際の線形部分族である。
同時に、CBのb=0境界を置換(1,3,2,0,4)で移した線形空間にも含まれる。
置換は旧座標添字iを新添字perm(i)へ送る規約である。
全包含は同じcanonical150係数空間で検査した。

八係数空間は8次元で、二つの表示旗の相対位置のGL5軌道は12次元。
その安定化群は両線形空間と交叉を保つので、付随する射影入射の次元は12+7=19。
したがってdim R≤19。
八パラメータの微分と全GL5作用のcanonical像は、指定QQ点でaffine階数20を持つ。
保存した20次小行列式は3623878656000である。
射影化で一つを除くとdim R≥19となり、厳密に19次元である。

## 2. 全階数二交叉を有限の66場合で覆う

[CAのproper旗模型](AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63.md)は、
全階数二境界をA_g(dim27)とA_f(dim25)で覆う。
[CBのproper模型](AMBIENT_CB_GLOBAL_FLAG_MODEL_20261103_v63.md)は、
全階数二境界をB_a(dim26)とB_b(dim27)で覆う。
[平方認識](AMBIENT_GLOBAL_BOUNDARY_SQUARE_RECOGNITION_20261103_v63.md)によりA_f⊂Nなので、
残るのはA_gとB_a,B_bとの交叉だけである。

各固定旗の空間は全てcanonical係数の線形部分空間である。
旗はF(1,3,4;5)、そのparabolicのWeyl部分群はz,wを交換するS2である。
二つの旗の対角GL5軌道は、標準的な相対位置の分類により

    S2 \ S5 / S2

の33個のdouble cosetで尽くされる。
B_aでは本来hlineを忘れられるが、全旗を用いる余分なP1は被覆を失わない。
したがって33置換×二つのCB端点、合計66の場合を検査すればよい。

各置換で、CAの19列と置換したCBの19列の150係数を並べ、
その差行列のQQ核から全線形交叉を求めた。
有限個の点の標本ではなく、各相対位置の**全交叉ベクトル空間**である。
全対の旗はこの33軌道に属し、proper模型は全境界点を持ち上げるので、これで全点を覆う。

結果は次の通りである。

- 64場合では全交叉がCA側のf∈hWを満たし、同じCA旗でNへ入る。
- 置換(3,1,2,0,4)、B_bの場合は七次元交叉で、CB側のfは常にx_B²の倍数。
  平方認識によりこの全族もNへ入る。
- 残る置換(1,3,2,0,4)、B_bの場合だけが第1節のRを与える。

各相対位置の族の飽和閉包は既約であり、CA/CBが閉かつ群不変なので、その閉包も実交叉に含まれる。
したがって全階数二交叉はN∪Rに入る。
一般階数四の全交叉は、基点を動かしてS61の全被約交叉定理を適用すればNに入る。
これで全射影交叉の台を覆う。

## 3. RがNに含まれないことを示す具体点

q4=q5=q7=1、他を零とする。生の構造は

    P=∂y∧V,
    V=4x²∂x−2xz∂z+(3h²+z²+2xw)∂w

となる。Vの四変数(h,x,z,w)でのdivergenceは8x。
canonicalな二ベクトルは∂y∧Cで、5Cは

    Ctilde=(-8hx, 12x², -18xz, 15h²+5z²+2xw)

というdivergence零の四変数二次ベクトル場である。
Eul∧v∧P_can=0の全175係数の行列は階数4、核は厳密にC∂y。
したがってこの点の内在的な定ベクトル因子は一意である。

Nの全階数二境界は、二つのproper閉像N_g=B_b,squareとN_f=A_f=B_a,squareで尽くされる。
次の二つの内在的な必要条件により、上の点はどちらにも入らない。

### N_fからの排除

A_fの生の構造はdivergence零なので既にcanonicalであり、二つの独立な定covector h,xを
常にCasimirとして持つ。GL5による移動でもこの条件は保たれる。
一方、上のcanonical Pの定covector contraction行列はQQ上で階数5。
従って非零な定covector Casimirは一つもなく、N_fへ入らない。

### N_gからの排除

N_gの任意の点は、一つの定ベクトル因子vと、四変数の生の場

    Y=−ell*(alpha*ell+beta*k)∂k+Hamiltonian

を持つ。ell,kは独立な一次式。
div Y=−beta*ellなので、四変数のtracefree代表C=Y+(beta*ell/5)Eulは

    C(ell)=lambda*ell²,
    C(k)=−alpha*ell²−4lambda*ell*k

を満たす。ここでは係数alphaを取り直しlambda=beta/5とした。
一意な定ベクトル因子が固定されると、これは商C5/Cv上で線形共役と倍率に不変な必要条件である。
canonical化がEulerの曖昧さを除くため、異なる生の表示を比較しても条件は変わらない。

上のCtildeについてell=A h+B x+C z+D wとする。
lambda=0ならCtildeの四成分の独立性からell=0となる。
lambda≠0ではCtilde(ell)=lambda*ell²のw²係数からD=0、次にh²からA=0、z²からC=0。
従って非零ellはxの倍数だけで、lambda*ell=12xとなる。
第二の必要条件はCtilde(k)+48xkがx²の倍数という条件に帰着する。
しかしk=A h+B x+C z+D wと書けば左辺は正確に

    40A hx+60B x²+30C xz+D*(15h²+5z²+50xw)

である。これがx²の倍数ならA=C=D=0で、kもxの倍数となり独立性に矛盾する。
従ってN_gにも入らない。

これでRの実点がNの全proper閉像の外にあることを証明した。

## 4. 被約成分の尽くしと限界

全交叉の台はN∪Rで覆われ、両方とも実交叉の既約閉部分である。
RはNに含まれず、次元27のNも次元19のRに含まれない。
よってこれらが全被約交叉の正確な二つの既約成分となる。

この19次元成分は全Poisson空間ではCA/CBに含まれる真閉部分なので、全成分数の下界57は変わらない。
全Poisson空間の尽くし、Rの一般scheme厚み、全交叉の大域primary decompositionは別の課題である。
Qでの局所交叉の正確なnilradicalとassociated primesは、階数四開の内部の既存定理として保持する。

## 再現

- `code/ambient_rank_two_flag_intersection_table_20261103_v63.py`：33相対旗×二端点の全QQ線形交叉。
- `code/ambient_global_intersection_completion_20261103_v63.py`：例外の全包含、19次元minor、非N点の内在的排除。
- `reports/AMBIENT_GLOBAL_INTERSECTION_COMPLETION_20261103_v63/`：全66場合、元の交叉基底、全入力hash。
- [独立な幾何監査](AMBIENT_GLOBAL_INTERSECTION_GEOMETRY_REVIEW_20261103_v63.md)。

全有限計算はQQ上で行い、標本の有限体ランクから尽くしを推論していない。
