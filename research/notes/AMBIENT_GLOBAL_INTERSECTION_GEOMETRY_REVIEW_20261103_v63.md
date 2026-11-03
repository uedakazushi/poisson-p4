# 全域 CA∩CB の二成分を尽くす旗入射の論理監査

2026年11月3日。これは全射影P4 Poisson空間の既知の被約成分CA,CBについての、
**global containmentとcomponent countingの概念的監査**である。
66ケースの各係数行列・minorの独立再計算は別の監査に委ねる。
本ノートは有限表の数値的正確性を標本から推定しない。

## 入力と結論の関係

次の既に証明された幾何入力を用いる。

* CAとCBの全proper旗模型。全generic Poisson rank≤2境界は
  CAではAg27∪Af25、CBではBa26∪Bb27である。
* rank4開では全被約交叉がN27に等しい。
* Nの全proper19次元固定旗模型と、全rank≤2境界の認識
  `N_g0=Bb,square`、`N_f0=Af=Ba,square`。

さらに有限表
`code/ambient_rank_two_flag_intersection_table_20261103_v63.py`
の次の **完全な係数結果** を入力とする。

1. 33個のdouble cosets×二endpoint、計66ケースを全て計算している。
2. 64ケースでは全linear intersectionのCA側fがhで割れる。
3. 一例外ではCB側Fが恒等的に平方である。
4. 最後の例外のlinear spaceは8次元、paired-flag orbitは12次元で、
   実際のcanonical150係数像のQQ微分階数は20である。
5. 最後の例外には以下の具体的な点が入り、別の正確な証明によりNに入らない。

これらの有限係数入力が成立すれば、以下の論理により無条件に

    (CA∩CB)_red=N27 ∪ R19

が従い、右辺が正確な二つの既約成分となる。
従って残る検証は有限係数入力の独立監査であり、別の未知の境界族を仮定していない。

## 1. どのrank2点も有限表へ入る

任意のP∈CA∩CBでgeneric Poisson rank≤2とする。
CAの全proper模型から、PはAgまたはAfに入る。
Afは既にNに含まれるので、新たに必要なのはP∈Agの場合だけである。

CBの全proper模型から、PはBaまたはBbへ持ち上がる。
Baの最小の旗底は(L,W)だけだが、有限表では **hlineも選んだfull flag底**を使っている。
この冗長な旗選択はP1の全射properな追加であり、覆いを壊さない。
したがって両endpointとも同じfull flag多様体G/Pを使ってよい。

同時に二つの旗を動かすG作用のorbitは、parabolic double cosetsで分類される。
ここでPのWeyl部分群はz,wの交換S2であり、P\S5/Pは33個である。
よって全ての旗対は有限表の33個のpermutation代表のいずれかへ移せる。

固定した旗対のsourceは、二つの **実canonical線形空間**の射影交叉そのものである。
raw代表をEuler項を無視して同一視したのではない。
その全linear intersectionを計算した表だから、これで全閉点を尽くす。

## 2. 各ケースの次元上界と既約性

permutationをw、共通安定化群をH=P∩wPw^−1とする。
paired flag orbitはG/Hで、その次元は表のdである。
交叉linear spaceをLw、次元をkとすると、該当sourceは

    G×^H P(Lw)

で、次元はd+k−1である。
像の閉包はG×P(Lw)という既約空間の像の閉包なので既約である。
この上界はtableの一点での微分階数から推定したものではない。

最後の例外はd=12,k=8なので上界19である。
実canonical150係数のaffine微分階数20はprojective像次元の下界19を与える。
したがってこの閉像Rは正確に19次元である。
Ag上の全source点はrank≤2であり、その条件は閉だからRの全点もrank≤2である。

## 3. 65ケースがNに入る理由

64ケースでCAのfがhで割れるなら、対応する全linear spaceは **同じCA旗で** Nの19空間に入る。
NはG不変な閉集合なので、その全G-saturationの閉包もNに入る。
これは一般点だけの認識ではない。

残る一例外は、元のCA旗ではfがhで割れないが、CB側Fが全係数で平方となる。
この場合は別旗を許す正確な平方認識 `Bb,square=N_g0` を用いる。
閉proper像の等式なので、F=0や重なった線を含む全linear familyがNに入る。
そのため単にgeneric座標変換を特殊点で割って使う必要はない。

## 4. 最後の例外はNに吸収されない

具体点を

    P=∂y∧V,
    V=4x²∂x−2xz∂z+(3h²+z²+2xw)∂w

と取る。
その唯一のintrinsic constant factorはL=C∂yであり、canonical代表は∂y∧C、
商U=C5/Lの上で

    C(h)=−8hx/5,
    C(x)=12x²/5,
    C(z)=−18xz/5,
    C(w)=3h²+z²+2xw/5

である。

Afの全raw表示はdivergence zeroなので既にcanonicalで、二つの独立な定数covector h,xを零化する。
canonical化はGL5同変であるため、この性質は全Afで保たれる。
一方このPの定数covector contraction mapは階数5で、零化するcovectorが一つもない。
従ってP∉Afである。

Ng0の場合、unique factor Lでquotient Uを同定すると、表示

    Y=−ell*(alpha*ell+beta*k)*∂k+X_H

のcanonical fieldは `C=Y−div(Y)*E_U/5=Y+(beta*ell/5)E_U`となる。
したがって独立な二つの線形形式ell,kが存在し

    C(ell)=lambda*ell²,
    C(k)=−alpha*ell²−4lambda*ell*k

を満たす必要がある。これは同じLを保つGL5作用の商GL4に関して不変な条件である。

ell=a h+b x+c z+d wと置く。
lambda=0はCの四成分の線形独立性に反する。
lambda≠0なら、C(ell)−lambda*ell²のw²係数からd=0、h²係数からa=0、z²係数からc=0となる。
従ってellはxの倍数であり、lambda*ell=12x/5である。
しかし任意のk=a h+b x+c z+d wについて

    C(k)+(48/5)x*k
      =8a hx+12b x²+6c xz+d(3h²+z²+10xw)

がx²の倍数になるにはa=c=d=0が必要で、kもxの倍数となる。
これは独立性に反する。したがってP∉Ng0である。

Nの全rank≤2境界はAf∪Ng0だから、P∉Nが従う。
この点のexact coefficient certificateは
`code/ambient_rank2_exception_not_N_20261103_v63.py` と
`reports/AMBIENT_RANK2_EXCEPTION_NOT_N_20261103_v63/report.json` にある。
このsourceは「点が最後のtable familyに入る」という別の係数入力を勝手に代替しない。

## 5. 全被約交叉の二成分

rank4部分はNに入り、rank≤2部分はproper模型と有限表によってN∪Rに入る。
逆包含も明らかだから全閉点で一致し、被約閉部分として等しい。
Rは既約19次元、Nは既約27次元で、上の点によりR⊄Nである。
次元差からN⊄Rでもある。
したがってこの二つが正確な全既約成分となる。

これは既知の二つのPoisson成分の **交叉の成分数**であり、
全Poisson空間の既約成分数や下界57を変更する結果ではない。
全交叉のgeneric-rank≤2点でのscheme structureは別問題であり、
指定rank4点Qでの平方零nilradicalの定理とも対象を区別する。
