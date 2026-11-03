# generic B18でのactual N完成環：odd idealの等号と全81関係

2026-11-03。対象は**被約target N自身**のgeneric B18点での完成環である。
source ramification schemeや、共通29座標へのpartial projectionとは異なる。
actual CA/CB交叉Dの全完成環、CA/CB個々の完成型、B18の全特殊点は本稿の範囲外である。
下記の無条件化にはactual odd2jetの独立PASSを必須とし、未受理の標本からは結論しない。

## 1. actual I=I+(nu)をNakayamaで証明する

既証のactual finite imageをA=R⊕I*t、B=R[t]/(t²−nu)とする。
Rはregular28次元の実target底である。
J=I+(nu)については、source critical quotientとregular ambientのminimal presentationから

    (R,J)≅(C[[g19,nu,s1,...,s6,w,v]], (nu,s1,...,s6,wv,v²))

という**全formal ideal pair**が独立受理済みである。
従ってJ/m_R Jは9次元であり、七つの一次classと、二つのcore二次classwv,v²に分かれる。

新しいactual odd2jetの計算は、I内の七つの一次classがJの全七一次classを張ることを示す。
さらにactual odd係数からその七一次基底をsubtractした二元は、Jの全一次部分のkernel plane上で

    −4wv, −6v²

という独立な二次制限をもつ。
m_R Jの二次項はこのplaneに制限すると零になるので、これら二classもJ/m_R Jで独立である。
ゆえにI→J/m_R Jは全射であり、J=I+m_R Jとなる。
有限R-module J/IへNakayamaを適用するとJ/I=0である。

結論は**I=J**、特に**nu∈I**である。
これはmu(I)≥9だけから等号を推測したものではなく、より大きい既知ideal Jの全minimal classを
actual Iの元が張ることを使った有限moduleの証明である。

nuを固定したままR/(nu)の座標をRへ持ち上げることで、

    R=C[[g19,nu,s1,...,s6,w,v]],
    I=(nu,s1,...,s6,wv,v²)

と取れる。nu方向の高次展開を打ち切ってideal型と呼んだものではない。

## 2. actual全完成環とfinite81関係

nu=t²をBで消すと、normalizerはB=C[[g19,w,v,s6,t]]であり、actual imageは

    A=C[[g19,w,v,s6,t²,s1t,...,s6t,t³,wvt,v²t]] ⊂ B.

これがgeneric B18での全cone N完成環である。
非零scalar座標を最初に固定して同じ議論を行うと、projective版のg変数は18個になる。
実際、選択Rには値1のscalar座標row106が含まれ、sourceのflag parameter tはtensor scalarで不変、
odd coefficient idealは単元scalar倍だけ変化する。従ってそのsmooth scalar因子を独立に除ける。

有限presentationはR0=C[[nu,w,v,s1,...,s6]]、九つのmonomial

    b=(nu,s1,...,s6,wv,v²)

と九変数Y_i=b_i*tで与える。
kernelは次の**全81関係**で生成される。

1. i<jの36本のTaylor pair syzygies:
   lcm(b_i,b_j)/b_i*Y_i−lcm(b_i,b_j)/b_j*Y_j.
2. i≤jの45本のproducts:
   Y_iY_j−nu*b_i*b_j.

productsで任意の多項式はR0+ΣR0*Y_iへ戻る。
その残るkernelはrank2基底1,tの独立性により、R0^9→Iのfirst syzygiesだけである。
monomial idealのfirst syzygiesは上の全Taylor pair syzygiesで生成される。
従って単に81関係がnormalizerで零になるだけでなく、これらが**全kernel**を与える。
complete power series版でもfinite R0-moduleへの完成により同じ結論が成り立つ。
19または18のsmooth変数を加えてcone/projective全環を得る。

## 3. dimension、depth、重複度とdomain

R/I=C[[g19,w,v]]/(wv,v²)はdim20・depth19である。
したがって0→I→R→R/I→0からdepth_R I=20、A=R⊕I*tからdepth A=20となる。
embedding dimensionはdim R+mu(I)=28+9=37であり、dim A=28である。
projective完成環はdim27・edim36・depth19となる。

odd minimal generatorの二乗はnu*b_i²∈m_R³に入る。
従ってm_R*Aはm_Aのreductionで、rank_R A=2からHilbert--Samuel重複度は2である。
CMではないので、CM minimal-multiplicity boundを満たす必要はない。

Aはregular domain Bのsubringで、t=(t³)/(t²)がfraction fieldに入る。
BはA上有限でnormalなので全normalizerである。
特に**actual N完成環はdomain**であり、embedded associated primeを持たない。
後述するconductor quotientのembedded primeをA自身へ移してはならない。

## 4. conductor、subintegrality、local bad locus

normalizer内のconductorは

    C_B=(t²,s1,...,s6,wv,v²)*B.

target A内のcontractionは、九つのeven I生成元と九つのodd I*t生成元を合わせたidealである。
R-moduleとしてC_A=I⊕I*tであり、B/A≅(R/I)*tとなる。
conductor quotient A/C_A≅R/Iにはembedded nil class vがあるが、A自身はdomainである。

t²とt³はAに入り、tはAに入らない。
B=A[t]はsquare/cube criterionによるelementary subintegral extensionである。
従ってnormalizationは有限universal homeomorphismで、seminormalizationはB、Aはseminormalでない。
特にこのformal germでは、いかなる次数でも二つの相異なるgeometric normalizer点や
reduced formal curve arcが同じactual N像を持つことはない。
非被約fiberのscheme長2やinfinitesimal ramificationが消えるという意味ではない。

conductorの被約台はsmooth C[[g19,w]]、cone次元20・projective次元19である。
外側ではA=Bがregular、conductor上ではIのheightが2以上でprincipalになれないので、
有限moduleのCM→freeの議論によりCMでなく、またnormalでもregularでもない。
従って**このformal germ内**のSing、Nonnormal、NonCMはいずれもconductor supportに等しい。
global Nのbad locusを分類した結果ではない。

* code: `code/ambient_B18_N_full_completion_20261103_v63.py`
* evidence: `reports/AMBIENT_B18_N_FULL_COMPLETION_20261103_v63/`

独立actual2jet PASSをgateとし、全81関係とnormalizer map、別Macaulay2 elimination入力を固定する。
