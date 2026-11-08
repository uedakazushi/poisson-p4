# moving-simple残二源：HBは既知99、HAは一つの実軌道

2026-11-08、`/root/quartic_review`。純証明、新しいoperator又は点計算無し。
旧whole-source算術 `PLANE_NONREDUCED_M3_MOVING_SIMPLE_SEPARATE_REVIEW.md` の
全normalizer逆及び全W十八affine解を入力とする。
HB_0,t=0の**全七源は既知Ga9-distinct99成分に入る**。
HA_1,t=2/15の**全六源は一つのactual PGL軌道に入る**。
後者の軌道次元の下界、成分性又は高次元親帰属は本純証明には含めない。

全座標は(a,b,c,z,w)、Π=Q+∂c∧W、y=(a,b,z,w)。
旧full-source座標同型は13P1＋FULL9P2からQ四＋W十八へのdet−1の同型である。
従って本稿のWはhomogeneous P2 kernelだけでなく、元の全affine particularも含む。
全normalizerのunit、零補集合、元α壁は旧独立証明のまま保持する。

## 1. HBの全W七と実D

t=0では

    Qab=−aw/2+bz, Qaz=0, Qaw=zw,
    Qbz=zw/2, Qbw=w²/4, Qzw=0.

旧reportの完全kernel順でβ0,...,β6を取り、m=β5とする。
全particularはW0=(a²,−2ab,az,2bz)で、全解は

    Wa=a²−4β4 bz+β4 aw−12m bw+4β3 zw+β0 w²,
    Wb=−2ab+(4β1+3β6)bz+β1 aw−β4 bw+β2 zw+β3 w²,
    Wz=az+β4 zw+m w²,
    Ww=2bz+β6 zw.

全七係数は無条件にPoissonである。これは旧350affine方程式の完全解であり、
新しい必要条件を十分条件と読み替えたものではない。
実DはD0に七つの次式をβ_i倍で加えたものになる。

    D0=−2ab wz²−2b²z³+c w²z²/2,
    Dβ0=w⁴z/4,
    Dβ1=a w²z²+4b wz³,
    Dβ2=w²z³,
    Dβ3=2w³z²,
    Dβ4=−2b w²z²,
    Dm=−a w⁴/4−3b w³z,
    Dβ6=2b wz³.

この七方向は、aw⁴,w⁴z,aw²z²,bwz³,w²z³,w³z²,bw²z²の係数から独立。
したがって固定Qと既Euler gauge上の全W七→Dはinjectiveである。

## 2. HBの全六parameter除去と零補集合

a'=a+αz,b'=b+δw,z'=z,w'=wはliteralにQを保つ。
c'=c+ell(a,b,z,w)もQを保ち、Wだけを変える。
必要なEuler W→W+linear(z,w)·(a,b,z,w)で旧gaugeを一意に戻す。
全mapは三角、det1、逆は明示的substitutionであり、actual projective actionである。

canonical W(m)はβ_i=0 (i≠5)とする。
canonical Dをinverse shearで置換したとき、c-translationで変えられない係数は

    bwz³: 2(α+2δ),
    w⁴z: m(α+12δ)/4,
    aw⁴: −m/4.

目標全Dにおいてs=2β1+β6とすると、m≠0上で

    δ=(β0/m−s)/10, α=(6s−β0/m)/5

により前二係数を合わせられる。
残りのaw²z²,bw²z²,z³w²,z²w³の四係数は、
既知α,δの項を引いた後、ellの四係数が各−1/2倍で独立に合わせる。
全W→Dがinjectiveであり、群作用が全affine Poisson sourceを保つので、
これで元tensorそのものがW(m)へ戻る。単なるdivisor一致からの推論ではない。

diagonal(a,b,c,z,w)→(s a,s b,s²c,z,w)、tensor scalar1/sは全Q/W0を保ちm→sm。
s=−1/(20m)でm=−1/20へ戻る。unitはmだけで、全inverseは有理で明示的。
m=0は全Dがzで割れる。さらに全七Poisson sourceはaffine空間であり、
m≠0は稠密なので、既知閉成分への包含が分かれば**全m=0も閉性で含まれる**。
逆向きの境界推論は用いない。

## 3. HB canonical tensorは全十括弧で既知99と一致する

m=−1/20のWは(a²+3bw/5,−2ab,az−w²/20,2bz)。
実可逆座標変更を

    (A,B,C,D,E)=(z,w,−a/2,b,−c/12)

とする。known99 reportのcanonical_eigen_coordinate_bivectorをKとすると、

    K = −240 g_*Π + Eul∧L,
    L=(0,0,48A,−72B,−8C).                            (1)

Kの全非零括弧は次で、省略ABは零。

    KAC=48A², KAD=48AB, KAE=32AC+B²,
    KBC=−72AB, KBD=−12B², KBE=−40AD−8BC,
    KCD=72AD+48BC,
    KCE=−48AE+6BD+32C²,
    KDE=72BE−88CD.

(1)はこの十括弧への直接代入で確認出来るfull150同値である。
したがって既受理のGa9-distinct99成分にcanonical Πが入り、§2よりHB全七源が入る。
元Q-source及びfirst-sourceのactual normalizersを戻しても同じ閉PGL安定成分への帰属が成り立つ。

## 4. HAの全W六と全六係数除去

t=2/15、k=−9/14で

    Qab=w², Qaz=0, Qaw=kz², Qbz=z², Qbw=2zw, Qzw=0.

全particularは

    W*=(a²−27bw/7, −2ab, −2az−3w²/5, 2aw/5+27bz/5).

実D*は

    −(3/5)w(w²+3az)²−(27/14)b w²z²−(27/5)abz³−(9/14)c z⁴.

旧完全W六のD imageは、D*に以下の独立六方向を加えたもの。

    r0[−9awz³/14−81bz⁴/196−9w³z²/56]
    −9r1z⁵/14−9r2wz⁴/28−3r3w²z³/14
    +r4[−az⁴−w²z³/3]+9r5bz⁴/14.

awz³,z⁵,wz⁴,w²z³,az⁴,bz⁴から独立で、従ってこの全W六→Dもinjective。
literalにQを保つactual shearsは

    a'=a+kδz, b'=b+γz+δw, c'=c+ell(a,b,z,w), z'=z,w'=w.

inverse substitution後のD*のawz³係数は−54δ/35、
w³z²係数との比は1/4で、r0方向と一致する。
従ってδ=5r0/12で全r0を合わせられる。
続くw²z³係数はγにより27/14倍で動くので、既知δの項を引いて全r3/r4残差を合わせる。
残az⁴,bz⁴,z⁵,wz⁴はellの四係数が各9/14倍で独立に合わせる。
符号はc'=c+ellのinverse substitutionにより決まり、いずれも非零定数。
全逆はpolynomial、parameterの零やdiscriminantを除く必要は無い。
injective D mapとactual group actionを合わせて、**全W六はW*の一実軌道**と結論する。

D*はprimitive c-linearで、先頭係数z⁴と残項のgcdはw⁵項のため1。
従ってC上絶対既約。∂cD*がz⁴の非零倍、z=0上の∂wD*がw⁴の非零倍なので、
被約特異支持は唯一の平面z=w=0。
一般normal cubicは−27a z²(aw+bz)/5で、multiplicity21である。
この存在確認と全軌道還元から、成分性又は軌道次元の下界はまだ推論しない。
旧G image≤23は上界としてのみ保持する。

## 5. 次の有限義務とscope

本稿は旧input/kernel/sourceを変えず、追加数値演算を行わない純証明である。
HBについてfull10 equalityと全parameter inverseは別作者が直接検証可能。
HAはcanonical W*一点の元全150variation/350Jacobiと全actual群像を一回だけ検証すれば、
成分性又はより大きな親への必要normal quotientを判定出来る。
新first symbolやnormal profile系列は不要であり、その次の有限jobは別契約にする。

moving-double degree(1,0)、generic-triple、m二、他zero cubicをこの包含へ混ぜない。
