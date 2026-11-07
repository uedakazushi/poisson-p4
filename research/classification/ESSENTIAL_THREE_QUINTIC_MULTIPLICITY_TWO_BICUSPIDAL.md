# 最大重複度二・二尖点free quinticの全一候補への還元

2026-11-07、PRIMARY E。対象は有理cuspidal quinticで、最大尖点重複度二、尖点が二つのもの。
free (2,2)から有理・全unibranch・全μ=τ=12が従う既済定理を入力とする。
本稿は有限degree5 linear seriesの全係数を扱い、Nambaのリストに依存せず、
A8+A4の一つのPGL3型だけが残ることを示す。
Saito核・Poisson全周囲・成分性・既知帰属の新計算は行わない。

## 準備

重複度二のunibranch cuspはA_(2k)、δ=kである。
全δ=6なので、二尖点の候補は

    A10+A2,   A8+A4,   A6+A6.

高いδの尖点をp、もう一方をqとする。
k≥3の尖点の接線との交点重複度は四である。
実際、三ならPuiseux pair(2,3)でA2、五なら(2,5)でA4、六以上はBézoutに反する。
A4の接線重複度は四又は五、A2では三である。

pとqを結ぶ直線はpの接線にはなれない。
A8+A4又はA6+A6ならqの接線にもなれず、残りの交点一つはp,qと異なる。
A10+A2に限り、qの接線である場合を別に残す。

以下、正規化のpをs=0、qをt=0、targetでp=[0:0:1]、q=[1:0:0]に取る。
同時零を持たない三degree5 sectionは、imageが直線でなければdegree五のbirational mapとなる。
これはmap degreeとimage degreeの積が素数五であることによる。

## 1. 結線の残交点が相異なる全chart

残交点をs=tに取ると、結線のsectionはY=s²t²(s−t)である。
Xがpでorder二以上、Zがqでorder二以上であることと、Yの倍を引くtarget shearから、
全sectionは

    X=s⁵+A s⁴t+B s²t³,
    Y=s²t²(s−t),
    Z=t⁵+C st⁴+D s³t²                            (1)

に正規化できる。
pでX+B Yのorder三係数はBなので、pがA6以上ならB=0である。
さらにpの接線重複度四からA≠0となる。
同様にqでZ−D Yのorder三係数はDであるから、qがA4以上ならD=0である。
qがA2の場合にはD≠0も残す。

## 2. A8+A4は唯一のnormal form

B=D=0でs/t=rとし、pのaffine座標をx=X/Z、y=Y/Zとする。

    x=r⁴(A+r)/(1+C r),   y=r²(r−1)/(1+C r).

直接

    x−A y² = r⁵[T+(C−A)r]/(1+C r)²,
    T=AC+2A+1.

A8に必要な最初の奇数次五を消す条件はT=0である。
K=C−Aとすると、さらに三次のyの項で偶数次六を消して

    x−A y²+K y³
      =K r⁷[(C+3)−3r+r²]/(1+C r)³.

従って次の奇数次七が消える必要条件はK(C+3)=0である。
C=AならT=(A+1)²=0なのでA=C=−1となるが、
このときX,Y,Zはs−tを共通因子に持ち、degree五のcurveではない。
残る枝はC=−3、T=0からA=1である。

従って唯一の候補は

    (s⁵+s⁴t : s²t²(s−t) : t⁵−3st⁴).               (2)

三sectionはbasepointを持たず、imageはdegree五である。
pではδ少なくとも四、qでは下の奇数次五の係数が
AC+1+2C=−8≠0なのでA4、δ=2となる。
genus budget六からpは正確にA8で、他の特異点は無い。

## 3. A6+A6はbasepoint枝だけとなり空

再びB=D=0とする。
pの最初の奇数次五が消える条件はT=AC+2A+1=0であり、
qの同じ条件はT'=AC+2C+1=0である。
差を取るとA=C、従って(A+1)²=0となる。
上と同じA=C=−1の共通因子枝しかなく、degree五curveは存在しない。

## 4. A10+A2：結線の残交点が相異なる場合

ここではB=0、Dは自由に残す。
den=1+C r+D r³ とし

    x=r⁴(A+r)/den,   y=r²(r−1)/den.

A≠0である。
まず奇数次五を消す条件T=AC+2A+1=0から

    C=−2−1/A.

K=C−Aと置くと

    x−A y²+K y³
      =r⁶[(K+ADr+Dr²)den+K(r−1)³]/den³.

括弧内のr係数はE=K(C+3)+ADである。
A10に必要な奇数次七も消すと

    D=−K(C+3)/A=(A+1)²(A−1)/A³.

括弧内のr²,r³係数をL2,L3と書けば

    L2=D(AC+1)−3K,
    L3=D(2C−A)+K.

偶数次八をL2*y⁴で消した後の奇数次九の係数は

    L3+(C+4)L2 = 2(A+1)²/A⁴.

A10ならこの奇数次九も消えなければならないのでA=−1となる。
しかしこの場合C=−1,D=0で、三sectionはs−tを共通因子に持つ。
従ってこのchartではA10+A2は不可能である。

この比較は有限degree5 sourceのPuiseux型を決定するための係数比較であり、
Poisson一法線の高次jetを反復した計算ではない。実計算は全て上記有限有理恒等式である。

## 5. A10+A2：結線がqの接線である残chart

この補集合ではYのdivisorは2p+3qなのでY=s²t³である。
全sectionを同じtarget shearで正規化し、pの型からorder三係数を零にすると

    X=s⁵+A s⁴t,
    Y=s²t³,
    Z=t⁵+C st⁴+D s³t²,
    A≠0, D≠0

となる。D≠0はqの重複度二を保証する。
pのchartでは

    x=r⁴(A+r)/den,   y=r²/den,
    den=1+C r+D r³.

最初の奇数次五を消す条件はAC+1=0であり、C=−1/Aとなる。
続いて

    x−A y²−C y³
      =r⁶[(C+ADr+Dr²)den−C]/den³.

奇数次七を消す条件C²+AD=0からD=−1/A³となる。
その後の括弧のr²係数はD(AC+1)=0だが、r³係数は2CD≠0である。
従って奇数次九が非零であり、A10ではない。

結線がpの接線という補集合はBézoutで不可能であり、これで全chartを扱った。

## 結論と次のE義務

最大重複度二・二尖点の全free quinticは、(2)のA8+A4型一候補へ入る。
MoeのC5やNambaの引用に依存せず、全parameterと結線接触の補集合を処理した。
この候補の全18Saito核、全cotangent源、ambient supportと成分帰属はまだ実行・採用していない。
全磁気除去定理により、固定Fの完全なE義務はこの一つの曲線・一軌道だけへ還元される。

最大重複度二の一尖点A12、および三尖点・四尖点の候補は別に残す。
本稿で新成分を追加しない。
