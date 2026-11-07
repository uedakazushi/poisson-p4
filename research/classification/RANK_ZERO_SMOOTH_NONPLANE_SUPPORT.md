# smooth非plane rank零曲面：clean・横断型を仮定しないG還元

2026-11-07、rootによる幾何証明案。既受理の個別source証明を統合する。
別作者の確認前には新global claimにしない。全E帰属又は特異曲面の分類は主張しない。

## 1. 滑らかなquadric上で零となる既約quintic Poissonは存在しない

S⊂P4を滑らかなquadric surface、πを一般rank四Poisson、π|S=0、
実Dを既約被約五次とする。
Sのspanは一意のhyperplane H=(ell=0)、I_S=(ell,q2)である。
座標をell=x4とし、q2はx0,…,x3だけのnondegenerate quadraticに取る。

E=Λ²T_P4とするとKoszul列は

    0→E(-3)→E(-1)⊕E(-2)→E I_S→0。

E(-3)=Ω²(2)はBott又はEuler計算でacyclicであり、全sectionは

    Π=ell L+q2 A

と表せる。Lは全linear-coefficient bivector、Aはconstant skew matrixである。
Euler gaugeのkernelを捨てた部分源を選んだのではなく、全H0(E I_S)のsurjectionを使う。

Sはnonplaneなのでsecant補題からDの一般S位数は二。
rank零first jetはFrobeniusであり、conormal idealはabelianとなる。
従ってH={ell,q2}_Πはsymbolic square I_S^(2)に入る。
prime complete intersectionなのでI_S^(2)=I_S²。
projective local generatorsell/h,q2/h²への変更による補正もI²に入るため、
このhomogeneous bracketを使ってよい。

Hの次数は三。ell=0へ制限すると、degree四のq2²項は現れ得ず

    H|H=q2 v(q2)=0、v=A^sharp(dell) はconstant base vector。

q2はnonzero、かつ四base変数でnondegenerateなのでv=0となる。
従ってΠ^{ell,j}は全てellで割れ、HはPoisson invariant hyperplaneである。
H上ではπが三次元接空間に入るため、そのtop wedgeは零。
よってellが実五次Dを割り、既約性に矛盾する。

この排除にはordinary double又はblowupのvertical成分不存在を仮定しない。
quadric coneではvがquadricのradical方向に残り得るので、この証明で除外しない。

## 2. 全smooth nonplane曲面の有限な場合分け

π|S=0でSが滑らかな非plane曲面成分の場合、既存のnormalization-length boundからdegree≤6。
Sがhyperplaneに入るなら、degree≥3は既済prime CI square次数論で不可能、
degree二は§1、degree一はplaneなので射程外である。
以下はSがP4で非退化の場合だけを扱う。

q(S)>0なら、小次数smooth irregular曲面補題によりSはdegree五stable楕円scroll。
全bundleのnoquadric・3regular分解・Bott計算からH0(E I_S)はcanonical q52の一次元だけ。
従ってπはT40に属す。この部分はclean residueを仮定しない。

q(S)=0かつdegree≥5は、linear-normalityとEuler完全源によって実Dが零又はq2²linearとなり不可能。

degree三はminimal-degree分類でsmooth cubic scrollだけ。
全symbolic-square線形制限・完全十二源Jacobi・別作者ideal lift/leading cover証明により、
この全rank零Poisson像のPGL射影次元は二十三以下。
rank二以下の余剰源を排除せず、必要な上界だけを使う。

degree四でlinearly normalでなければ、唯一smooth projected Veroneseでありrank零sourceに非零Poissonは無い。
linearly normalならEuler完全源H0(E I_S)=I_S(2)⊗Λ²C5を使う。
r=h0(I_S(2))が零又は一なら、再び実Dは零又はq2²linearで不可能。
r≥2なら、独立な二quadricのCIはdegree四でSを含み、unmixed性と同じdegreeからschemeとしてSに一致する。
従ってSはCI(2,2)である。
既受理のsmooth CI22 parity証明では全sourceのrankは二以下なので、既約rank四Dは無い。
このsmooth統合には、特異CIに対する新しい全pencil還元を先取りする必要すら無い。

非退化surfaceのdegreeは少なくとも三であり、これで全場合を尽くした。

## 3. 結論と正確な未処理集合

smooth非plane曲面S上で零となるrank四πが既約被約五次Dを持つなら、
πはT40又はPGL像次元二十三以下のcubic-scroll源に入る。
したがって、未知の二十五〜二十七次元G優越成分の一般元は、
smooth非plane曲面上で零にはなれない。

Sの特異点、非normal性、quadric cone、rank零plane、rank二restrictionは射程外である。
特にsmooth非planeのvertical・非node枝も処理したことと、
全非正規quinticのG/Eを尽くしたことを混同しない。
