# R19のgeneric B16点：quadric normalizerを持つ全完成像

2026年11月3日。被約Rとその被約境界の局所環だけを扱う。
全Poisson空間または実CA/CB scheme交叉の境界厚みは扱わない。

## 定理

B16の非空稠密PGL5軌道上で、被約Rの完成局所環は正確に

    A = C[[u1,...,u16,t,r,c,p²,pr,pt]]
      ⊂ Bnorm = C[[u1,...,u16,p,t,r,c]]/(r²−25t²−2cp)

である。右辺はnormalizationで、三次元rank4 quadric coneと十六正則因子の積である。
normalizerを滑らかとは仮定しない。

Aの次元19、embedding dimension22、depth18、multiplicity3であり、analytic domainである。
normalization quotientは

    Bnorm/A = C[[u1,...,u16,a=p²]]*p.

conductorはnormalizer内では(t,r,c)Bnorm、target内では(t,r,c,pr,pt)Aである。
normalization fiberは一点で長さ2である。
B17の完成像はC[[u16,a]]、B16の完成像はC[[u16]]である。
被約境界B18のA内idealは(t,pt)である。

## 1. generic点と全proper fiber

指定点を

    P0=∂y∧D0,
    D0(h)=D0(x)=0, D0(z)=−x², D0(w)=h²+z²

とする。全physical GL微分のaffine rankは17、projective orbit次元は16である。
既存のB16次元は16なので、この軌道はB16の稠密開である。

L=C∂yはintrinsic factorとして一意である。
x=0に制限するとD(w)=h²+z²が非零であり、定factor候補の商vはw方向に限られる。
次にD(z)=−x²でv=0を得る。
定covector Casimirsはspan(h,x)である。
全spectral持ち上げはt=0上にあり、XがD(z)=−x²を割る必要からX=Cxに限られる。
D(w) mod xが非零なのでK=span(h,x,z)に一意に定まる。
よってproper sourceの全geometric fiberは一点であり、近傍で有限化できる。
既存のnormal spectral sourceはその近傍の有限normalizationとなる。

指定source点はt=r=b=c=d=0,e=−1で、相対quadricのvertexに位置する。
このことがsmooth normalizer型のpinchを仮定できない理由である。

## 2. 全有限source sliceとimage代数

四変数(p,t,r,c)について、spectral族を

    b=0,d=−2p,e=−1,G=h²+z²,
    r²−25t²−2cp=0

と置く。さらにg=I+pE_xhを作用させる。
E_xh²=0なのでinverseはI−pE_xhであり、全150係数は有限多項式となる。

v_iをtarget係数から中心値を引いたものとし、a=p²、b1=pr、d1=ptと置く。
以下は厳密な恒等式である。

    t=v128/2,
    r=−v112−13*v128/2,
    c=−v22,
    a=v105,
    b1=(v15+v108)/2,
    d1=(v108−v15)/26.

逆に全150出力はt,r,c,a,b1,d1の多項式になる。
normalizerの関係式2cp=r²−25t²を使った全出力のexact reductionを保存した。
従ってsliceの全完成像は正確に

    A0=C[[t,r,c,p²,pr,pt]] ⊂ C[[p,t,r,c]]/(r²−25t²−2cp)

である。p²を既に含むためnormalizerはA0の有限拡大である。

source variable pを除去するideal

    (p²−a, pr−b1, pt−d1, 2cp−r²+25t²)

のlex Groebner basisから得た十一式をreportに保存した。
elimination theoremにより、これはt,r,c,a,b1,d1による**全target kernel**である。
単なるpullback-zero式の列挙ではない。

## 3. full R完成像であること

projective chartをu110=−1とする。
指定点で独立な十六projective GL方向を選び、その有限積で局所GL chart H16を作る。
これらとsliceのt,r,c微分はtargetでrank19を持つ。
p微分のtarget像は0だが、source旗X=x−p*hは動くため、sourceでは非零である。
全tensor係数と旗Pluecker係数を合わせた二十列のrankが厳密に20であることを確認した。

quadricを課す前のsource ambientは滑らか20次元である。
従ってH16×四源sliceからこのambientへの射はetaleとなる。
相対quadricのpullbackはunit倍のr²−25t²−2cpである。
これを課したnormal sourceの完成は、十六正則因子と上のquadric coneとの積になる。

全H作用後のtarget出力はA0[[u16]]に入る。
一方、十六orbit cotangent列と六つのminimal algebra生成元

    t,r,c,p²,pr,pt

の係数ベクトルを合わせたrankは22で、minorは6144である。
従ってambient完成環からA0[[u16]]への写像はcotangent上全射で、complete Nakayamaにより全射となる。
有限normalizationへの完成注入と合わせ、実際の全R完成像は記載したAに等しい。

## 4. quotientとconductor、depth、multiplicity

normalizerはA+A*pである。
tp,rp,cp=(r²−25t²)/2は全てAに入るため、quotientはpにより生成され、t,r,cに殺される。
t=r=c=0へ制限したnormalizerはC[[u16,p]]で、target像はC[[u16,p²]]である。
よってこれら以外の関係はなく、正確に

    Bnorm/A=C[[u16,a]]*p

を得る。従って記載した二つのambientでのconductorが定まる。
normalizerはquadric hypersurfaceなのでCM19次元、quotientはCM17次元である。
正確な列からdepth A=18となる。
六つの横断minimal generatorと十六正則変数からembedding dimension22を得る。

multiplicityは有限birational拡大のnormalizerで、延長ideal

    I=(u16,t,r,c,p²)

について計算してよい。a=p²、v=c+a、q=r²−25t²と置くと

    a³−2v*a²+v²*a−q²/4=0.

したがってaはQ=(u16,t,r,c+p²)のideal integral closureに入り、QはIのreductionとなる。
normalizerのQによる商はC[[p]]/(p³)で長さ3である。
normalizerはCMなのでmultiplicityは3、targetのmultiplicityも3となる。
normalizationのtarget maximal idealによるfiberはC[p]/(p²)で長さ2である。

## 5. 既存generic B17像との整合と範囲

a≠0でnormalizerはp=±sqrt(a)の二つの局所branchへ分かれ、
各branchではc=(r²−25t²)/(2p)を消せる。
それらはt=r=c=pr=pt=0の滑らか17次元部分に沿って交わる。
これは既存のgeneric B17 codimension2二枝モデルと一致する。
a=0で両normalization点が合流し、source自身もquadric vertexとなる。

boundary Cartier divisorのsource上の式はt=0である。
Bnorm=A+A*pだから、そのtarget contractionはtBnorm∩A=(t,pt)Aであり、
この商はnormalizerのt=0 quotientへの像なのでreduced domainである。
これが被約B18の全完成像である。

* `code/ambient_R19_B16_completed_image_20261103_v63.py`
* `reports/AMBIENT_R19_B16_COMPLETED_IMAGE_20261103_v63/report.json`

全QQ計算0.78秒。有限150出力、十一式の全kernel、22-rank algebra cotangent、
20-rank ambient source chart、multiplicity reduction恒等式を保存した。
本定理はB16の指定した稠密開軌道についてであり、B16の全退化点の分類ではない。
