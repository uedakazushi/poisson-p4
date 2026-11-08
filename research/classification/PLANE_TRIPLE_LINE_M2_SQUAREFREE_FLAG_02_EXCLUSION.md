# R=3L・m二squarefree・flag02は全空：二つの固有値比の矛盾

2026-11-07、root純証明、独立読解前。新CAS・matrix・point・jet計算は無い。
actualD被約既約quintic、rank二Poisson plane、零cubic f=a³、normal位数二、
generic squarefree normal quadraticでprimitive root degree flag(0,2)を仮定する。
normaldouble、flags00/01、f=a²b、m三等は別残余である。

## 1. 全connectionとraw fields

`PLANE_NONREDUCED_M2_SPLIT_ROOT_ENTRANCE.md` の全極予算から、

    G2=gamma*a*u*(A*u+a²v), gamma≠0,
    M=[[1,0],[A,a²]],  A homogeneous quadratic, gcd(A,a)=1,
    eta=theta1−theta2=−2 dlog a.

二rootのscaleはactual normal GL2とgammaに吸収している。全primitive Aと全codim二basepointを保持。
beta=diag(theta1,theta2)をMで戻すと、全diagonal entriesは同一thetaで

    beta11=beta22=theta, beta12=0, beta21=−d(A/a²).

thetaはdlog a+d(P2/a²)、P2全quadratic modulo Ca²で、全nonlog exact係数を許す。
tracefree homogeneous base-field gaugeの全P0は

    P_base=P_f,
    P_base,u=uX,  P_base,v=uZ+vX,
    P_uv=k*u²+kappa*u*v+k02*v².

X,Zはtracefree linear fields。全higherは忠実P1三十・P2九を保持する。
connectionとfieldの正しいtwistは、Omega=da db dc、E=a∂a+b∂b+c∂cとして

    iE iX Omega/f=dlog f−3theta,
    iE iZ Omega=−3f beta21=3a dA−6A da.

第一式のclosed性及びdivX=0からX(f)=0、従ってXa=0。
第二式とtracefree条件を全三係数で解くと

    Za=0, Zb=−3 A_c, Zc=3 A_b.

余分なEuler multipleはtracefree条件で零になる。
よってXのmatrixはa-row零で、(b,c)block Hのtraceは零。a-columnの二translationsも除外しない。

## 2. 全actualG2の二係数を比較する

全Euler-signed Pfaffianのv²係数は3f*k02、targetは零なのでk02=0。
uv係数はX=Yからmixed cross項が消え、正確に3f*kappa。
従って **gamma=3kappa、kappa≠0**。このunitはtarget squarefreeから導いた。

u²係数は

    3f*k − det(E,X,Z) − Z(f).

Za=0なのでZ(f)=0、上のZ式からdet(E,X,Z)=3a*X(A)。
これがtarget gamma*a*Aと一致するので

    X(A)=a²*k−kappa*A.                       (1)

これは全P0の必要恒等式であって、flatconnectionだけから存在を仮定していない。

## 3. 全higherに必要なbad Cとcurl

fullP1のu-columnのuで割れない部分はv²*C、Cはconstant base vector。
P2にnormal行は無い。C=0ならfullu-rowがu倍となりactualDにu因子、既約性に反する。
従ってC≠0が必要。

既元projective Jacobiの一normalweight係数から得る普遍curl恒等式

    (2kappa*Id+4DY−DX) C=0

はRの被約性を使わない。P_fの項はgradient(C(f))でcurl零、Zとkはu因子のため寄与しない。
今Y=Xなので

    (2kappa*Id+3DX) C=0.                     (2)

Xa=0とkappa≠0によりCa=0。従ってHはnonzero eigenvalue−2kappa/3を持つ。
trace H=0からeigenvaluesはlambda,−lambda、lambda≠0であり、

    lambda=±2kappa/3.                        (3)

Hはdistinct eigenvaluesでdiagonalizable。
a-columnのtranslationはH可逆なのでbase shear(b,c)→(b,c)+a*constantで除去できる。
このshearはa及びfを保ち、inverseもregularなactual GL3である。
さらにactual GL2でH=diag(lambda,−lambda)へ。全AとCも同時に変わり、primitivityは保たれる。

## 4. primitive Aが要求する異なる比

式(1)をa零で取る。gcd(A,a)=1なのでAbar(b,c)は非零binary quadraticである。
H(Abar)=−kappa*Abar。
Sym²(H)のeigenvaluesは2lambda,0,−2lambda。
kappa≠0だから零eigenvalueは使えず、

    lambda=±kappa/2.                        (4)

(3)と(4)はkappa≠0で両立しない。
従ってC≠0は不可能。C0では既約actualDも不可能なので、指定flag02全higher sourceは空。
nonsemisimple/zero Hは(2)から既に除外し、未知parameterをgenericにしていない。

この固有値比証明は原projective一式とactualG2の二係数だけを使う。
全Poisson非被約scheme、全E又は他flagsの完了を主張しない。
独立レビューでは上のfield/character factor3、actualuv=3f*kappa、curlの係数2対3を特に確認する。
