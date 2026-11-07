# rank二plane・R三重直線・m三squarefree leadingの全first connection

2026-11-07、主義務G。producer手証明。ここでは全first symbolを二つの実正常形へ還元する。
higher30/9の全G上界は後続finite jobで判定する。m二、m三重根、m四及びR=2L+Mは別補集合。

## 1. 非被約Rにも使えるparallel tensor入口

S=P2⊂SingD、π|Sがrank二、実Dが既約被約quintic、S沿いの一般重複度が三とする。
πSの零cubicをR=3Lと仮定する。
generic rank二点のWeinstein分解とD∈I_S²から{I_S,I_S}⊂I_S²。
従ってN*=O(-1)⊗C²にはflat cotangent-algebroid connectionがあり、
U=S\L=A2ではordinary flat connectionとなる。
これは非被約Rに沿い高位極を持ち得る。log又はsimple poleとは仮定しない。

最初のnormal symbolはG3∈H0(O(2))⊗Sym³C²である。
Hamiltonian associated-gradedの証明によりγ3=G3/fはSym³N*⊗detNのparallel tensor。
この証明は`PLANE_REDUCIBLE_CUBIC_CONNECTION_ENTRANCE.md`のparallel部分と同じであり、
f被約性を使うpole-orderの主張はここへ流用しない。
本noteではgeneric binary cubic G3がsquarefreeである全open stratumだけを扱う。
parallel tensorのdiscriminantはUで消えないので、U上のpolynomial Δ(G3)はunit、すなわち非零定数。

## 2. A2上の三rootの全factorization

U上でG3のroot incidence⊂U×P1はfinite étale degree三である。
複素解析空間A2はsimply connectedなので、そのfinite étale coverは三つのcopyへ分かれる。
各copy→Uはdegree一のfinite mapであり、normal U上のisomorphismとなる。
従って三rootはU→P1のglobal algebraic sections。
Pic(A2)=0なので各rootはprimitive polynomial vector (pi,qi)で表され、
G3=h∏(pi z+qi w)と書ける。

rootの相異性から各pair determinant pi qj−pj qiはUで消えないpolynomial、従ってnonzero constant。
overall hもdiscriminantがunitであることからnonzero constantとなる。
各factorのbase polynomial degreeをdiとすると、polynomial ringがdomainなので
leading homogeneous productの消滅は無く、d1+d2+d3=deg_base G3≤2。

nonconstant factorがただ一つなら、他二つはconstant independent vectorsである。
その二つとのdeterminantが両方constantなら残るvectorもconstantとなり矛盾。
よって次の全二場合だけである。

- 三factorが全てconstant。
- degreeは順序を除き(0,1,1)。

後者でconstant factorをzへ直す。
他二factorのw係数はzとのdeterminantなのでnonzero constantsであり、それぞれ1へscale出来る。
二factor間のdeterminantもconstantなので、両者はw+ell zとw+(ell+c)z、c非零constantとなる。
ellはnonconstant affine-linear function。normal diagonal及びconstant unipotent changeでc=1へ直せる。
一つのchosen polynomial matrixから全sourceを推論せず、全root sectionsとdegreeの和から尽くした。

## 3. 二つの全first-connection正常形

base coordinates(a,b,c)でL=(a=0)、affine chart a=1を取る。
homogenizationとactual base/normal linear transformationsにより、全G3は次の二形で覆われる。

1. a² g3(z,w)、g3は全constant squarefree binary cubic。
2. z(aw+bz)(aw+(b+a)z)のnonzero scalar倍と、その同じconnectionのparallel space。

第二型のbはellのhomogenizationをnormal constant shearで調整したもので、b mod a非零。
元の任意ellへはactual base linear inverseとnormal GL2 inverseで戻れる。

三rootが相異なるbinary cubicのtwisted GL2 Lie stabilizerは零である。
projective stabilizerは有限で、scalarはweight m−2=1で非自明に作用するからである。
従ってparallel γ3から**全**normal connectionが一意に復元され、nonlog接続を追加parameterとして捨てない。
第一型のconnectionはcanonical scalar connection βcan=(1/3)dlog f。
第二型はchart a=1でscalar part零、tracefree part−N db、N=z∂w。
全mixed columns及びD2零を戻したquadratic raw liftの代表は

`P0_scalar=3a² ∂b∧∂c`、

`P0_unipotent=3a² ∂b∧∂c+3az ∂c∧∂w`。

chart a=1ではπS=3∂b∧∂c、{c,w}=3zなので∇∂b(w)=−zとなり上のconnectionに一致する。
normal-normal quadraticの追加三係数はD2=0によって一意に零へ戻る。
両P0は全projective Jacobiを満たしrank二であることを後続finite certificateでも確認する。

固定P0上で許される全regular parallel leading spaceは、scalar型でa²Sym³(z,w)の四次元、
unipotent型で

`<a²z³, az²(aw+bz), z(aw+bz)²>`

の三次元である。
後者はa²g3(z,w+(b/a)z)をregular quadratic coefficientへ戻す条件がg3のw³係数零であることから従う。
squarefree条件はこのspaceのopen条件で、finite dimension判定では全spaceをsupersourceとして保持してよい。

## 4. 実source次元と有限target

固定planeでLの選択は二次元、πSのnonzero amplitudeは一なのでscalar P0全sourceは三次元。
unipotent型ではnonzero nilpotent Nのcone二次元とell mod aの非零二次元に、
Nとellの逆scale一冗長性がある。従って全first symbolは六次元以下。
このcountはfull actual normalized sourceの復元であり、固定anchorのdimensionを流用しない。

全higher sourceは既存global normal filtrationのP1三十、P2九で尽くされる。
P1に全350のweight一Jacobi及びD3が上の全parallel spaceへ入るlinear条件を課す。
任意のP1上のP2 fiberは空又は固定weight二operator kernelのtranslateなので、
全higher dimension dはこれら二kernel次元の和以下。
plane六とtensor projectivization一を戻すimage upperはscalarで8+d、unipotentで11+d。
したがってd≤16又は13なら、それぞれ全specified stratumのG≤24が従う。

一回のfinite jobで両全P0、full350×30/9、全24 D3係数、parallel-space membership、
complete kernels/minorsを検証する。足りなければexact whole linear supportとactual rank四条件を保存し、
大きなrank二部分を理由に成功扱いせず、元Jacobi/essential-variable geometryへrouteを変える。
R=2L+M、m二、m三repeated normal root及びm四はこの二モデルの結論へ含めない。
