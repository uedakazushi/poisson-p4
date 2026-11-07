# m4 squarefree special j：全 pencil と全 source の独立review

2026-11-07。主義務G。producerと別作者による全幾何・独立Fraction算術。**PASS**。
rank二planeのzero cubic Rがtransverse conic+line、tangent conic+line、nonconcurrent triangleの三型であり、
実Dのplane沿い位数四、leading binaryquarticがgeneric squarefree、I≡0又はJ≡0とする。
**この全枝のprojective PGL像dimensionはそれぞれ≤24,23,23**。
既受理の非特別j枝と合わせ、この三型の全m4 squarefree枝がG非優越となる。
concurrent、generic重根、m2/m3、非被約R、全Eは別残余である。

## 1. I0の全 stable pencil

quarticをg=au⁴+4bu³v+6cu²v²+4duv³+ev⁴とし、
I=ae−4bd+3c²、J=ace+2bcd−ad²−eb²−c³を用いる。
Gの全linear coefficient imageをW⊂Sym⁴C²とする。
Iは五次元上nondegenerate quadratic formなので、I|W=0ならWitt indexからdim W≤2。
rank一は§3のscalar枝へ送る。

rank二ではgeneric squarefree性からJ|Wは非零binarycubicである。
その各rootは非零quartic I=J=0であり、triple+simple又はfourth powerである。
これはdiscriminant恒等式Δ∼I³−27J²と、multiple rootをu=0へ置いた形
u²(Au²+Buv+Cv²)のI=3(C/6)²からも直接分かる。

root quarticがu⁴なら、もう一つの元Bの係数はe=0、−4bd+3c²=0。
d=0はc=0となりwhole pencilがu³を共通因子に持つため、generic squarefreeに反する。
d≠0ではv→v+tuでcを消せる。I=0からbも零になり、u⁴成分を引けば
whole pencilは〈u⁴,uv³〉である。そのJは非零constantを除き−L₁L₂²。

root quarticが4u³vなら、I-polar条件からd=0、Bのbを引けばae+3c²=0である。
J(sA+tB)=−e s²t−4c³t³。
e=0はc=0で全unstable。e≠0,c=0は前のmixed型へ変数交換で戻る。
e,c≠0では三rootが相異なり、diagonal GL₂とBのscaleで

`A=4u³v, B=−3u⁴+6u²v²+v⁴`, `J=−s²t−4t³`

へ戻る。必要な平方根choiceはC上有限であり、whole supportを減らさない。
これでstable I0 pencilの全分類を得る。**perfect-cube Jの枝は存在しない**。
三相異型のdiscriminant supportは三つのconcurrent linesなので、nonconcurrent Rには入らない。
mixed型は二つの異なるused linesを持つ。

## 2. J0のcommon-kernelとcompressionの両方

JはHankel symmetric matrix
M=[[a,b,c],[b,c,d],[c,d,e]]のdeterminantである。
全singular symmetric3×3 linear space Zを考える。rank二の元をcongruenceでdiag(1,1,0)へ直す。
B=[[p,q,r],[q,s,t],[r,t,h]]についてdet(A+zB)の一次・二次係数から
h=0、r²+t²=0。したがってlinear map B→(r,t)のimageは等方line以下である。

image零なら共通kernel e₃を持つ。
image非零ならt=irと選べ、det B=r²(p+2iq−s)=0がZ上恒等零なので、p+2iq−s=0。
よって全formは共通二plane K=〈(1,i,0),(0,0,1)〉上で零、すなわちcompression型である。
rank一以下しかないlinear spaceは安定quarticを含まない。
このargumentはsample pencilではなく、Z全体へのpolynomial identityを使っている。

Hankel五空間内で、common-kernelの三条件はquadratic qによるSym²→Sym⁴の乗法のtransposeであり、
q≠0の乗法は単射なので三条件は独立。したがって次元≤2。
compressionでもK=〈q₁,q₂〉に対するq₁²,q₁q₂,q₂²は独立なので三条件は独立、次元≤2。
後者の独立性は、依存があれば非定数有理式q₁/q₂がC上の二次方程式を満たすという矛盾で分かる。

dim W=2のstable型は厳密に次の二つ。

- common-kernel quadraticがrank二：q=uv、b=c=d=0、W=〈u⁴,v⁴〉（secant）。
- compression二planeのannihilator quadraticがrank二：K=〈u²,v²〉、a=c=e=0、W=〈u³v,uv³〉。

各quadraticがrank一の場合は、全pencilがtriple-root因子を持ちunstableである。
上の二stable型はendpointのroot multiplicityが四重と三重+単純なのでGL₂で同一にはならない。
従ってsecant-onlyという元の推測は採用せず、**compressionもwhole sourceとして保つ**。
両型のdiscriminantはL₁³L₂³に比例し、二つのused linesは相異なる。

## 3. 全nonconcurrent型から一意な接続へ

既受理のHamiltonian associated-graded補題により、γ=G/fは
Sym⁴N*⊗det Nのparallel sectionであり、Δ(G)の零supportはRの成分に含まれる。
log性は仮定しておらず、supportがR全体に等しいとも仮定しない。

rank一ならG=L g₀、g₀はconstant squarefree quartic。
定数四rootを保つprojective Lie stabilizerは零なのでtracefree connectionは直ちに零となる。
scalar部分はdlog(L/f)+2β=0、すなわちβ−βcan=−dlog(L³/f)/6へ一意に固定される。
I,Jの非零性を要しないため、旧独立scalar reviewの**同じ全source**の上界24,23,23を適用できる。

rank二は少なくとも二つの異なるused linesを要する。
conic+lineはline componentが一つしかないので全てrank一へ戻る。
I0三相異型はconcurrent supportのためtriangleにも入らない。
triangleではfinite permutationで二used componentsをx,y、unusedをz、f=xyzと出来、全残余は

`Axu⁴+Byv⁴`, `4Axu³v+4Byuv³`, `Axu⁴+4Byuv³`

の三つだけである。A,Bの全非零振幅を保つ。
以下の有限supersourceはその零境界も許すので、dimension上界を弱めるだけである。

normal covectorのdiagonal connection βu,βvを(dx/x,dy/y,dz/z)の係数で書くと

|型|βu|βv|
|---|---|---|
|J0 secant|(1/8,3/8,1/2)|(3/8,1/8,1/2)|
|J0 compression|(0,1/2,1/2)|(1/2,0,1/2)|
|I0 mixed|(1/6,1/3,1/2)|(1/2,0,1/2)|

となる。各Euler sumは一。
uᵏv⁴⁻ᵏのtwisted weightは(k−1)βu+(3−k)βvなので、parallel条件は
base coefficientのlog derivativeがdlog f−(k−1)βu−(3−k)βvであること。
全五kをexact代入すると、上の二monomialだけがbase一次coefficientを持てる。
これが全leading spaceであり、chosen amplitudeを固定したsliceではない。

四相異rootを持つquarticのtwisted GL₂ Lie stabilizerは零。
したがってdisplayed flat connectionがparallel条件を満たす以上、実接続もそれに等しい。
off-diagonal/nonlog parameterを仮定で捨てたものではない。

## 4. 正しいraw P0と別chartによる独立確認

τ=β−dlog f/3に対しVτ=diag(τz−τy,τx−τz,τy−τx)(x,y,z)を使い、
Pfのmixed columnsを−uVτu、−vVτvへ補正する。
raw {u,v}を零と置かず、actual D₂=0を課す。
k uvのactual divisorへの係数は3fuvなので一意に決まり、三型のkは1/8,1/4,1/6となる。

別の独立チェックとしてEuler gaugeでPxu=Pxv=0へ戻し、x=1 chartを取った。
{y,z}=3yz、{y,u}=3βu,z yu、{z,u}=−3βu,y zu（vも同様）となる。
D₂=0の直接Pfaffian計算は

`{u,v}=3(βu,y βv,z−βu,z βv,y)uv`

を強制し、chartでの係数は3/8,3/4,1/2。
raw kとの差はEuler gaugeそのものであり、接続のfactor三と符号を独立に確認している。
三P0の全150係数、全350 Jacobi及びactual divisor零を新算術で再展開した。

## 5. 全30/9 sourceと弱い三角上界

新checkerはproducerもarithmetic helperもCASもimportせず、五変数sparse Fraction計算を別に実装する。
weight一raw三十六方向に対しcanonical三十とEuler六の合併rank三十六を確認し、
weight二は全raw九方向を保持した。
全350 Jacobi weight一とD₃=0の全24行から374×30行列を再構成し、
保存全entry、六kernelの全行代入・独立性、nonzero24minorを確認した。

weight二は全350 Jacobiと、許された二monomial以外の全13 leading-shape行から363×9行列Lを作る。
その全entry及びminorを再計算し、独立eliminationで完全kernelも新reportへ保存した。
P1を自由六parameterのまま保ち、Q(P1)の全363行・全二十一quadratic coefficientも照合した。
必要なboundにはquadratic obstruction、initial ideal、GB又は追加jetを使わない。

|whole source|dim P1必要核|rank L|各P1上のP2 fiber上界|fixed P0全source上界|
|---|---:|---:|---:|---:|
|J0 secant|6|9|0|6|
|J0 compression|6|7|2|8|
|I0 mixed|6|8|1|7|

各任意P1でL(P2)+Q(P1)=0のfiberは空又はker Lのtranslate。
rankはP1に依存しない同一operatorのrankであり、generic rankを境界へ流用していない。
更なるhigher Jacobi、D既約性等はこのsupersourceを狭めるだけである。

固定P0の全sourceを保つ実connected subgroupは、base triangle torus二、normal diagonal二、
normal-to-base shear六からなる十次元群である。
normal GL₂全四を引くと非scalar connectionを動かすため誤りであり、ここでは使わない。
後者shearはnormal associated gradedでidentityなのでP0及びleading spaceを保つ。
PGL₅×X→projective tensor spaceのfiberへのこの群の作用は、PGL₅因子への右作用が自由なので自由である。
個別tensorのstabilizerに関係なくfiber dimension≥10。

したがってtriangle三全source像は24+6−10=20、24+8−10=22、24+7−10=21。
rank一scalar triangleの既受理上界23と合わせ≤23、conic+lineはscalar24,23である。

## 6. 固定された証明書

finite job gate PASS後、0.893秒、最大RSS16,028KiBで全算術が終了した。
新たな成分・全点帰属又はfull nonlinear source equalityを主張しない。

|artifact|SHA256|
|---|---|
|code/classification_plane_m4_special_review.py|3c43fe85efcae39dbd0d2f8228271ce5fa291830315e62c9169ac5924b7408f8|
|reports/CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_REVIEW_20261107.json|f2b6ae48d87bb7f88d18ecca5883238e955b275fe0ee2ab91727ceba0d02c82d|
|reports/CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_REVIEW_JOB_20261107.json|e17cf9552094480c1798889bdb5aeaed7fe2af484f0c247f768cd0111d1d1a3b|
|producer CLASSIFICATION_PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_20261107.json|191019ce0ab5f6e734e888300132d55d8d3c610e62a56930464d25a0e947d742|
|producer PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_PENCILS.md|829ad9c841494924eb009bba818bba9718a8a8b090ac735bb7948e26d7b86ba7|
|producer PLANE_REDUCIBLE_CUBIC_M4_SPECIAL_BOUND.md|ed84ad40ce1f42df47e90d7c7216cac31ea366b2422f37f04495ec4d35560c51|

全依存path/hashはjob及びreview reportへ保存した。producer report内のimmutable source snapshotを参照し、旧source/reportを変更しない。
