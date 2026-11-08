# R=2L+M・m二constant double：全irregular源と有限ratio入口

2026-11-08、PRIMARY G/E、`/root/quartic_review` の新純手証明。
**独立review前。** 新CAS、source/operator/minor、jet又はGBは実行していない。
全G/E閉鎖を主張せず、全源を包含する正常形と、その後の最弱有限判定の範囲を固定する。

## 1. 全diagonal/off-diagonal pool

座標(a,b,c,u,v)、restriction cubic f=a²b、normal 位数二、generic double quadratic、primitive root degree零。
actualF は非零被約既約 quintic、Poisson rank四、plane restriction rank二を仮定する。
全moving degree一 root の既空性は別証明を用い、本稿では再開しない。
root u を選び、G2=γh u²、γ!=0、h=a^m b^(3−m)、m=3,2,1,0。
全raw36/Euler6 P1と全9 P2を保つ。

H=diag(1,−2,1)、A=a∂c、B=b∂c と置く。volume convention は既 constant-root note と同じ。

    αH=3d(c/a), αA=dlog(b/a), αB=d(b/a).

全closed ordinary diagonal formsは period と二 exact irregular potentials を保持し、対応する全tracefree fieldsは
λH+pA+qBで尽くされる。
double leading relation θ−σ=dlog(f/h) から、両fieldsの exact 部分は同じで

    δ=3(2−m),
    Y=λH+pA+qB, X=Y+δA.

最初からnilpotentだけに制限していない。
全raw P0は base P_f、mixed uX と uZ+vY、normal k u²+κuv。
normal uv target零より

    κ=det(E,X,Y)/(3f)=δλ.

off-diagonal flatnessは ω=(f/h)r のclosed性。
r=−iE iZvol/(3f)、trace Z=0なので、これは exactly Z(h)=0。
従って全Zはh=a³/b³で各五dim、h=a²b/ab²で各三dimである。
logarithmic formsだけを残す操作はしていない。

全u² coefficient条件は

    3fk−det(E,X,Z)−Zf=γh.                          (1)

R=p+δ と書けば、左側の引き算部分を独立展開して

    det(E,X,Z)+Zf
      =Za[(R+2)ab+qb²+3λbc]
       +Zb[(1−R)a²−qab]−3λab Zc                 (2)

を得る。下表の全coefficient比較はこの式から行う。

## 2. 元full bad-vectorはsemisimple diagonalを一枝だけに制限する

bad u-row v²C、C=(x,y,z)は非零が必要。
C零ならfull u-rowがu倍で、全P2九はnormal columnを持たないためactualFはu因子を持つ。
元全projective式を直接展開すると、次の全linear conditionsを得る。

    λ(δ−3)y=0, λ(δ+3)x=0,
    (2p−δ−2)x+q y+λδ z=0,
    qx=0, λδx=0,
    (p−δ+2)y=0,
    (p+2)x+2q y+λ(δ+3)z=0.

λ!=0なら全δ∈{−3,0,3,6}でx=0。
δ!=3ではy=0となり、残二式からz=0なので空。
δ=3 (h=ab²) だけ、p=1、y!=0、q y+3λz=0 が可能。
λ=0なら常にx=0。q!=0ではy=0、q=0ではp=δ−2のときだけyも許す。
z方向は全periodで残る。このjumpをgeneric値で捨てない。

h=ab²で全Z=ν diag(−2,1,1)+rA+sBと書くと、(1)は

    λν=0, k+ν(p+4)+λr=0, γ=3(νq+λs)!=0.

λ!=0の残枝ではν=0、Zはradical nilpotent、κ=3λ!=0。
次節のgeneralized ratioでこの全枝も空になる。
従って以後の**surviving necessary sourceではλ=κ=0が結論として従う**。

## 3. generalized ratioの全source補題

ρ=v/u、K(ρ)=k+κρ、D(ρ)=Z−ρ(X−Y)をK0=C(ρ)上で読む。
row divergence はRu=κu、Rv=−2ku−κv。
canonical strict Casimir liftを用いると、M=Π#dρは

    M(u,v)=(3K(ρ)/5)(u,v),
    M(base)=(D(ρ)−2K(ρ)I/5)base+T1(u,v), M(F)=0.

T1は全mixed higherから来るrational homogeneous degree一vectorで、ρはM-invariant。
K(ρ)I−D(ρ)可逆なら、全translationをnormal degree一のaffine K0-mapで消せる。
normal degree n、base degree5−nのoperatorは

    D(ρ)|Sym^(5−n)+K(ρ)(n−2)I.                    (3)

n=3/4/5の全operatorも可逆なら、exactly F=γu²h(translated base)。
hはreducible cubicなのでC(u,v)上のpositive-base-degree factorizationをGauss lemmaで戻し、既約性矛盾。
ρ又はK(ρ)の反転はcoordinate fraction fieldだけであり、source parameterの未記録な除外ではない。
K(ρ)=0 identicallyとなるk=κ=0の枝は別に保持する。

特にD(ρ) nilpotent、K(ρ) nonzeroなら全n!=2 operator可逆で空。
前節のλunit/h=ab²枝はZ−ρδAがradical nilpotent、K(ρ)がκunitによりnonzeroなので全空。
このためsemisimple diagonalの排除に追加source計算は不要。

## 4. λ=0の全leading coefficient table

以下は(1)/(2)の全monomials比較であり、Zの自由成分を先に零にしていない。

### h=a³, δ=−3

    Z=(0,A0*a+B0*b+C0*c,D0*a+E0*b−B0*c).

全必要条件は

    γ=A0(p−4)!=0, C0=0, q*B0=0,
    3k=(4−p)B0−qA0.

A0 unit。relative normal scaleでA0=1。
B0!=0ではq=0、actual base c-shearsでD0,E0を消せる（2B0² minor）。
base diagonal/tensor/relative normal scaleでB0=1も同時に出来、pは変えない。
全normalformは Z=(a+b)∂b−c∂c、k=(4−p)/3、p!=4。
D(ρ)は全ρでeigenvalues0,±1なので、(3)の全spectral/translation failuresから

    p∈{5/2,11/2,1,7,−2,10}

の六fixed P0だけが必要残余。

B0=0ではk=−q/3、D(ρ)は全ρでnilpotent。
q!=0は前節で全空。
q=0ではk=0、D0をnormal shear又はbase c-shearで消す。
E0零/nonzeroをactual diagonal/tensor/normal scalingsで0/1へ分けると、

    X=(p−3)A, Y=pA, Z=a∂b+εb∂c,
    k=κ=0, p!=4, ε=0/1.

ε0ではcanonical two affine fieldsのlinear spanはa∂b,a∂c。
p−4 unitからspanは全K上二dimで、commuting affine-pair lemmaによりF_b=F_c=0。
従ってこの全p源は既essential-three八成分へ全点被覆される。
ε1の全p!=4源は未閉として残す。

### h=a²b=f, δ=0

    Z=νH+rA+sB,
    qν=0, γ=3(k−pν)!=0.

ν=0はkunit/nil ratioにより全空。
ν!=0ではq=0、actual c-shiftでsを消す。
relative normal/base diagonal/tensor scalesでν=1、r=ε=0/1へ戻せ、pは固定される。
全modelは X=Y=pA、Z=H+εA、γ=3(k−p)!=0。
(3)の全eigenvalues（nil central part εAは消さない）から全必要kは

    k∈{0,−1/2,−2,1,4}, p!=k, ε=0/1.

k=0はfield lemmaで反転せず、そのwholep源として残す。
これが十のwhole oneparameter p-linesである。

### h=ab², δ=3

残λ=0枝は

    Z=ν diag(−2,1,1)+rA+sB,
    νq!=0, k=−ν(p+4), γ=3νq.

actual c-shift/normal shearでrを消す。
base a/b比、relative normal scale、base c-scale/tensor scaleによりq=ν=1、s=ε=0/1へ戻せる。
全sourceのinverseはν/q/s非零chartでの逆diagonalと逆shearで、s零は別に保持する。
全modelは Z=diag(−2,1,1)+εB、k=−p−4。
全ρのD(ρ)はeigenvalues−2,1,1、(3)の全failuresから

    p∈{−4,−7/2,−2,−5,−8}, ε=0/1.

十fixed P0を残す。p=−4のk零境界も含める。

### h=b³, δ=6

    Z=(A0*a+B0*b+C0*c,0,D0*a+E0*b−A0*c).

全条件は

    γ=−qB0!=0, C0=0,
    A0=−(p+8)B0/q, k=−(p+8)²B0/(3q).

q/B0 unitを使うactual diagonal/relative normal/tensor scalesでq=B0=1。
normal shearがD0を動かし、base c-shearでE0も消せる。
p+8=0ではB0 unitのbase shearを使い、p+8で割らない。
全modelは

    Z=(-(p+8)a+b)∂a+(p+8)c∂c,
    k=−(p+8)²/3, γ=−1.

D(ρ)はeigenvalues0,±(p+8)。p+8!=0上の全(3) failuresとp+8=0境界を合わせ

    p∈{−14,−11,−19/2,−8,−13/2,−5,−2}

の七fixed P0を残す。

## 5. 実normalizerと逆の明示

全base diagonalはold(a,b,c)=(αa,βb,χc)、tensor scalar σ=χ/αと取る。
P_fをexactに保ち、X/Yのpとδは固定、qは(β/α)qへ移る。
relative normal v'=e*vはZ,kをe倍、κuvのκは変えない。
これらは非零定数のlinear mapsで、逆は各定数の逆を取るだけ。

normal shear v'=v+τuはZ→Z+τδA、k→k−κτ。
base c'=c+r a+s bもdet一で、X/Y radical時にはX/Yを変えない。
全higherはこれらの実PGL mapsで同時にpushし、その後必要なら忠実Euler gaugeへ戻す。

* h=a³は最初e=1/A0でA0=1とする。
  B0unit上で c'=c+r a+s b はD0→D0+B0r+s、E0→E0+2B0s。
  s=−E0/(2B0)、r=−(D0+s)/B0で消す。
  α=χ=1、β=1/B0、e=1/B0はA0=B0=1を同時に保つ。
  B0=0ではs又はnormal shearでD0だけを消す。
  q=0,E0unit上はα=χ=1、β²=1/E0、e=βでE0=1へ移し、両square-root choiceを許す。
  E0=0はそのまま残す。
* h=f,νunitはc'=c+(s/(3ν))bでoff Bを消す。
  runit上はe=1/r、α=β=1、χ=r/νでZ=H+A。
  r=0ではe=1、χ=1/νでZ=H。
  両者ともkはk/νへ移りpを変えない。
* h=ab²,νqunitはc'=c+(r/(3ν))aでoff Aを消す。
  α=1、β=1/q。sunit上はe=q/s、χ=s/(qν)、s=0上はe=1、χ=1/ν。
  これでq=ν=1、s=1/0、k=−p−4となる。
* h=b³は最初からqB0unit。
  c'=c+r a+s bはD0→D0+2A0r、E0→E0+A0s+B0r。
  A0unitならsでE0を消し、normal shear δ=6でD0を消す。
  A0=0ならr=−E0/B0でE0を消し、同じnormal shearでD0を消す。
  α=χ=1、β=1/q、e=q/B0によりq=B0=1となる。

どのzero branchも分母で落としていない。表に残る各全p-lineと全fixed値のliteral inputは
`inputs/plane_double_line_m2_constant_double_thirtyfour_symbols_20261108.json` に保存する。

## 6. 次の正確なwhole義務

空又は既essential三へ送った枝を除き、必要sourceは

* h=fの十wholep-lines（p!=k）。
* h=a³,ε1の一wholep-line（p!=4）。
* h=a³六＋h=ab²十＋h=b³七の計二十三fixed P0。

全bad-C jump、k零、nil/semisimple、off-field零/nonzeroを上の別枝で保持した。
これらは必要全sourceであり、未検査候補を全Poisson存在又は全G閉鎖とは呼ばない。
後続計算を行うなら原全P1三十/P2九、affine P2 particular、全p壁、actual normalizer逆と像への含意を先に合意する。
新しいprobe、jet又はminor列をこの手証明の代替として始めない。
