# actual reduced normalizer tensor への完全比較 operator

2026-11-04。E4 は実零weight normalizer tensor、D0 は元 actual CA∩CB の
零weight subring である。両者を同一視しない。新しい係数展開、source/jet/点、
GB は行わず、受理済み proper↔physical unit interface と旧640 chain lift を
使って、E4_red への小さい**実係数付き有限 operator**を構成する。

## 1. reduced target と九つの実関係

受理済み nilideal の全 radical 証明により、同じ正則九変数環 C で

    Ered=E4_red=C/[d*L],
    L=(u1,...,u5,h,c*r,c*l,r*l),   d=c+r+l。

ここで c,r,l,h は actual core (a,d,x,t) からの有理 GL4 座標であり、u5 は
受理済み第二段五 pivot minus exact core root。全五つの d*u_j を残す。
C は同じ actual physical R0 上 free4 の代数で、basis は
(1,tauA,tauB,tauA*tauB)、各 trace root は自分の monic を満たす。

旧八正則関係 G8=(u5,Q0,Q1,f) と新九関係 R9=(u5,h,cr,cl,rl) について

    G8=D8x9*R9

という **実小多項式行列 D** を保存・全八恒等式を検査した。D の後三行は
nilideal packet に保存した Q0,Q1,f の I3-cofactor そのものであり、選択した
quad minor の逆を新たに使う必要もない。

## 2. 完全35 syzygy と全自由解像度

Z6=(u5,h) の Koszul complex と、(cr,cl,rl) の Hilbert–Burch complex を tensor
する。後者の二 syzygy は (l,-r,0) と (l,0,-c)。Z6 は三軸 quotient 上の正則列
なので tensor complex は exact である。自由 rank は

    (1,9,35,77,105,91,49,15,2)。

全 basis label・全 sparse differential を保存し、全次数の chain 積を QQ の
九変数多項式として厳密に確認した。最初の35列は、Z6 の Koszul15、Z6×quad18、
Hilbert–Burch2 である。d は正則環 C の非零因子だから Syz(dR9)=Syz(R9)。
最初の differential だけ d 倍すれば Ered の全自由解像度となる。
全 differential は C の maximal ideal に入るので、これは C 上で minimal。
R0 への free4 restriction は exact で、rank は四倍となるが R0-minimal とは
主張しない。特に新関係 Qrad4x36 の全 syzygy は36x140で与えられる。

## 3. 元 physical12 関係への literal adapter と transpose

受理済みの同じ unit source matrix N12 は

    N12*F_old12=(dG8,zero4)^t

を満たす。従って E8 を最初八行への injection とすると

    A12x9=N12^-1*E8*D8x9,
    F_old12=A12x9*(dR9)。

新しい operator は A の全108 C-entry の有限式を保存する。N12^-1 は旧
proper↔physical source unit system の同じ inverse であり、新しい任意係数行列
ではない。全 physical48 relation columns と同じ trace basis に対して

    Qold4x48=Qrad4x36*Res4(A^t),
    block(j,i)=mu(A[i,j])。

transpose は native matrix の外側にだけ掛ける。mu の内部を transpose しない。
これらの式は実八恒等式、受理済み N12 source identity、matrix inverse 恒等式の
合成から従い、一般点評価・係数近似・新しい cofactor 探索ではない。

各新 C-value は、actual proper u_i、GL4(kappa_core)、同じ N12 unit inverse、
有限和積から一意に定まる。四 R0 成分への lowering は受理済み16構造積と
unit norm-Hensel 仕様をそのまま使う。source96/group120未知数の branch と
非零 odd origin も保持する。未指定の正則函数や0/0商を新しく導入しない。

## 4. 全640旧関係の実 chain lift を新36列へ運ぶ

旧 actual source-bound B48x640 は UP=Qold*B を満たす。旧 r/g 函数は元の
regular-function root のまま opaque に保つ。新しい coefficient との composition
で分母を評価し直さない。operator は旧Bを十二の C-value block にまとめ、

    Brad=Res4(A^t)*B

の全640列を literal finite C-value circuit として構成し、各成分を読み出す。
したがって UP=Qrad*Brad は actual source map の恒等式となる。
旧B自体の構築を再実行したのではなく、既受理全640列を新adapterで運んだ。

同時に Qrad4x36 と全 Srad36x140 を同じ C/R0 basis で保存する。
新35-column complex を exact restriction したため、関係の一部分だけを
切り出した kernel criterion ではない。

## 5. 原点 Nil(D0) の exact homology model

P9x640 は受理済み honest C0 presentation、U4x9 はその同じ actual comparison。

    V=(Qrad,-U):R0^45 -> R0^4,
    X=((Brad,P),(Srad,zero9)):R0^780 -> R0^45。

全 chain 恒等式は前節と完全 syzygy から従う。H1 が ker(D0->Ered) に等しい
ことも exact：kernel pair の C0 class を取り、その class が source relation なら
P を引き、残った Qrad-syzygy を全 Srad で消せる。逆写像は Qrad の持上げを使う。

U の第一列は実際に basis1そのもの。-U0=-e0 の単位 target/source cancellation
では、他の44 source列の row0を unit列で消す。chain恒等式により対応するXの
pivot行は零となり、他行は変わらない。従って正確に

    R0^780 --Xred--> R0^44 --Vred--> R0^3

を得る。Xred の最初640列は (Brad36,P8)、後140列は (Srad36,zero8)。
Vred の全3x44 entries、Brad の全36x640 entries、P8の元函数ID、Sradの全列を
新しい finite operator に保存する。

元 D0->E4 は finite integral、Spec 上 surjective。その kernel は nil であり、
Ered は被約である。よって

    Nil(D0)=ker(D0->Ered)=H1(Xred,Vred)

が厳密に成立する。これは**完全な比較 homology への表示**であって、まだ
H1 の generator を取り出した cokernel presentation ではない。
Nil(E4) の全生成元が D0 から持ち上がるとは仮定せず、Nil(D0) が principal とも
仮定しない。Nil(D0)^3、socle、元の深さや finite T の長さを計算したものでもない。

新operatorは実D8x9、全native complex、全A12x9 entries、全Brad640列と完全Sradを
追加するものであり、既存の48列表示の単なる shape 言換えではない。
新 scalar node≤150000、polynomial terms≤150000、output≤32MiBを強制する。
展開された物理係数冪級数や新 source solve は生成しない。
