# 被約非concurrent零cubic・位数三 root covers の独立幾何確認

2026-11-08、PRIMARY G、別作者 `/root/quartic_review`。
開始時 HEAD は `136a2a3ccb459da913b5db2e274bb5389830b9a6`。
本稿は二つの凍結 producer 証明を独立に読解・再導出した純証明 receipt である。
新 CAS、source/operator/minor 計算、jet、anchor は実行していない。

* `PLANE_CONIC_LINE_M3_ROOT_COVERS_AND_EXCLUSION.md`
  (`9ab911206e4606d8c79bf1466aa4e8179e38011ea58f6fe463713069f46ad9b6`)
* `PLANE_TRIANGLE_M3_ROOT_COVERS_AND_FIRST_SYMBOLS.md`
  (`81cb35e2d159ede8f024e385d0c688af3609954b7195464277b43d10b4a79a91`)

**独立受理する結論。** rank 四 Poisson、被約既約 actual quintic、restriction rank 二 plane、
normal 位数ちょうど三の下で、横断 conic+line の全 source は空である。
接線 conic+line は S/N/D/M の四 first-connection profiles、非concurrent triangle は
C3-unused-line と四 triple の五 diagonal profiles に、全 source を必要還元できる。
残九 profiles の higher source の次元・存在・成分帰属は本 receipt の結論ではない。

## 1. 全 cover と、例外曲線を失わない次数制限

既に別作者確認された associated-graded parallel tensor と全 root-degree 入口を用いる。
squarefree の root incidence は U=P2\R 上で finite étale、repeated 21/3 は別補集合である。
横断 conic の U は {ab!=1}、a!=0 patch は (C*)²。
近傍 {|ab|<1/2} は contractible、その patch との交わりは C*×disc なので van Kampen は
a-meridian だけを殺す。従って π1=Z。
接線 conic は C×C*、triangle は (C*)²である。
従って全 SF permutation image は 1/C2/C3。S3 を使わない理由はこの abelian 性であり、
concurrent の S3 cover を排除したのではない。

primitive な実 incidence の class を (e,3)、e<=2 とする。
splitting cover の正規 compactification を取り、deck-equivariant に三 root maps を同時解消する。
p:T→P2、H=p*O(1)、N=deg p とし、各 root pair は L_i=phi_i*O(1) の生成 section pair とする。
各 DVR で pair の一成分は unit である。Gauss content lemma より三 pair の積にも unit coefficient がある。
従って pullback 実 G3 との比 h は全 prime divisor で pole を持たず、normality より正則であり、

    eH−L1−L2−L3 = div(h) >= 0,
    Σ H.L_i <= e H² = eN <= 2N.

これは単なる一般 root の次数推測ではない。ramification、元の content、全 exceptional divisors は
effective な右辺に保持される。connected C3 では三 H-degrees は deck 対称で等しい。
非定数なので正、従って各 degree は一又は二。

minimal resolution Y へ root pencil を push し固定曲線を除くと、movable Cartier class D を得る。
共通解消の exceptional curves は H-degree 零なので HD は変わらない。
push した二一般 members に共通曲線は無く、D は nef、D²>=0。
conic+line の cyclic cubic surfaces と triangle の w³=xyz は全て Du Val cubic で、
H=−KY、H²=3、χ(OY)=1。
Hodge と Cartier parity より

    0<=D²<=(HD)²/3,    D²−HD ∈ 2Z.

HD=1 は不可能、HD=2 では D²=0。
従って basepoint は無い。Stein factor の degree が二なら primitive fiber の H-degree 一となって矛盾する。
全 root maps は connected conic pencils F、HF=2、F²=0 で尽くされる。
解消前の basepoint、特異点上の exceptional fixed parts、normalization の非flat性を仮定で消していない。

E=H−F は E²=−1、HE=1、χ(E)=1、H(K−E)=−4 より effective。
push E は projective line ℓ と H-degree 零の exceptional parts である。
h0(F)=2 で、E の section を掛けた H0(F) は H0(H) 内の ℓ を含む hyperplanes の二次元空間全体に等しい。
よって root pencil は Σ の実 line からの hyperplane pencil の比である。
line の一般族だけを使う主張ではなく、全 line の列挙へ移れる。
使用した標準幾何入力は Du Val の crepant/rational resolution、surface RR、Hodge、Bertini/Stein のみ。

## 2. conic+line の全 C2/C3 補集合

C2 の rational fixed root の primitive degree d、residual quadratic degree e は d+e<=2。
d=2 又は e=0 は split、nontrivial C2 は d=0/1 のみ。
d=0 では G=u(Au²+2Buv+Cv²)、C,T=AC−B² は非零 R-supported。
全 higher source の constant-root 障害 q(Cbad)=0、Cbad!=0 から q=dlogL、従って C=Q。
T の全候補 L4/L²Q/Q² のうち nonsquare は L²Q のみ。
Q|B、constant normal shear により B=0 となり D profile を得る。
content を持つ残余も C=Q と nonsquare 条件が検査しており、零係数を SF と誤同定しない。

d=1 では F の discriminant は Q、三 linear coefficients は independent。
T=(F=0) は P1×P1、H=(1,1)、normal pencil class F0=(1,0)、ell|T の class は (2,1)。
resultant の R-supported 性から ell divisor の全 support は L/Q の逆像だけである。
横断 L の逆像と ramification は共に (1,1) なので (2,1) は作れない。
接線 L の逆像 A+B と ramification Qr に対しては 2A+B 又は A+Qr のみ。
前者 ell=c*u は nonprimitive、後者は唯一 ell=b*u+c*v。
restriction の kernel H0(O(0,−1))=0 なので同じ divisor の別 section を失わない。
これが M の全非対数接続であり、有理 frame (b*u+c*v,u) を tensor の PGL 同値には用いていない。

C3 cover は w³=LQ。横断型は 2A2、接線型は E6 のみである。
前節の line 列挙では projection center は Σ 外なので各 line は base line M へ同型に写る。
w³|M が linear cube となる条件より、ramification line L=w=0 と、横断の二交点における conic tangent
各三 lifts が全て。接線の場合は ramification line 一本だけ。
その pencil は w,L で、norm の content L を除く全 primitive candidate は Q*u³+L²*v³。
横断の extra line の代表は y=0,w+z=0、その primitive norm は

    xz*u³+3z²*u²v+3zy*uv²+y²*v³,
    discriminant = −27*y²*z²*(xy−z²)².

y は R 外なので全 extra 六 lines は SF parallel 性に反する。
同じ root maps を持つ degree 二 polynomial は primitive norm の constant 倍のみ。
残 candidate の一 pure coefficient Q は Cbad(f) の像に無い。
実際 Cbad(LQ)=kQ なら modulo Q で Cbad(Q)=0、Q nondegenerate から Cbad=0、k=0。
従って全 C3 source は空。

S/N/D に残る constant roots は全 curl (2k11 I+4DY−DX)Cbad=0 を満たす。
横断では k11=0、ker DV=C∂c、しかし q=dlogc は Cbad(c)=0 を要求するので全空。
接線では ker DV=C∂a で、この矛盾を外挿しない。
この curl は full P1 三十/P2 九の literal identity であり、k11!=0、m二 SF 又は log 性を使わない。

接線の全 regular parallel spaces はそれぞれ Q Sym³、Q*u³、span(L²*u³,Q*u*v²)、
(bu+Lv)(au²+2buv+Lv²) の一 line で正しい。
M では frame z=bu+Lv,w=u 上の有理 terms (A*z³+B*Q*z*w²)/L の polynomial 性が A=B を強制する。
N の零 nilpotent は S、D の両 amplitude の零境界も retained。全 triple/21/split SF は既 degree 入口で含まれる。

## 3. triangle の全二次・三次 covers

C2 branch は二 lines、C3 branch characters は (1,1,1) 又は (1,2,0) とその inverse/permutations で尽くす。
C2 moving degree 一では F=xu²+yv²。incidence surface は F2 で、H=C0+2f、ell class C0+3f。
vertex exceptional C0、二 branch fibers Fx/Fy、第三 line の section Z~H を全て含めると、
ell divisor は C0+bFx+cFy (b+c=3) 又は Z+Fx/Z+Fy のみ。
前者の ell=r(βu³+δu²v−αuv²−γv³) は支持条件から一 monomial となり nonprimitive、
後者は z*u/z*v で同じく nonprimitive。全 moving C2 は空。
constant fixed root C2 では C=xy/xz/yz、T=xyz²/x³y/xy³ の全候補を保持する。
C=xy の平方 coefficient は actualD3 pure-image span(xy,xz,yz) 外。
C=xz/yz に残る profile は curl diagonal (9/2,3/2,−3) が invertible で空。

C3 all-three cover w³=xyz の全 lines は三 ramification lines のみ。
norm は xy*u³+z²*v³、square coefficient が同じ actualD3 image 外なので空。
C3 two-line cover の正規化は P(1,1,3)、resolution は F3、
H=C0+3f、HD=b、movable D=aC0+bf は b>=3a。
HD<=2 より a=0、D=f 又は 2f。全 root maps は φ(s/t)、degree 一又は二である。
degree 二の deck involution τ と φ(r)=φ(ζr) の collision を比較すると、
tracefree τ=[[a,b],[c,−a]] の polynomial ζc r²−a(1+ζ)r−b は 0/∞ 以外に root を持てない。
invertibility を保持すると b=c=0、τ(r)=−r のみ。
従って φ は normal PGL2 を除き r²、leading は x²*u³+y²*v³ で空。
degree 一では h=x/y が同様に空、unused z の全 span(xz*u³,yz*v³) だけ残る。
normal basis changes は constant、rational cover frames を actual PGL と混同していない。

## 4. triangle の全 repeated branches と非対数部分

既 constant-root 障害後の split 000/001 は、記載された curl が invertible 又は bad-axis 条件と両立せず空。
21 の全 02/10 は既 active-component 障害で空。
00 の全 period q では bad axis と curl が q=(0,1/2,1/2) の三 permutations を要求する。
他方の constant root には全三値で bad vector が存在せず空。
01 の全 marked residues は H=L の ρ∈{−1,2}、H!=L の ρ∈{0,2} を列挙し、
記載された四 curl を代入して全空である。

triple では全 quadratic monomials h を保持する。
q_i=0 axis の curl は (3/2)(3−h_i)(q_j−q_k)、3−h_i!=0 より再び q=(0,1/2,1/2) permutations。
off-diagonal は最初から log とはしない。全 flatness χ(Z)=0、χ_i=(q_i+1+h_i)/2!=0 より
Z 各 row は corresponding coordinate 倍、従って diagonal。残 r は closed logarithmic になり、
二つの horizontal closed log forms の wedge と p−q!=0 により r=constant*(p−q)。
constant unipotent inverse で消せる。
q の零 axis を x とすると、残 y/z 交換を除く h=x²,y²,yz,xy が全四型。
対応する p と全五 first-symbol JSON の X/Y/K11 は character 式と整合する。
firstsymbol は toric logcanonical strict Poisson、D2=0 を満たす。
固定 plane 上の source upper B=triangle choices6+anchor amplitude1+normal conjugacy2=9 も正しい。
leading amplitudes は higher source に残し、今後の fiber upper がなければ全 G 上界へ進めない。

## 5. 受理範囲と残余

本 receipt は whole necessary reduction と横断空性の純証明である。
接線四/triangle五 profiles の global K consistency、full higher compatibility、actual quintic の存在、
G 次元又は E membership は後続の別有限仕事である。
concurrent、nonreduced R、m二/m四、他 singular surfaces を含めた全分類は主張しない。
gate と入力 SHA 照合は custody/control-plane の確認であり、数学の機械証明とは呼ばない。
