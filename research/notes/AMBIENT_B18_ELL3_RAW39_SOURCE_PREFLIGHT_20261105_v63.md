# ell三乗の実 product base39：三次 source だけで正確に決まる

B=R0/(IA IB)、U=R0/(IA∩IB)、T=(IA∩IB)/(IA IB) とする。実両 conductor 曲線は D=QQ[[ell]] 上 **free3** であり、free2 ではない。従って U は free6。T は一般に J=IA+IB に消され、実 J=(z7,r²,rs,s²) は ell² を含む。

既存の実長さ B/ell=20、B/ell²=33 から、T の最小生成数は14、長さは21、型は(D/ell)^7⊕(D/ell²)^7。よって B/ell³ の長さは39、raw cokernel の長さは9+9+39=57である。これは right-exact 商の基底のサイズであり、C 自身の torsion 指数の上界とはしない。

union の既存 Weierstrass basis は (1,1),(r,r),(r²,0),(0,r²),(ell²,0),(ell r,0)。重みを0,1,2,2,2,2と付けると、r³=A2 r²+A1 r+A0 における ell 次数1,2,3の評価から r と ell は重みを少なくとも1上げ、zi は少なくとも2上げる。従って m^5 U⊂ell³U である。未知の四次係数は使わない。

0→T→B3→U3→0 と B3→B2 を考える。T は ell² に消されるので reduction は T 上単射である。m^5 B3 は U3 で0となるので T に入り、既存の actual raw33 では m^5 B2=0 だから reduction の核にも入る。従って **m^5 B3=0**。

同じ論法は三次 Taylor で切った conductor pair にも成立する。その ell² 商は既存の厳密な cubic-truncation theorem によって actual raw33 と等しい。両者で m^5=0 なので、四次以上の generator の省略による product の差は0となり、

(IA IB,ell³)=(IA_≤3 IB_≤3,ell³)

が全次数で正確に従う。beta3 の z 依存を既存の unit generator change で除いたときに出る四次 plane 項も同様に許される。

49本の monic zi(zj+beta_j) を消すと、1,z7 の plane module、21 pair と14 g の関係になる。m^5=0、ell³=0 の範囲で ambient は75、必要な plane 多項式倍は高々77行・938項である。この SOURCE packet は実入力を保存するが、rank36の計算や39基底の乗法表はまだ計算しない。

注意すべき点として、個々の branch quotient R0/IA とR0/IB moduloell³ には四次項が寄与し得る。product base39 の cubic-truncation から、branch の四次係数も0としてよいとは従わない。後の有限105 source の厳密評価で actual rawbA/bB ideals を求め、各9次元の branch quotient を作る方法ならこの問題を避けられる。この次の source値評価はまだ実行していない。
