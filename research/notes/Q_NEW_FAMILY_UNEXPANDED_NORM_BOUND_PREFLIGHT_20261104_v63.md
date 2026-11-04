# 有限10次代数の未展開norm：非零な方程式倍と次数上界9840のpreflight

2026-11-04。四つの∞poleやgeneric交点数20だけから定義次数を決めない。
ここでは既受理のG10／primitive C1／balanced U,Vを使い、**展開しない有限operator**として
非零なinvariant polynomial倍を定義し、保存supportだけで粗い上界を得る。
新しいnorm、resultant、polynomial積、点、prime、source高次数jobは実行しない。

K0=QQ(R,Z0,Z1,Z2)、L=K0[r]/(monic G10)、L0=lc_r(G)。
primitive C1=A_lin*s+B_constからs=−B_const/A_lin。
既存のbalanceによりactual eta=N12/F4³はU/Vで、U,Vはrと四base変数の多項式。
V≠0は既受理のprimitive chartから従う。したがって

    P(eta)=L0^134 det_10(eta*M_V−M_U)

は非零。これはfixed degree(10,134)のSylvester determinantのroot-product式なので、
leading paddingを含めてbase変数の多項式である。eta次数10の係数はL0^134 Norm(V)≠0。
actual sourceではeta=U/VゆえPは零となる。

support次数を `(r, R−Z0, 2(Z0+Z1+Z2))` で測る。保存全係数を走査すると

    G10: (10,0,10),  C1の二係数: (31,8,10)。

小さいJ,Fn,N79,N41のs-homogenizationから、U,Vのr次数は≤134、R−Z0次数は≤34。
etaへweight3を与えたとき `eta*V−U` の第三次数は≤44。
resultantの係数同次性によりPでは

    exponent_R−exponent_Z0 ≤ 134*0+10*34 = 340,
    2Σexponent_Z+3exponent_eta ≤ 134*10+10*44 = 1780。

R=(p3+4p4)/p0、Zi=pi/F4²、eta=N12/F4³を代入し、p0^340*F4^1780を掛ける。
各monomialのp0指数340+b−a、F4指数1780−2(b+c+d)−3eは非負なので多項式となる。
p0=F4²/400−4p4も代入すれば、六つの選択low invariantの多項式Hを得る。
そのβ-homogeneous次数は

    8*340+4*1780=9840。

選択low6はambient上で代数的独立で、このratio変換はF4を残してbirational。
従ってHはambient W16上で非零。actual NEW hypersurface上で零なので、既約定義式f_NEWは
Hを割り、定義次数は **9840以下** と分かる。これは具体的な未展開倍であって最小方程式ではない。

low6はHSOPでも全域有限でもなく、一般cover次数8のため、HにはNEW以外のhypersurface因子や
境界因子が含まれ得る。QQ因数分解とdense actual sourceへの恒等式代入で因子を分離する
有限strategyはあるが、今回その展開・因数分解は実行していない。
generic20交点と∞poleからの根拠のない次数等式、最小多項式の同定は主張しない。
