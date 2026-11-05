# 実 Cprime の指数 (1,2) と全 Hilbert 関数

既存の同一 source/DVR について Cprime は D=QQ[[ell]] 上有限、rank 2 であり、ell 零商長 4、ell² 商長 7 から D 加群型は D²⊕D/(ell)⊕D/(ell^e), e≥2 であった。これは非標準の D 加群同型であり、R0 加群や代数の分解ではない。ell=5(sA−rA) で、rA と区別する。

完全二次作用素の独立検証により真の HF の最初の三項は (3,3,3) である。既知 cubic closure の 27 行中 26 行を消す固定汎関数に対し、実一次 syzygy 第 49 列の三次 initial V2 u1 は 99/64 と非零になった。全 135 係数の V1 u1=0 と、既存の二次値からの pairing は独立に検証済みである。従って既知関係、ell³、およびこの実三次関係によって全 27 行が消え、m³(Cprime/ell³)=0 となる。

よって length(Cprime/ell³)≤3+3+3=9。一方、DVR 構造定理による同長は 7+min(3,e)≥9 である。等号が成立し、e=2 が確定する。実 Cprime の torsion 長は 3 であり、任意 n≥1 に対して ell^n 商長は 2n+min(n,1)+min(n,2) である。**これは Cprime の torsion であって、元 H0(D0) の長さではない。**

さらに m³(Cprime/ell³)=0 と自明な ell³Cprime⊂m³Cprime から m³Cprime=ell³Cprime を得る。

    ell³Cprime ⊂ ell m²Cprime ⊂ m³Cprime=ell³Cprime

なので m³Cprime=ell m²Cprime である。m を掛けて帰納すると j≥2 で m^{j+1}Cprime=ell m^jCprime、j≥3 で m^jCprime=ell^jCprime が従う。これらは D-free rank 2 である。従って実 Cprime の全 m-adic Hilbert 関数は

    (3,3,3,2,2,...)、 Hilbert series=(3−t³)/(1−t)、 multiplicity=2

となる。Cprime 全体で m³=0 と言っているのではない。

指数が 3 未満と確定したため、真正な二次元 target は標準的に

    Cprime[ell] ≅ ker(ell:Cprime/ell³→Cprime/ell³)/ell²(Cprime/ell³)

から得られる。有限商の kernel 全体には free 部分の top が含まれるため、それをそのまま真正 torsion としてはならない。Cprime/ell³=Cprime/m³ の実九次元作用を構成すること、および genuine connecting の値は、次の独立した source/value gate である。

この帰結は新 source・gauge・jet を用いない。元 nil 三乗、nil generators、元 H0(D0) 長、元 socle の結論は出さない。
