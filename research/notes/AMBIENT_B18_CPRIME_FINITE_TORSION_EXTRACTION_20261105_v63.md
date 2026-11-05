# Cprime の真正 ell-torsion を有限商から取る正確な条件

2026-11-05。これは DVR 加群の純粋な補題で、edge の値を計算したものではない。D=k[[ell]]、C は有限 D 加群とし、Cn=C/ell^n C と置く。

    Gn = ker(ell:Cn→Cn) / ell^(n−1) Cn

には自然写像 C[ell]→Gn がある。DVR 構造定理で C を自由部分と D/(ell^e) の直和に分ける。自由部分では Cn の socle と ell^(n−1) Cn が等しく、Gn には寄与しない。e<n の torsion summand では分母が零で、元の socle が一個残る。e≥n の summand では Cn の socle が分母に入って消える。この計算から自然写像は全射で、kernel は指数 e≥n の summand の socle の直和になる。kernel は標準的には C[ell]∩ell^(n−1)C である。証明の分解は非標準的だが、この写像と商自体は標準的である。

特に、全 torsion 指数が n より小さければ C[ell]≅Gn である。単に ker ell(Cn) を真正 torsion とするのは誤りで、分母を除いても n 以下という仮定だけでは不十分である。

実 Cprime は既証の rank二・ell零商長四・ell二乗商長七から

    Cprime ≅ D² ⊕ D/(ell) ⊕ D/(ell^e),  e≥2

という D 加群型を持つ。従って G2 は無条件で一次元となる。真正 Cprime[ell] の二次元のうち、高い指数の socle line は G2 に入る際に消える。genuine connecting Tnil→Cprime[ell] と合成した Tnil→G2 は全射・rank一で、kernel は十五次元。これは N/ellN の十四次元 kernel を含むが、それと等しくはない。

また

    length(Cprime/ell³Cprime) = 7+min(3,e)

なので、長九と e=2 は同値であり、e≥3 なら長十である。**e=2 が証明された場合に限り** G3=ker ell(Cprime/ell³)/ell²(Cprime/ell³) が真正二次元 target を与える。

現在の e は未決で、上記の十五次元 kernel の実 basis や edge 値も未評価である。formal numerator V_full σ を ell で割る edge の評価では、必要な ell 次数を一段上げた元の formal source を先に保持する。有限 Artin 商内で ell を反転しない。target の m二乗零から、その Cprime への埋込みの filtration strictness を推論しない。

依存は受理済み Cprime の generic rank二、ell零商長四、ell二乗商長七、genuine connecting16 の正確列である。新しい source/jet/gauge/積は用いない。この補題の独立監査は次の引継ぎ課題とする。
