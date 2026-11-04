# B-rankbadを含む有限free次数112のcurve cover

norm出力を既知単元だけで整理すると、厳密に
R=content(R) * a³¹ C⁴ * hprim、hprim=h/10⁸となる。
これは31個の保存gcd因子が最初の4回aC/100、残り27回aであることと、
全R係数の再構成で確認した恒等式であり、未知因子のfactorizationではない。

hは225項、a次数14、e次数28で、最高e係数は非零のQQ定数。
E=h/LCe(h)とするとEはeのモニック28次式になる。既存Hはdのモニック4次式。
R0=QQ[a,a⁻¹]上で

R0[e,d]/(E,H)

は基底e^i d^j（0≤i<28, 0≤j<4）を持つ有限free rank112代数である。
まずEでR0[e]を割り、その商の上でモニックHを使うため、E/Hの既約性や
reduced性は不要。これはpair quotientがfree112との主張でもない。

実際のFULLU B-rankbad点では選択minor f=g=0である。既存のa-unit membership
からH=0。さらにR=0かつaCが単元なのでh=0、従ってE=0。
よって全てのB-bad点はこの有限curve coverに含まれる。B-badの実在・空性や
全域rank3は未証明であり、cover全体の点がB-badであるとは主張しない。

この代数にはまだFULLUのdFBCJの逆元を加えていない。追加局所化でbase上の
有限性が保たれるとは自動的に言えない。special fibreによるNakayamaを使うには、
全局所化因子の単元性又はexactなclopen分離を別に検証する必要がある。
