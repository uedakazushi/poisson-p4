# normalized E4 nil 二乗の三支持：独立監査

対象は `E4=S9/(dJ)` の `M=Nil(E4)^2` のみ。元 D0 の nilideal やその像の二乗とは区別する。producer は実行・import せず、独立 Fraction 多項式演算と既監査の generic Artin source を用いる。

実 Q0,Q1,f を受理済み CI source と一致させ、五つの symbolic 直線を代入する。保存済み N/P の Laurent 模型について Q0,Q1 の二次までの零性と f の先頭係数を照合する。新しい source jet の解法や数値点は使わない。

N/P の nilideal は局所的に d*(S/J) である。乗法 d は module 同型であって環同型ではない。そのため nil 二乗は core の d 倍に対応する。N の core は kappa[v]/v³、d は v² の非零倍、P は kappa[v]/v²、d=v。従って M の generic 長さは両方で 1 である。C では d は単元で、受理済み HF=(1,2,1) より maximal ideal 二乗が長さ 1。R の双対数、L の被約点、主成分 (d) では二乗は零となる。

`Ass(M)⊂Ass(E4)` と既受理の全六 associated primes を使う。以上の局所評価から Ass(M) は正確に N,P,C。これらは全て一次元直線であり原点は associated でないから、M は一次元 CM module となる。特に H0(M)=0。これは E4 への部分 module 包含から H0(E4)=0 を用いても従う。

独立に `c=d,p=a/2,n=t−a/2,H=a−d+x` と逆変換を検証し、行列式 1/2 を確認する。三直線は n,p,c 軸となる。三つの monomial prime の交叉を LCM によって直接計算すると `(H,cp,cn,pn)`。原座標へ戻すと、五つの u と共に

`I=(u1,…,u5,a−d+x,ad,d(2t−a),a(2t−a))`

を得る。squarefree monomial prime の交叉なので reduced である。

各 generic localization の M は単純な長さ 1 module であるため、I*M はそこで零となる。もし I*M が非零なら、その associated prime は部分 module として Ass(M) のいずれかとなるが、三点の localization は全て零で矛盾する。従って I*M=0。逆に三点全てで M は非零なので Ann(M) は各 prime に含まれる。よって Ann(M)=I は正確である。各枝の generic rank は 1、直線の重複度も 1 なので associativity formula により e(M)=3。

結論は normalized M の `dim=depth=1, e=3` と reduced annihilator である。M が三軸環そのものや free/canonical module であるとは主張しない。D0 の像が N/P/C の三点全てで平方零となるかは別の source gate であり、元 nilideal の上界 8 は導かない。
