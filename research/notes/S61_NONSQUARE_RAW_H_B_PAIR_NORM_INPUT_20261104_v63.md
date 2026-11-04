# モニック4次モデルの4×4ノルム入力

前packetのモニックHを使い、S=QQ[a,a⁻¹,e][d]/(H)の基底1,d,d²,d³に
関するgの乗算行列Mを作る。全4列についてg*d^j=H*qj+Σi Mij*d^iを
全係数で検査する。これはモニック除算だけであり、未知因子は逆にしない。
行列式・終結式、GB、点・有限体評価はまだ計算しない。

(f,g)=(H,g)なので、R=detMはpair-zero集合の射影を含む非零多項式となる。
非零性は既存のexact gcd(f,g)=1とGauss補題による。実際の係数は未計算。
S/(g)はSからの商なので係数環上有限であり、両最高係数消失によるescapeはない。
ただしFULLUでのdFBCJの追加局所化後も有限とはまだ主張しない。

次にRを計算する場合はv=adj(M)e0を作り、v(d)g−RをモニックHで割って
exactなBezoutを保存できる。Hの既存linear combinationを代入すればf,gからの
Rのliteral membershipになる。Rの未知因子は単元とみなさない。U-supported
因子を除く場合にもexact gcdと再構成が別途必要である。
