# 全 weighted cubic 最大 minor の有限体 gcd

この記録は p=32003 の厳密な多項式行列計算である。
QQ への持ち上げ、全点での階数、S61 の exhaustion は単独では主張しない。
元 135 式から再構成した全 63×91 三次行列に、既知 a,F,B,C 単元の 21 pivot を施した
42×70 行列を使う。後の pure-kernel 生成元だけの行列とは区別する。

保存された first 36 minor と、行 index1 を19、列 index22 を35 に置き換えた隣接 minor を計算する。
有限体 primitive Schur は全 content 割算を厳密に検査し、各 pivot の指数 −(n−2) を
signed factor bag に記録する。反対符号の因子は gcd と exact division で相殺する。
最後の全指数は非負である。有理定数の代わりに Fp* の scalar normalization を行うので、
行列式の非零定数倍は一つの保存 point determinant から復元する。

first minor の個別次数は (a,e,d)=(161,213,79)、全次数310であり、
lex(d,e,a) leading term は係数8492の a^101 e^129 d^79 である。
既知因子の重複度は、全 bag を繰り返し exact division して
**a^37 F^10 B^12 C^26** と得る。

両 bag から既知単元を除いた後の pairwise gcd 支持は、21項、92項、6項の三多項式だけになる。
それぞれ FLINT の完全 multivariate factorization で既約因子一つ・指数1となり、
factor expansion と元多項式の一致も検査する。
三つの既約因子について、それぞれ次の有理零点を保存する。

- (a,e,d)=(15781,20631,1240),
- (26448,2,3),
- (5990,2,3).

各点で因子値は零、aFBC は非零、全 42×70 行列は階数36である。
非零 36 minor の行・列 index と determinant 値を report に保存する。
従って三つの既約因子はいずれも全 maximal minor の共通因子ではない。
全 minor の共通既約因子は、計算した二つの minor を割るので、既知単元を除けば
必ずこの三つのいずれかでなければならない。
以上から **Fp[a,e,d,(aFBC)^−1] における全 maximal minor の gcd は単元**となる。
これは階数低下 locus が余次元2以上であることを意味し、階数低下点が無いという主張ではない。

主 producer は元 135 式から source matrix を fresh 再構成し、既存 work cache を使用しない。
有限体へ写す全有理係数について分母が p と素であることを確認する。
二つの小行列式、gcd、既約分解、三つの零点と full-rank witness は毎回再計算する。
全 native source、入力 hash、生成物 hash と実行 log を保存する。
QQ の core 被覆への橋は、別の水平次数・B 整除と垂直先頭係数・aFC 整除に依存する。
