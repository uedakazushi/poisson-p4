# 全三次行列の一つの最大 minor：水平次数と既知 B 因子

対象は `S61_BINARY_HYPERPLANE_CUBIC_MATRIX_20261102` の元 135 式から作る
全 weighted cubic 行列である。既知単元による 21 pivot 後の 42×70 行列から、
保存された 36 行・36 列を取った行列式を P(a,e,d) とする。
本証明は **deg_d P=79** と **B^12 が P を割ること**を扱う。
全 core の階数被覆、点ごとの階数、S61 の affine exhaustion はここでは主張しない。

既知因子は

- F=500a+90e²+3e−24,
- B=(50a+18e²−9e−2)d+60e−40,
- C=100a+18e²−9e−2

であり、係数体を K=QQ(a,e)、素数を p=32003 とする。

## d 無限遠の厳密次数

36×36 の各成分の d 次数について Hungarian 双対整数 r_i,c_j を保存する。
全成分で r_i+c_j≥deg_d A_ij、総和は 107 である。
従って M_ij(δ)=δ^(r_i+c_j)A_ij(1/δ) は多項式で、
`det M=δ^107 P(a,e,1/δ)` が成立する。個々の双対整数が負であっても、
各成分の δ 指数は非負である。

M(0) の QQ[a,e] 上の核を厳密に求め、pivot 座標列を補完した定数基底変換を行う。
核の free block は非零多項式の対角行列であり、補完の行列式は K の単元である。
核列を掛けた全成分を QQ[a,e,δ] で計算し、δ による割算をすべて検査する。
この手順の nullity は順に 14,9,5 で、最終 leading 行列の階数は 36 となる。
最後の行列式は (a,e)=(2,3) mod p で 21773 である。
従って δ 消滅次数は正確に 28、P の d 次数は 107−28=79 である。

さらに、読み込んだ有理係数の分母はすべて p と素であり、各核の free diagonal
多項式は mod p で非零であることを別に検査する。
V=Z_(p)[a,e]_(p) 上で全基底変換が可逆で、最終 leading determinant も V の単元となる。
これは単なる有限体標本の次数観察ではない。

## B に沿う消滅次数

QQ[a,e,d] の primitive Schur を 18×18 まで実行する。
各 row/column content の割算と各 Schur 積は厳密で、signed factor bag を保存する。
P と残 18×18 行列 A の行列式との有理関数関係において、bag の B 消滅次数の
総和は正確に 2 である。負指数を捨てず、全 18 factor の B 除算を検査する。

h=50a+18e²−9e−2、d=(δ−60e+40)/h と置く。
各行に十分な h 冪を掛けて多項式化する。h は K の単元で、δ=B である。
ここでは最初の一階 jet M0+δM1 だけを保持する。高次項は、以下の最初の
leading 核と一回の δ 割算後の定数項に影響しない。

M0 の階数は正確に 11 で、右核 K0 と左核 L0 は各 7 本を厳密に求める。
右核の free diagonal と選択した 11 pivot 列から補完基底を作ると、
7 核列は δ で割れる。一回割った後の leading 行列は

`[ M0 E , M1 K0 ]`

である。M0 E の階数は 11、L0 はその余核全体を与えるので、その階数は
`11+rank(L0 M1 K0)` に等しい。

この 7×7 行列を元の leading データから直接積算する。非零 row/column content
と非零有理定数の除算は K 上の階数を保存する。
簡約した 7×7 行列の 3 本の独立な QQ[a,e] 核を、全成分の積が零となることまで検査する。
また mod p の非零 4 minor があるため、階数は正確に 4 である。
従って一回割った行列の leading corank は 3 となる。
残 determinant はさらに δ³ で割れ、元 18×18 determinant は δ^(7+3) で割れる。
P では bag の 2 を加えて **B^12 | P** を得る。
最終の巨大行列を再び展開する必要はなく、leading corank の補題だけを用いる。

有限体の別証明 `S61_BINARY_HYPERPLANE_CUBIC_MODULAR_20261102` は P mod p の
B 重複度を正確に 12 とする。B の d-leading coefficient h は V の単元なので、
その結果と組み合わせれば QQ での重複度も 12 となる。
本 report 自体には、独立な QQ 側の必要条件である「少なくとも 12」を保存する。

## 再現と限界

主 producer は元 135 式を再構成し、保存された first-minor index と Schur matrix hash を照合する。
巨大な work cache は入力としない。生成物は TemporaryDirectory 内に置き、
全 source/input hash、kernel/rank metadata、生成物 hash、native log と実行時間を report に保存する。
`--output` は必須で、既存 report を上書きしない。

探索時の 36×36 B-kernel と次の 18×18 kernel の時間切れは、完成した証明には使用しない。
最後の証明は左右核から得る 7×7 階数に置き換えた。
未知の core 曲線を除外するには、この結果に加え、別の有限体 gcd 証明と垂直係数の
整除・先頭係数検査が必要である。
