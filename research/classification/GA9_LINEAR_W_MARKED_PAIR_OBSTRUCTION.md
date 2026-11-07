# Ga9 linear-w と A3+A2 marked quartic の閉包の区別

2026-11-07、PRIMARY E。既済の全 Ga9 linear-w 源は一つの23次元射影軌道の閉包に入る。
その退化因子は、単純超平面 H と平面 E6 quartic F の組である。
本稿は、この点が既知103番目の24次元成分に入らないことを、marked pair の幾何から示す。
新たな Poisson 成分の存在・次元や E 被覆の完了は主張しない。
以下は今回手で確認した証明で、独立の数学査読は別ノートに記録する。
旧計算証明書・producer bytes は変更しない。

## 入力と命題

103番目の一般退化因子の既済記述は
`GA9_CUSPIDAL_QUARTIC_COMPONENT.md` と
`CLASSIFICATION_GA9_CUSPIDAL_QUARTIC_COMPONENT_20261107.json`
（SHA256 `4b9d73769906c97c1ab5d90c658750445d3688a7e499ba47cc2324aee305fe27`）による。
既約 quartic は A3 tacnode と A2 cusp を持つ有理平面曲線であり、
H は tacnode の二枝の共通接線で、それぞれの枝と order2 で交わる。
H と quartic の交点はこの一点だけである。

Ga9 linear-w の全源・軌道閉包の既済記述は
`HYPERPLANE_QUARTIC_NIL211_SHORT_TRANSLATION.md` §8 と
`NIL211_SHORT_TRANSLATION_REVIEW.md` による。
canonical 退化因子は非零定数倍を除いて

    D0=h F0,   F0=4h³w−8h²y²+y⁴.

F0 は点 p0=[0:0:1] だけに E6 cusp を持つ既約有理 quartic である。
正規化上、p0 の唯一の逆像で h,y,w の消失位数は順に4,3,0である。
H0={h=0} はその cusp tangent である。

**命題。** 上記 A3+A2 quartic と marked common tangent H の PGL 安定な族の閉包は、
marked pair (H0,F0) を含まない。
従って、103番目の全 Poisson 成分は Ga9 linear-w の密な23軌道を含まない。
これは unmarked quartic の E6 と A3+A2 の隣接関係を否定する命題ではない。

## 仮定した退化を DVR 上の平面族へ移す

反対に閉包に含まれると仮定する。
曲線選択と有限 base change により、R=C[[t]] 上の族を取り、
generic fiber を上記一般 marked pair、special fiber を (H0,F0) としてよい。
必要なら有限拡大後に R を取り直す。

線形因子と四次因子の積の写像は射影空間の積からの proper 写像である。
special D0 の線形因子は H0 だけで、四次因子は既約なので、
generic 因数分解は有限 base change の後に R 全体へ延びる。
因子の係数はそれぞれ primitive に取る。
generic と special の双方で D の essential-variable rank は3である。
勾配の係数行列は R 上で一定階数3なので、その二次元 kernel と三次元 quotient は自由直和となる。
正則な target 座標を選び、族を一つの P²_R 内の H と F の族にしてよい。
これは t の極を持つ座標変換を許した議論ではない。

quartic 族 X⊂P²_R は平坦で、generic と special は共に integral である。
全空間も integral であり、R は excellent なのでその正規化

    ν : Xbar → X

は有限である。以下で、Xbar が実際に同時 P¹ 正規化になることを示す。
この平坦性を仮定として省略しない。

## この族で同時正規化が成立する理由

Xbar は integral normal surface である。
t は零因子でないため R 上平坦であり、ν が有限なので proper である。
normal surface は Cohen–Macaulay なので、その Cartier special fiber は
一次元 Cohen–Macaulay で embedded point を持たない。

X0 の generic point は滑らかである。そこで quartic のある相対偏微分が非零だから、
全空間 X も regular で、ν はその近傍で同型である。
従って Xbar0 には X0 を birational に覆う generically reduced な成分が一つ、重複度1である。
他の一次元成分は有限写像で X0 の点に縮まれないし、X0 全体を覆う第二成分も
generic point 上の同型性から存在しない。
embedded point が無いことと合わせ、Xbar0 は reduced integral curve である。

generic quartic の δ は A3 の2と A2 の1で尽きる。
算術種数3から、その正規化 Xbar_K は genus0、幾何的に P¹ である。
必要な有限 base change の後に P¹_K としてよい。
proper flat family で Euler characteristic は一定なので

    χ(O_Xbar0)=χ(O_P¹)=1.

Xbar0 は X0 と birational で、その正規化は E6 quartic の正規化 P¹ である。
正規化の完全列は χ(O_Xbar0)=1−length(Q) を与えるため、Q=0である。
従って Xbar0 自体が P¹ である。
両 fiber が smooth である平坦族 Xbar→Spec R は smooth で、
三つの相異なる special section を持ち上げて source の P¹_R を trivialize できる。
引戻し ν*O_X(1) は両 fiber で次数4なので O_P¹_R(4) となる。

ここまでの議論では、正規化の平坦性・special fiber の被約性・同時性を
一般の曲線族の性質として仮定していない。この integral rational quartic 族で直接証明した。

## 正則な section の標準形

generic tacnode の二逆像を p_t,q_t、cusp の逆像を ζ_t とする。
有限 base change により三点を section として取る。
これらは全て special E6 cusp の唯一の逆像へ近づく。
実際、特異点の極限は特別曲線の特異点であり、その正規化逆像は一点しかない。
source 座標 r を正則に移して

    p_t=b(t),  q_t=−b(t),  b∈tR,  b≠0 in K

とする。ζ∈tR で、generic fiber では ζ≠±b である。
source の無限遠はこの三 section を避けるように取る。

target の tacnode section を [0:0:1] に固定し、H={h=0} とする正則な target frame を選ぶ。
normalization 上の H section は二逆像それぞれに order2 を持ち、次数4を使い切るから

    h(r)=(r²−b²)²

とできる。係数の scalar は special h の order4から単元である。
y(r) は二点 ±b で消える次数4以下の多項式なので

    y=(r²−b²)(A r²+B r+C).

special E6 の vanishing sequence (4,3,0) より、
h の正則な倍を y から引いた後の r³係数 B は単元である。
この正則な target 変更と B による scaling により

    y(r)=(r²−b²)(r+c),  c∈tR.

となる。w(r) は R[r] の次数4以下の多項式で

    w0(0)≠0

を満たす。

## generic cusp の微分条件が矛盾する

q=r²−b² とおく。cusp の逆像 ζ は ±b と異なり h(ζ)≠0 なので、
正規化 map [h:y:w] の微分が零になる条件は

    (h' y−h y')(ζ)=0,   (h w'−w h')(ζ)=0.

最初の左辺は直接展開して

    h' y−h y'=q²(r²+2cr+b²).

従って ζ²+2cζ+b²=0 であり、q(ζ)=2ζ(ζ+c) である。
generic fiber で b≠0 だから ζ≠0、q(ζ)≠0 だから ζ+c≠0。
二番目の条件へ h=q²、h'=4rq を代入し、
generic fiber 上でのみ非零の q と 2ζ を消去すると

    2 w(ζ)=(ζ+c) w'(ζ).

両辺は R の元なので、この恒等式は special fiber にも延びる。
ζ,c∈tR と w∈R[r] より右辺の極限は零だが、左辺は2w0(0)≠0である。
矛盾した。命題が従う。

## Poisson 成分へ適用する範囲

103番目の一般点はこの marked pair を持ち、全成分はその一般源の PGL 飽和閉包である。
退化因子の写像は、退化因子非零の開集合上で正則である。
Ga9 canonical 点では退化因子は非零であるから、仮にその点が同成分に属すれば、
上で排除した DVR marked-pair family が得られる。
よって Ga9 の密な linear-w 軌道は同成分に属さない。

Ga9 全源の磁気零などの境界が103番目の成分に入る可能性はこの議論で排除しない。
密軌道の非包含は、全源の同成分への帰属を否定するための最弱の結論である。
他の24/25次元親成分への帰属又は新成分性は未解決のままである。
