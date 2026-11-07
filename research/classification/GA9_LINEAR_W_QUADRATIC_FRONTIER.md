# Ga9 c=0 linear-w全源：既済二次理想の支持と停止点

2026-11-07、PRIMARY E。Ga9 c=0のpower1/power3は新24成分への全実族で既に被覆した。
本稿は残るlinear-w全源だけを扱う。新23成分の最弱判定を一度確認し、
**二次理想がm-primaryでないため、この経路を終了する。三次には進まない。**

## 既存の全周囲証明書を再利用する

入力は不変保存済み
`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL211_ZERO_C_LINEAR_W_QUADRATIC_20261107.json`。
その既済全150変数・全350式の結果は

    rank D=96,  dim ker D=54,
    rank(actual PGL/Euler/scalar)=rank(whole source)=48,
    normal数=6, 射影軌道次元23, 射影接次元29。

指定点は元の全十一parameter源のcanonical値を持ち、六磁気項も保持する。
旧kernelと保存normal六本を用い、実作用48と合わせたrank54・非零54小行列だけを確認した。
これによりnormalが真正補空間であることを確かめた。
**D、二次Jacobi列、旧GBを新たに計算していない。**
全parameter源がこの密な23軌道の閉包に入るという既済の全源定理もそのまま保持する。

## 保存された二次七式と十六式は同じ理想を生成する

normal座標は旧証明書の順序 `e0,…,e5` をそのまま使う。
元の全二次余核は七式で、旧計算はlex basis十六式を保存している。
今回、保存十六式の各々を元七式の多項式結合として、
逆に元七式を保存十六式の結合として、次数二〜四の有限Macaulay線形算術で検査した。
全乗数を公開reportへ明示し、各多項式恒等式も直接展開して零を確認した。
これは既存のbasisを新たに探索するGB計算ではない。

## 被約二次支持は二平面と一本の直線

保存式はまず `e4^4` と `(e2+2e4)^2` を含むため、複素被約支持上でe4=e2=0。

e5=0なら残式のe3²からe3=0となり、他には制約が無い。従って全二平面

    P={e2=e3=e4=e5=0}

が残る。
e5≠0なら、保存式の

    e5²(e3+27e5/20),
    e5³(e1+27361807e5/380),
    e5(e0+26019e5)  （e2=e4=0,e3=−27e5/20の後）

から、残る全支持は

    L=s*(-26019, −27361807/380, 0, −27/20, 0, 1)

という直線に入る。PとLの双方へ元七式・保存十六式を代入し、全て恒等零と確認した。
e5=0とe5≠0以外の補集合は無いので、全複素被約支持は正確にP∪Lである。

したがって二次理想のheightは4、被約次元は2で、原点のmaximal idealにprimaryではない。
特に `(1,0,0,0,0,0)` が非原点の零点であり、新23成分を認定する二次判定は失敗する。
これは点や素数を増やして検査した結果ではなく、保存された全理想の支持計算である。

## 得られる上界と、得られない結論

実48方向を切る真正102次元線形切片でDの独立96式を形式陰関数で消去すると、
normalは六変数になる。その全残式の初期理想は、今回の二次七式を含む。
よって真のnormal支持の次元は二以下、Krull不等式から全局所raw次元は48+2以下、
射影局所次元は25以下となる。

高次初期式はP又はLをさらに縮め得るため、PやLを実Poisson変形とは呼ばない。
二次支持と真の全tangent coneの同一性も主張しない。
この上界だけから新23・24・25成分のいずれも数えず、既存成分への帰属も推論しない。
直線に対応する全150係数の一次方向はreportに保存したが、Jacobi曲線の探索はしていない。
三次・高次の延長、normal方向の点列probe、固定N軌道極限ルートの再開を行っていない。

## 停止済みD5曲線の局所障害を保持する

別途検討した係数族は

    F_eta=h²(h+eta*y)w+C*y4+beta*h*y3+gamma*h²*y²,
    C≠0、eta≠0

で、HはD5 quarticのcusp tangentとなる。
chart w=1でhFの主weightは、weight(h,y)=(3,2)に対して11。
主部は `eta*h³y+C*h*y4`、次weight12の項は `h4+beta*h²y³` である。
主部のJacobian quotientでは

    h4=−(4C/eta)h²y³,

かつh²y³は非零のweight12 classなので、その次項classは
`(beta−4C/eta)h²y³` となる。
標準的なsemiquasihomogeneous normal-formの適用の下でτ=μ=12に必要なのは
`eta*beta=4C` である。
従ってC(0)≠0、betaが零点で正則というこの係数族は、eta=0のGa9 E6 canonical形へ
free-quintic smoothingとして到達しない。
この局所normal-form適用は今回の二次支持certificateの入力ではなく、独立再現も未実施である。
全D5模型・全非線形族の不存在とは拡張せず、同曲線の係数・次数変更を繰り返さない。

## 保存物

共有runnerのisolated mode `--ga9-linear-w-quadratic-support`、
`jobs/ga9_linear_w_quadratic_support.json`、
`CLASSIFICATION_GA9_LINEAR_W_QUADRATIC_SUPPORT_20261107.json`。
追加前の38関数はmainを除いてAST不変、旧source/report bytesは全て不変である。
全源のE帰属は未解決のまま、次は全軌道に効くglobal inverse認識又は別のPRIMARY Eへ進む。
