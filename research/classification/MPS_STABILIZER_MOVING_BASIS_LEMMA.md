# Normalcone と unipotent 安定化場を混ぜる実 PGL 族

2026年11月7日。未被覆の double_A/Ga9零c源のための道具。
この補題は全族の正則性を証明する。現在の候補planeから全源帰属が得られたという主張ではない。

## 1. 全poleを先験的に消す構成

五座標を base 三個と fiber 二個に分け、base重み0、fiber重み1とする。
既知閉PGL安定集合 K に属する実normalcone族を

    Qε=Q0+εQ1+…+ε³Q3

とする。各Qjは二次斉次二ベクトルで、fiber-scaling重み j−1 を持つと仮定する。
これは既知二次tensorを、二つのnormal座標ε倍とtensor ε倍で引き戻す通常normalconeに成立する。
Lagrangian plane条件により負冪が消えていることは別に確認する。
固定のEuler鉛直項を引いてQ0だけをh-Casimir表示へ直すことも許す。

線形場Nがbaseを固定し、二fiberにbaseの一次式を加えるとする。
従って行列として N²=0、fiber重み−1である。
さらに **[N,Q0]=0** を仮定する。
ε=δ²と置き、実PGL変換 `gδ=I+cN/δ` を作用させる。
Lie微分をL_Nとすれば

    gδ_*Q_(δ²)=Σ_(j,k) c^k δ^(2j−k) L_N^k Qj/k!.          (1)

二次斉次bivectorのfiber重みは最小−2である。
従ってj≥1に対し `L_N^k Qj=0` for k>j+1。
残るδ指数は `2j−k≥j−1≥0` であり、j=0の負冪は安定化仮定で零である。
よって(1)は**全δで正則な実多項式族**で、端点は正確に

    Q0+(c²/2)L_N²Q1.                                      (2)

追加項のfiber重みは−2なので、純粋な二fiber方向のmagnetic二次式である。
全δ≠0で元の既知族とPGL/scalar同値だから、閉性により(2)もKに属する。
これは固定tensorだけのunipotent軌道極限でも、Poissonのformal jet延長でもない。

Q0が既済の全六磁気源のzero-magnetic点で、translation cokernelが一次元の場合、
(2)の追加磁気式をその**完全cokernel**へ射影する。
非零なら既済normalizerによるgeneric磁気軌道へ入り、その軌道が全源に稠密であることが確認済みなら、
全源の閉包もKへ入る。これらの適用仮定なしにmagnetic項の非零だけで帰属を主張しない。

## 2. Lagrangian plane schemeとの関係

固定planeの近傍を六変数graph Mで記述し、二つのnormalの括弧の制限を
六つの二次式係数 C(M) とする。magnetic fiber translationの一次作用は、
このLagrangian条件の微分 dC と一致する（固定可逆frameを通して比較する）。
NがQ0を保存すると、対応するgraph方向 B_N は ker(dC) に入る。
式(2)のmagnetic cokernel classは、この方向の二次Lagrangian障害に対応する。
したがって、Lagrangian schemeがその点で滑らかなら、B_Nは実Lagrangian曲線に延長でき、
二次障害はim(dC)に入るので、essentialな磁気classは零になる。

非零classを探すには、単に別planeの数値を増やすのではなく、
対応する非被約又は非滑らかなLagrangian局所構造が必要である。
translation rank5・cokernel1の double_A では、fullGrassmann tangent次元は1となる。
もしLagrangian schemeの局所次元も1なら、その点は滑らかで上のclassは必ず零。
局所次元0でtangent1の非被約孤立planeのような場合は、この理由では排除されない。
この必要条件は一般のmoving-basis arcや異なるweight flagを排除しない。

## 3. 今回確認したT37 planeは対象が異なる

旧cell013の一つのrational plane族は、t=1で制限因子が `H²C3`、`H=A−2B` となるが、
C3のcusp `[H:B:C]=[0:−6:1]` がH上にある。
double_Aのcuspは唯一のdoubleHの外にあるので、この源とは異なる。計算jobは追加していない。

未処理cell012の第二代数curveの断面 a4=0 は全五共役点

    a5=t, 31104t⁵=1,
    (a0,a1,a2,a3,a4,a5)=(18t³,−1/(18t),−18t²,6t²,0,t)

である。sourceの座標重み `(0,3,2,4,1)` による
`oldXi=t^wi vi` を使うと、T37 tensorの全単項式のt指数が5の倍数となり、
全五点は共通のQQモデルへ厳密に降下する。
そのrational planeは

    v=(A,18A−B/18−18C,6A+C,B,C).

実normalconeを全係数で構成すると制限因子は、非零定数を除いて

    (5184A−B+864C)
    * (1679616A²+2592AB+1119744AC−B²−432BC−186624C²)².

これは H*Q² であり、double_A の2H+C3、Ga9零cのH+既約F4のいずれでもない。
従って式(2)の対象source認識・磁気class計算へは進めず、この有限locusを停止した。
missingcurve全体の因子型、全Lagrangian plane又は全MPS成分の非所属を結論していない。

公開のexact producer/report/jobは `classification_double_a_t37_algebraic_plane`、
`CLASSIFICATION_DOUBLE_A_T37_ALGEBRAIC_PLANE_20261107.json`、
`jobs/double_a_t37_algebraic_plane.json`（0.235秒）。
全入力はpublicのT37完成族とexplicit plane・五次Galois関係であり、/tmpの旧出力は必要ない。
全double_A源とGa9零c三源のE帰属は依然未解決である。
