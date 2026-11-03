# W1 の actual CA/CB/N 接空間と交叉の独立検証

2026-11-03。W=1の同じ物理点で、source二次混合項と各sourceのkernel曲線一つの三次係数を
独立に再構築した。結果は `dim T_CA=T_CB=38, dim T_N=36, dim T_D=37`。
被約交叉はこの点でNであるため、nilradicalのcotangent像は正確に一次元である。
全非零Wへは既受理のschemeを保つGL/scalar移送を使う。
全Dの環型、depth、nilpotency index、nilの台の次元、新しい埋込成分は主張しない。

## 元のsourceと座標の独立再構築

producerを呼ばず、標準有理数と150係数の添字だけを用いる。
CAの固定20空間、NのC2=0固定19空間、CBのf3係数・H16係数・ratioをraw式から作る。
Euler補正とh/x交換も独立の係数操作で実行する。
CBの `g=−ratio*f` は最初から掛け、基点ratio=1でf方向の微分が変わること、
ratioとfの二次混合項が残ることを含める。
基点の全150係数は

    P_W=P0+x²∂z∧∂w

から直接構築し、三つの固定sourceで同一になることを確認する。
proper旗の9方向は既受理のactual source証明書に固定されたものを使い、
そのGL一次列を独立した添字Lie演算で作る。
全150×29（Nは150×28）source Jacobianの全要素が一致し、rankは28/28/27である。
保存kernelに掛け戻して零であることと、kernelが正確に一次元であることを確認する。

## 全mixed二次とkernel三次

source方向vを固定方向速度delta、GL行列M、CBの二次項Qに分ける。
CBでは `Q=−v_ratio*Σv_fi*G_i` であり、この項を省略しない。
その真のsource曲線は

    exp(sM) · (P_W+s*delta+s²*Q)

である。有限5×5行列の指数と逆指数、二次空間座標の代入、exterior作用を別に実装する。
各source基底e_jについて二次係数のpolarization

    q2(k+e_j)−q2(k)−q2(e_j)

を計算し、全kernel mixed列の150係数を保存配列と照合する。
さらに各sourceのkernel曲線一つだけを三次まで計算し、その三次係数を

    L_M³(P_W)/6 + L_M²(delta)/2 + L_M(Q)

という独立添字Lie式とも比較する。四次sourceは計算しない。

これらが真のtarget接方向に入る理由も必要である。
target定義式Fをsourceへ戻すとmixed係数では `df(k)=0` によりFの二次部分の寄与が零となる。
またkernel曲線の変位はO(s²)なので、Fの非線形部分はO(s⁴)であり、
その二次・三次係数はいずれもdFを消す。
一般sourceの任意の三次係数を接方向とする議論ではない。

## 全接空間であることと交叉

二次までの全列のrankは37/37/35、各kernel三次を一つ加えたrankは38/38/36となる。
選択した独立列と単位minorの全有理数値を別に検査する。
既受理の全非零Wの全完成環はCA/CBがAc8、NがAc7であり、
cone埋込次元の**上界そのもの**が38/38/36である。
このactual finite source・全kernelの結果により、作った下界が全target接空間に達する。
normalizerの微分だけを全接空間と誤認した計算ではない。

独立した有理数rank検査で

    rank(T_CA+T_CB)=39,
    T_N⊂T_CA, T_N⊂T_CB

を確認する。したがってscheme交叉の接空間は両接空間の共通部分で、次元は
`38+38−39=37` である。保存された37方向も両空間に入り、独立であることを検査する。

基点をh=y=z=w=0,x=1で評価したbivectorの行0,2,3,4のminorは1であり、rank4である。
generic rank2成分Rには入らない。独立受理された全被約交叉 `N∪R` の定理により、
この点の被約germはN、従ってその接空間次元は36である。
actualDからD_redへのcotangent全射のkernelは
`(Nil(D)+m²)/m²` なので、その次元は正確に37−36=1となる。

既受理の `r^-1 diag(1,1,r,1,1)_*P_W=P_(W/r)` はCA,CB,Nを個々に保ち、
和idealで定義するscheme交叉Dも保つ。よって同じ主張を全W≠0に移せる。
W=0への特殊化や全D型への拡張は行わない。

専用コード：`code/ambient_W1_actual_tangent_intersection_review_20261103_v63.py`。
flat receipt：`reports/AMBIENT_W1_ACTUAL_TANGENT_INTERSECTION_20261103_v63_REVIEW.json`。
上限60秒・2.2GB。producer再実行、GB、新しい高次source計算は行わない。
