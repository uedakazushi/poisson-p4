# 正weight6軸の実8変数模型・全initial4の独立検証

2026-11-03。対象は正軸eta≠0の物理的canonical点と、その全Poisson完備局所環である。
二つの凍結producer Q_POSITIVE_WEIGHT6_QUADRATIC_MODEL と Q_POSITIVE_WEIGHT6_INITIAL4 を検証する。
接空間・二次錐・真のinitialの次数4まで・実局所次元の上界を区別する。
有限8変数多項式idealと全形式idealの等号、全初期idealの生成、実枝・実次元等号は結論しない。

## 1. 独立に構成した物理的モデル

reviewerは五物理変数の二次bivector150係数から、25個のdivergence線形係数と
全350個のstrict Jacobi係数の定数二次積表を独立に作った。
producerのlinearization、bilinear、initialの行列生成関数は呼ばない。
divergence行列rank24と保存C150×126の全kernel包含・独立free行から、canonical126方向の完全性を確認する。

GL25のpushforward微分とscalar一方向も物理的係数から再構成し、保存150×26行列と一致する。
GL rank22、scalar込みrank23、23orbit cut minorの非零を確認する。
canonical空間の残る103方向との微分の直和分解から、形式逆写像定理による
23滑らかaffine因子、projectiveでは22滑らか因子を除く操作が正当となる。
非零pointのscalar方向はorbit+scalar空間に含まれ、射影化でちょうど一本を除く。

物理的定数二次積表から元Jacobi350×150の微分を全QQで構成し、103切断上の保存Jacobianと照合する。
95×95逆行列の厳密積と、8独立kernel方向を確認するのでrank95は上下から確定する。
従って全canonical接空間31、projective接空間30、横断埋め込み次元8である。

選択95本は元の式そのもので、線形blockが可逆である。
形式陰関数定理で95dependent座標を八自由座標の全形式級数として一意に解ける。
残る255元式の全残余を保つので、全次数の**形式presentation**として正しい。
有限jetを全次数の多項式模型と見なす議論ではない。

## 2. 全350係数・次数0から4の独立native QQ replay

独立に生成した物理的二次積表をC++/FLINTへ渡す。
入力はP0,R1,R2,R3,R4の八変数斉次係数で、degreeを分けた正確な畳み込みを行う。
各元Jacobi350行について、全八変数QQ係数を次数0,1,2,3,4でそれぞれ比較する。
全1750個の物理的多項式と、divergence25×5次数を検証する。
方向標本や有限体の階数で代用しない。次数4を超える項だけを省く。
raw方向の係数から計算する別engineなので、producerの二次補正・三次/四次bilinearのsignとindexも検査する。

二次補正は−B^(-1)Q_selected、全残余はQ+J_dependent*H2である。
選択95行が全次数2–4で消えることも確認する。
これにより保存二次5式、全三次・四次残余が実8変数形式idealの正確なjetであることを検証する。

## 3. 真のinitial3/4の完全性

五つの独立な二次leading行Aを選び、残る345行からその定数組合せを引いてorder≥3のGを作る。
この可逆定数行操作は形式idealを保存する。
degree3の全sourceはS1*A2とG3であり、reviewerは全385×120行列を独立に作り、rank48を確認する。

保存kernel337×385について零積と独立rank337を確認し、385−48=337から**全kernel**と分かる。
同じkernelを(S1*A3,G4)へ掛けた全carryを再計算する。
degree4空間はS2*A2+S1*G3+carryで尽くされる。
下次数の倍のrank196、全rank204、新しい四次8類を確認する。
三次新類も8本で、実Hilbert functionは次数4まで正確に(1,8,31,72,126)である。

完全性の理由は次の通り。degree2でA2は独立なので、order4以上のideal元のAへの定数係数は零。
degree3が消える条件はAへの一次係数とGへの定数係数が上記全kernelに入ることに等しい。
その次次数のcarryと、Aへの二次係数・Gへの一次係数以外にdegree4 sourceはない。
各元式の最初の項だけを集めた計算ではない。

## 4. 小さいGBと実上界の範囲

保存八変数QQのMacaulay2入力を新しい出力先で再実行する。
二次idealのradicalは(z1,z7,z5(z4+2z6))、二つの5平面の和集合である。
元五式の包含と、z1^3,z7^3,(z5(z4+2z6))^2の逆包含証明からradical等号を確認する。
二次ideal自体の次数3は、被約二平面の次数2と同一視しない。

二次5本・新三次8本・新四次8本の有限斉次idealはdim3、degree6。
全実initial idealはこれを含むため、実横断次元≤3、projective滑らか22因子を戻して
**全projective Poisson局所次元≤25**を得る。
dim3やdegree6を実局所環の等号とすること、全枝の存在・個数・尽くしは主張しない。
非零eta全体は保存torus作用で一つのprojective軌道なので、上界は全eta≠0へ移る。
eta=0のQa点は別の局所問題である。

review code: `code/q_positive_weight6_small_model_review_20261103_v63.py` と同名`.cpp`。
receipt: `reports/Q_POSITIVE_WEIGHT6_SMALL_MODEL_20261103_v63_REVIEW/`。
