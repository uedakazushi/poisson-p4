# CA generic V境界：八actual odd制約の後で三次は非零

2026-11-03。最初のSymPy版の60秒timeoutは結果に使わない。
保存したnumeric normalizer frameと、独立に構成したsparse150係数Lie Hessianを使う新V2計算は2.60秒・50MBで完了した。
CBのQQ(V)計算は本稿に含めない。

## 成功した有限計算

sourceはordered formal exp(s_i E_i)のG20とfixedCA補fiber9である。
実20 orbit target座標を固定し、29 quadratic-base targetと選択八odd outputをpaired sourceでevenにする。
source kernel coordinateのodd部分はepsilonだけに固定する。
八actual odd制約はraw row27,28,29,35,60,65,67,68で、orbit20と合わせたregular minorは100/6561である。
選択57×57 blockを厳密QQで反転し、全source曲線を別の5×5有限座標pushforwardで再評価する。
全150出力が一致し、mixed Hessianのkernel二乗も独立に一致した。

選択29 base出力はeven、八odd generatorの出力もevenだが、
他の出力に

    raw63: (800/27)*epsilon³,
    raw64: (−40/27)*epsilon³

が残る。
従って八actual regular odd coefficientを除いた一変数nuのidealは、三次の時点でnu単元を持つ。
ここでtracefree t=epsilon、nu=epsilon²である。
raw28 Y=0を固定してordinary cuspと誤認する計算ではなく、八actual odd制約を解いてから読む。
非零三次が出たためp⁵,p⁷の再実行は不要である。

補助計算では全150 actual odd係数のR29 firstjet行列も再構成し、rank9と非零九次minorを保存する。
全normalization Fitting critical cotangentのrankも9である。
これらの計算の独立reviewと、以下の適用鎖の受理までは、generic Vの全完成型を台帳へ昇格しない。

## generic Vへの適用鎖（独立監査用）

1. CA source critical V-axisの全monic表示はP0を含む実代数曲線で、global G作用でもcriticalである。
   そのG飽和閉包の次元はP0での全critical次元20以下である。
   V1のGL orbit微分rank20（scalarを含む）により、その閉包の次元は20、V1 orbitは稠密開となる。
2. P0では全proper normalization fiber長2であり、G不変のfinite・fiber長≤2開近傍がある。
   その開近傍は上の稠密orbitに会うため、V1 orbit全体でもfiniteである。
   V1のsource微分corank1により選択source点のfiber長は2以上なので、全fiberは長2一台に確定する。
3. source normalizer Bはregular29次元で、20 orbit方向はtargetでも独立である。
   全target critical schemeはcotangent rank9、一方orbit20を含むのでsmooth20となる。
   regular actual base R29上のJ=(I,nu)は、商がregular20次元なので九生成regular idealである。
4. actual Iのfirstjet rank9はJの全九classを張る。NakayamaでI=Jとなり、R座標を
   I=(s1,...,s8,nu)と取れる。
   actual completionは

       C[[g20,s8,t²,s1t,...,s8t,t³]] ⊂ C[[g20,s8,t]].

   20 orbit方向を使うformal action/sliceはtarget成分のG不変性から全次数で同値であり、
   未検査のimage/base-change交換を行わない。
   非零scalarを除いたprojective版ではg変数は19個である。

この適用が受理されれば、CAのgeneric V境界はordinaryな指数1のpinched型となり、
以前のV枝長3・cusp(2,7)という候補は採用できない。
ただしB18原点でのfullCA idealや、全三変数cubic f3(U,V)をこの一点から決定してはいけない。
U-axisでのU³ witnessは混合UV²項を検出しない。

* 成功V2: `code/ambient_CA_V1_odd_order_v2_20261103_v63.py`
* evidence: `reports/AMBIENT_CA_V1_ODD_ORDER_V2_20261103_v63/`
* actual odd rank9 supplement: `code/ambient_CA_V1_actual_odd_rank9_20261103_v63.py`
* evidence: `reports/AMBIENT_CA_V1_ACTUAL_ODD_RANK9_20261103_v63/`

初回failure、cheap checkpoint、newHessian helper、V2成功の全sourceとハッシュを別々に保存する。
