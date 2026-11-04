# 最小混合normの一個determinantと旧117次条件とのGCD

2026-11-04。承認されたm1一個の16×16 determinantと、既存primitive a117とのQQ[a] GCDだけを計算する。16-vector solve、他の混合係数、新しい点・素数・GBは使わない。

入力はMIXED_NORM16_INPUTで保存した全16本のmonic除算。H=h1はe-monic16、m1はe-monic19、乗法行列MはQQ[a]上。numeric DでN=D*Mを整数化し、独自scalar-polynomial Bareissで全exact除算を再乗算する。

detからQQ contentと既知単元a冪だけを除く。未知a因子は除去しない。旧a117とGCDを取り、両商とのliteral productを保存する。可能ならcoprime商二つのBézout恒等式を同じ32MiB/60秒内で保存し、巨大なmatrix cofactor bundleは作らない。

actual BbadのG!=0部分はh1=k1=m1=0を満たすので、両normの零点条件を満たす。GCDがQQ単元ならこの部分は空。GCD非定数なら必要a集合が縮小するだけで、実在は結論しない。G=0のQ/L例外はどちらの場合も別に保持する。

各代数閉体点でHはmonic16であり、Hとm1が共通根を持つなら、有限16次代数でのm1乗法は可逆でない。従ってdetN=0。この集合論的必要coverには巨大なfull matrix adjugateや元h1/m1 Bézoutを展開する必要はない。全source除算とdeterminant定義が有限証拠を固定する。

GCDのBézoutはunivariate条件の最大性を検証するためのもので、m2/m3へ自動切替する指示ではない。零determinantならその結果を保存して停止する。

## 今回の確定結果

承認された一回はnative5.8867秒、全体8.6415秒で完了した。全1240回のBareiss exact除算を再乗算した。detNは非零、a次数98、a-order35。QQ contentとa35だけを除くとprimitive63次式になり、旧primitive117次式とのGCDは**primitive13次式**である。未知a因子は除去していない。

両exact quotientを保存し、coprime quotient二つから得る整数Bézoutを全多項式係数で確認した。保存S/Tは104/50項、最大係数109077/108388bitである。巨大係数でも実artifact合計5772222 bytesで32MiB以内に収まった。full16-vector solveやm2/m3 normは実行していない。

従って、今回受理するprimary結論はG!=0の旧badに対する**13次必要a-cover**である。G=0のQ/L例外、actual badの実在・空性は保持する。単なるa根数と明示triple-quotientのmodule lengthとの違いは、同じHを使う二乗法行列のFitting idealを別途確認してから述べる。
