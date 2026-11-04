# Actual E4の二次模型の独立監査

2026-11-04。対象は凍結された `AMBIENT_B18_ZERO_NORMALIZER_QUADRATIC_20261104_v63` である。既存FULL9の57未知数×55係数を読み、新しいsource求解・高次jet・旧r/g合成・点・素数・GBは使わない。actual source、unit系、normalizer tensorの既受理reviewをhashで結ぶ。

二族のsource DAGに必要な祖先はCA 9851、CB 10092ノードである。producerの配列jet演算を使用せず、次数≤2の単項式をkeyとする独立有理係数dictionaryで和・積・有限幾何級数によるunit逆元を実装する。それぞれ元の48式とactual F9九式を全55係数で照合する。F9は族ごとの八つのzero odd生成元から七つと、二つのpositive odd生成元を選ぶ同じ実source函数である。両trace基底は別のまま保持する。

九physical R0 readbackの一次行列と保存逆行列の両積を確認し、det CA=3/128、CB=96を独立有理消去で確認する。保存された二次までの逆写像について、x(Y(x))=xとY(x(Y))=Yの両方を全55係数で検査する。unit inverseの一意性により、これが同じphysical R0への正しい座標変換である。nuと十二のphysical出力のeven/odd pairを、この共通座標へ別の多項式合成で戻す。元sourceから読み直した全係数が証明書と一致することを確認する。

十一変数(x27,x108,x109,x112,x113,x117,x118,x119,x61,tauA,tauB)で、二つのmonic式と十二output差を独立に再構成する。odd係数の定数項が零なので、二次出力差にはその一次部分だけが寄与する。二族のtauを同一視せず、右辺はそれぞれnuAとnuBを使用する。

十四式のlinear rankは2である。monic二式と(x27,x61)のminorは−20/9。まず一次置換を行い、選択二式の二次残差にこの2×2逆行列を掛けた二次補正を消去変数へ加える。この補正を含む全十一置換、および全十四式の補正後jetを照合する。選択二式は二次まで零になり、他の式も全てorder≥2となる。

残る九変数は

`x108,x109,x112,x113,x117,x118,x119,tauA,tauB`。

全二次係数の独立QQ行階数は5で、その線形包は正確に

`L*(x108,x109,x113,x118,x119)`、

`L=5*x112−20*x117+24*tauA+6*tauB`

である。全十四式について保存されたlinear quotientとの積を確認し、五つの標準積との両側spanをrankで確認する。

actual E4の全十四式という入力は、別に受理されたnormalizer tensor模型に依存する。unit linear eliminationを全次数で行った後は、生成元は全てorder≥2である。そのidealの二次initialは、各生成元の二次部分の定数係数線形包と一致する。非定数cofactorや相殺は新しい二次方向を作れないため、このrank5は二次initial全体を表す。embedding dimensionは9となる。

二次initialの共通因子Lから、全次数idealがLで割れるとは結論しない。actual N pair liftによる全次数hypersurfaceの同定、四変数への残り函数、E4のsocle/depth、D0からの単射性、D0のdepthはこの監査の範囲外である。生成元の二次情報が正確でも、閉点embedded部分を排除する根拠にはならない。

専用reviewコードはsource cacheを再評価するだけで、producerやsource solverをimport・実行しない。全input/evidence hash、実行時間、メモリは同名REVIEW/report.jsonに記録する。
