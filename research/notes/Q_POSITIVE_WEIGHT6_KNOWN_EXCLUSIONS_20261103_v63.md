# a=8/9の正weight6軸は非零点でactual T01, CA, CBの外にある

2026年11月3日。全global37模型の閉weight6 sliceで見つけた二本の軸のうち、
正の軸を対象とする。負の軸a=−8/9は既存のT40九変数graphと全150係数で同定済みである。
ここでは正の軸の非零点が三つの既知被約成分に入らないことだけを証明する。
他の既知成分への所属、新成分、全Poisson局所環の枝・非平坦性は結論しない。

## 1. 対象は実際の全次数Poisson曲線

global37入力は `a=8/9,z33=eta`、他の全正weight入力を0に置く。
旧物理座標x0,...,x4の150係数で、非零なのは

    64=−256eta/3, 83=−64eta, 90=8/9, 105=−2/9,
    116=64eta, 121=4/3, 137=2/3, 140=1

である。従って上三角の非零tensorは

    P12=−(256/3)eta*x0*x4,
    P13=−64eta*x1*x4,
    P14=(8/9)x0²,
    P23=−(2/9)x0²+64eta*x2*x4,
    P24=(4/3)x0*x1,
    P34=(2/3)x0*x2+x1².

QQ[eta]上でdivergenceと全350 projective Jacobi係数が零であり、
actual quinticは恒等的に `−16*x0^5/81` である。
eta=0は既存のQa=Q_(8/9)である。

固定の非零a正規化frameによってS61の全61出力は、one-basedで

    u4=−2048eta/27, u16=256eta, u25=−3/2,
    u36=27/8, u38=81/16, u57=−256eta

となる。原134式・全150再構成・projective Jacobiも確認する。

## 2. actual T01全体で消える一つの多項式

既存の `TORIC_ZERO_FIBER_COFRAME_QUADRICS_20261102` は、
primitive一次部分16座標とclosed projectionを施した二次部分45座標に、
actual T01全体で消える274本のweight4関係式を与える。
これらはdiscriminantで局所化した認識式ではなく、dense coframe上の多項式恒等式と
Sp4不変性から被約像閉包全体で消える多項式である。

本packetでは全274本を保存済み15本のhighest vectorから再構成し、
四つのroot作用の814個の非零像が同じ空間に入る全QQ係数検査と、
六変数dense coframe上の全係数零性を新たに再実行する。
coframeがSp4右作用について同変で、その六変数chartとSp4がGL4の稠密開を覆うので、
選んだ式が全actual T01で消えることを再確認している。

本曲線の一次部分は既にprimitiveで、perspective補正alphaは全4項が零。
closed jetのzero-based非零座標は

    l5=−2048eta/27, l11=−3/2, l15=27/8,
    q26=192eta, q33=−64eta, q41=81/16,
    q45=−64eta, q56=−192eta.

ここで添字は合計61個のclosed-jet座標内の番号である。
producerの分母払いに使う4Piと、実際のclosed座標Piを混同しない。
単位coframeで一次を10倍・二次を100倍とする規約を有限tensorから別に較正する。

highest weight(1,1)の第0式Rは403項を持つ。
この曲線に代入した際に残るのは正確に次の四項だけである。

    R|curve = (301/1920)l5²*l11*l15
             −(17/240)l5²*q41
             −(97/480)l5*l11*q26
             −(1/3)l5*l11*q45
            = −(1163264/135)eta².

従ってeta≠0ならRは非零となり、全非零軸はactual T01の外にある。
全274式では正確に39本が非零となるが、排除には上の一本で足りる。
eta=0で全て零になることは既知のQa所属と両立する。

正則冪零Nのcommuting-powers表示を直接探す方法では、row0=0,column3=0を課した
linear stabilizerが二次元で、そのNは全てN²=0となった。
その方法の失敗自体はT01排除の根拠ではない。上のglobal多項式が排除の証明である。

## 3. 全proper flag-sourceを使うCA/CB排除

S61内のCA/CBの**全被約成分**は、独立監査済みのclosed proper Sp4 flag-sourceで覆われる。
preimageのline flagはintrinsic K_betaのkernel内に入る。
本曲線では、標準4座標順(x,y,z,w)、二形式の23成分を除く五行で

    K_beta = [ 0, 0, 0, 0;
               0, 0, 27/8, 0;
               3/2, 0, 0, 0;
               0, 0, 0, 0;
               0, −2048eta/27, 0, 0 ].

保存する三次minorは非零有理定数倍のetaである。
従ってeta≠0でrank3となり、**全projective flagの候補はe_wだけ**である。
e_wをreference e_yへ送る固定symplectic permutationを全61出力に有限作用させる。
その後のCAとCBの全closed beta-image constraintsには、両方とも
`2048eta/27` と非零定数 `3/2` が現れる。
ゆえに唯一の候補も両sourceの像に入らない。

したがってeta≠0で全actual CA/CBのいずれにも属さず、その共通被約Nにも属さない。
これは一つの旗chartだけでの失敗やsource微分rankからの推測ではない。
全proper sourceの必要なkernel条件により、未検査の旗候補は残らない。
またCA/CBからの排除をconstant-translationの未証明なheuristicに依存させていない。

## 4. 保存とscope

* `code/q_positive_weight6_known_exclusions_20261103_v63.py`
* `reports/Q_POSITIVE_WEIGHT6_KNOWN_EXCLUSIONS_20261103_v63/report.json`
* 同ディレクトリの `certificate.json.gz`

証明書は物理150係数、全350零、S61曲線、closed jet、global Rの403項、
274式の全評価、K_beta、唯一旗への有限変換、両base constraintsを保存する。
負の軸との同定や閉sliceのprimary decompositionは既存packetを変更せず保持する。

正の軸がT02/T03など他の既知成分に入るかどうかは未解決であり、
本結果から既知成分数やglobal lower boundを増やさない。
