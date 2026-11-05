# 真正 C[ell] の R0 作用：補正代表による独立監査

主 SOURCE `55507052...` に対する独立報告 `db17520a...` を固定した。
七次元 D 商の r/z 作用をそのまま採用せず、真正 G=kerE23/imE³ へ作用を移した。

低い二つの kernel 座標の ell像は F の三次に入り、
c*r ell²+d*ell³ の既知形式演算子で表せる。
次数二の Y2=c*r ell+d*ell² を引き、さらに残る次数四の ell像を消す
次数三 Y3∈F を選ぶ。未知 Y3 の値を零にせず保持する。
mY3⊂m⁴C4⊂imE³ なので、その作用の genuine class は選択に依存しない。
高い二 kernel class も次数三の補正代表を持つため、m で消える。

最初の七 source について、delta2 と bAi2 の平面三係数、全42個を
旧 C2 と pure-ell 読戻しへ結び付けて確認した。
z_j*rawtarget∈m² と RHi=(0,delta_i,bAi)、bAi1=zi から ziH∈m³となる。
低代表の rH/ellH 部分、r²A/B 部分、選んだ free Y2 の ell²A/B・r ellA/B 部分へ
zi を掛けると m⁴へ入るため、genuine target 上で全 zi は零である。
これは任意の raw 次数二ベクトルへ同じ包含を主張する議論ではない。

r 作用は q(r*完成した補正代表) を計算して求めた。
d*ell² の free 成分は q で消えるが、それ単独を boundary とはしていない。
Y2 の r ell 成分を含めた四 affine operator を独立に組み直すと、
既存 HB 二列の行零・三の負に係数ごとに一致する。
新しい値を pair せず、既に独立検証済みの値から低二→高二の行列

    ( 3/5,     -6 ),
    (121/80, -121/8)

を得た。階数一、r²=0 で、ell は定義から零である。
Jordan 基底と逆変換の両積も確認した。
従って実 G の R0 加群型は R0/(z0,...,z6,ell,r²)⊕k²、
最小生成数三、socle 次元三、m²G=0 となる。

F の R0 安定性、全23作用表、未知 Y3 の係数、全形式 torsion lift、
元比較の connecting 値、元 H0 や nil 三乗を計算したとは言わない。
再現は `code/ambient_B18_C_genuine_R0_transport_source_review_20261105_v63.py`。
