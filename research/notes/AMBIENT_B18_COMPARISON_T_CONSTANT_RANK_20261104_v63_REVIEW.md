# 実Tのμ=5：独立constant消去・literal wedge・rank監査

2026-11-04。全source X0/X1と116×936提示の独立peerを固定し、作者の消去を再実行せず、
保存済みの固定rank証明書を別演算で検査する。新source、jet、rank探索、348行列は作らない。

全768列の定数を確認し、rank27の保存constant source変換で得る741列X′が零定数を持つことを
全て確認する。従って `(m X′)_1=m0 X′1` であり、mの未知一次は不要。
27個のunit source/target消去で、残る列を補正する係数の定数は零なので、selected27列の
未知top wedge定数は残るwedge定数を変更しない。

独立checkerは、literal d2の全84列から一次係数を取り、84本の微分座標が正確なidentity minorに
なることを確認する。各cycleのw0をこのanchorで抽出し、作者のskew差/2公式は使わない。
741列×24×8=142272係数を全て再合成し、対角・skew・r微分を含め一致を検証する。
Wは84×741、非零715係数である。

保存された正規化basisのdistinct leading row、全741列のfactorization、全固定引算・正規化を
検査する。basisをleading row順に並べるとunit triangularで、選択元列のfactor行列も
記録された非零divideを対角に持つtriangular行列となる。従って固定84minorのdet積が厳密に非零。
新しいGaussianやminor探索なしにrank84の両不等式が従う。

complete提示から27 unit対を消すと89 target、残るconstant提示はWと零5行なので
`mu(T)=89−84=5`。これは長さ5や深さの結論ではない。

## 次の小さいSchur判定の数学的条件

84 wedgeのunit blockをSとして消す。残る5行をY、upper blockをWとするとY0=0であり、
残り657列の一次は

\[
R_1=Y_{rest,1}-Y_{S,1}W_{S,0}^{-1}W_{rest,0}.
\]

Y0が零なのでW1はこの式へ寄与しない。先の27 unit消去の非定数係数も、remaining列の
X′0が零であるためY1とW0には影響しない。元d3はbottomが恒等零で、topは一次から始まる。
従って、その実残余関係の二次初期項は厳密に
`−Y_S1 W_S0^-1 d3_linear` であり、source X1とliteral d3一次だけで得られる。

もし残りlinear rank≤36なら、`dim mT/m²T≥4` なので長さ≥9となる。
もしrank37なら、この層は3次元。実linear関係の8変数倍と168本の上記実d3二次が
全5×Sym²(8)=180係数を張れば、`m²T=m³T`。m²Tは有限生成なのでNakayamaにより零となり、
長さは5+3=8である。この**十分条件**に未知二次関係を落とす問題はない。
反対に二次spanが足りなくても、未入力の実関係が残余を消し得るため、非零の高次層を結論しない。
本peerはこの理論を認めるが、次段の実linear/二次rankは計算しない。
