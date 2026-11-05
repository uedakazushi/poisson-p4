# 元 connecting の plane 次数と残存 sigma²：独立検証

SOURCE 検証は PASS した。独立報告 SHA256 は
`d91a97f75207a0c9d1369f31d4d84d9892ab5a97893ef1914799c71215815661`。
五 head の値や、新しい physical/gauge/jet 係数は計算していない。

原 Q48 の各 physical 関係は
Fj=delta+bA tauA−bB tauB で、tauA²=nuA、tauB²=nuB である。
基底1,tauA,tauB,tauA tauBとの積の三非定数成分を独立に使い、
既存12行の C2 と pure ell² readback から全56列の低次数を復元した。
全 origin、全 pure-plane first は零で、A 成分の平面二次も全列零である。
ここでは z 方向の全 first まで零とは主張しない。
第一七 U は実 bA source alias、末尾 U は gA の厳密な零七座標 tail であり、
両者の平面二次までが零であることを別々に確かめた。

B と H の平面二次行列は、三 pivot physical18,60,67に対応する保存行列と一致した。
q7 の真正四自由座標が読む19個の分子 coefficient functionalを取り出し、
sigma が ell 非依存である条件の下、componentwise divisor を全列挙した。
元の product routing と全て一致する。

必要な sigma jet は sigma0、sigma1_r、sigma1_z0、sigma1_z4、
sigma1_z5、sigma2_rr のみである。sigma³ が必要になる経路は正の ell 次数を持つ
V の一次平面係数に掛かり、これは実 readback により零である。
一方、B の高座標 N(r²ell²)−(96/275)N(r³ell) では

  (V2_B,ell²−(96/275)V2_B,r ell) sigma2_rr

が残る。この adjoint は元 source 列0,16,20の三つで非零であり、省略できない。
「sigma³ が消える」ことを「sigma² も消える」ことへ拡張していない。

最高未解決部分は A/B の V4 sigma0 の r³ell,r²ell²、
coherent な V3 sigma1、混合 V2 sigma1、およびこの B2 sigma2_rr である。
旧 r/g は実正則関数参照のまま保持する。source4 全体や sigma² の cache が
既に存在するという仮定も置かない。

U9 と OLD640 multiplier の namespace の限定訂正は同時に検査した。
U9 の末尾は gA505938、OLD640 の末尾は f63=7790であり、全七 telescope を保持する。
FIVE_FORMAL の過大な説明だけを訂正し、凍結された実参照や数式は変更しない。

再現コード、実行資源、境界降下・形式除算の証明は同時に保存した
FIVE_FORMAL_SOURCE_REVIEW ノートを参照されたい。
この結果は source/type/order の証明であり、接続値・元 H0・nil 三乗の新結論ではない。
