# 実 Kbar 七基底と typed mod-ell cycle：独立 SOURCE 検証

独立検証は PASS した。報告 SHA256 は
`a5ec34be9b391a8dc9b6a34ce242c962046249301faab7c43deea42cc5cfe71a`。
実HEAD29/30値の独立peer bffdb3…、旧HEAD23値、正直なUCT、
FIVE formal cycle sourceを固定した。新しい応答値は計算していない。

旧Toldの実QQ十一基底は、五headの後に
`rh3,z5h4,rh4,z5rh4,r²h4,r³h4` を並べたもの。
これはexact HF `(5,3,2,1)` の受理済みfiltered basisと一致する。
Kbar=K/ellKの七基底は
\[
z5h4,\quad z5rh4,\quad r^2h4,\quad r^3h4,\quad rh3,\quad U,\quad S
\]
であり、
\[
U=h3-a_Uh1-c_Urh4,\qquad
S=h2-h3/4-a_Sh1-c_Srh4.
\]
全係数は確定済みの有理数で、
`a_U=b3/b1`, `c_U=(p3−a_U p1)/d`,
`a_S=(b2−b3/4)/b1`, `c_S=(p2−p3/4−a_S p1)/d` と定義する。
新しく確定したp2,b2,p3,b3と旧p1,b1,dから、独立に11×7行列を再構成した。

接続写像の未計算値p0,b0,p4,b4は全て形式変数として保持した。
七列のh0/h4座標が零なので、その未知四成分の寄与も恒等的に零となる。
形式変数を零に置かず、full edge×Kの全28成分を独立に確認した。
補集合 `(h4,rh4,h0,h1)` のedge determinantは
`d² a0 b1=11979/2260608800≠0` で、未知四成分に依存しない。

十一座標の方程式を直接解いて明示的逆変換を構成し、
両向き121成分の積が単位行列となることを確認した。
従って七列はQQ一次独立であり、正直なUCTの七次元kernel全体の基底をなす。
この補集合は**QQベクトル空間としてのもの**である。
元のR0加群UCTの分裂や、七個の自由R0基底を意味しない。
元UCTの非分裂証明e4abd…と矛盾せず、その結論を変更しない。

56成分cycleへのsourceも、定数originsだけに読み替えなかった。
FIVEの元node0/1は実mbar/nbar、node2は行22,23,29,30,55を選ぶEであり、
node6は全formal式 `sigma=(I56−nbar*mbar)E` である。
実旧十一基底を作る5×11のRbar monomial mapと上の定数11×7行列を合成すると、
56×5→56×11→56×7のtyped operator root4を得る。
全shape・source順序・元八変数指数を確認した。
受理済みfull768関係のdescentと正直なUCTを使うので、
単なるhomology tensorや欠落したTor像に頼っていない。

このpacketは実mod-ell barclassesと有限typed cycle operatorを与える。
56成分の係数は展開せず、full-R0のK lift、全Kのpresentationや作用行列、
元H0・nil三乗・nil指数は決めていない。
独立SOURCE所有上界21,966、実行0.0304秒、55秒・2 GB・RAW32 MiB内で終了した。
再現コードは `code/ambient_B18_C_Kbar_actual_seven_basis_source_review_20261105_v63.py`。
