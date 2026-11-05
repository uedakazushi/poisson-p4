# honest M7 の分裂と kernel source：独立純粋検証

独立 SOURCE peer は PASS した。報告 SHA256 は
`7ff33ee4f06498bef805e8d60f4ee7a80f24db74e537cd4cfb3859fddb8a2ffd`。
原値とその独立 VALUE peer を pin し、八個の未計算高座標は全て不定元のまま残した。
新しい source、係数値、jet は求めていない。

実 M7 の基底は h0,h1,h2,h3,h4,rh3,rh4。
受理済みの r 作用と真正 target G4=J2(r)⊕k² を用いて、
全28 equivariance 恒等式を確認した。
S=span(h4,rh4,h0,h1) についても全28の作用恒等式を確認し、
S は M7 内で実 R0 部分加群 J2⊕k² である。

既に値が証明された d=−36/1681、a0=2/25、b1=1331/9216 により、
delta|S の行列式を全24 permutationから再計算した。
残る高座標による項は消え、d² a0 b1=11979/2260608800≠0となった。
従って delta|S は実 R0 同型で、その逆は section を与える。
結論は **M7=S⊕ker delta** である。

kernel の三列 u,rh3,s も別に構成した。既知の非零 QQ 定数 d,b1 だけを割り、
四つの真の未計算係数 p2,b2,p3,b3 を保持する。
全12個の delta*K、全21個の rK 恒等式、および source 行2,3,5の minor1を確認した。
kernel の次元は3なので、この三列は完全な基底で、型は J2(r)⊕k である。

ITold の既知四類 z5h4,z5rh4,r²h4,r³h4は既に実 Kbar に入る。
head3 の低二座標が零なので delta(h3) は真正 socle に入り、rh3も Kbar の実類である。
その M7 像が非零なので、ITold の四類から独立である。
残る u,s の lift と合わせた七つの bar 類を、未計算の真の四係数を使って表示できる。
全形式 K-lift を展開したという主張ではない。

S の安定性は M7 での r²h4=0を使う。一方 Told では r²h4≠0である。
従ってこの証明は Told→G の section、Kbar の ITold による extension の分裂、
Kbar の全 R0 型、元 H0・torsion 指数・nil 三乗を決定しない。

固定の小さい形式行列のみを計算し、0.0105秒、所有上界22,216で完了した。
コードは `code/ambient_B18_C_connecting_kernel3_split_review_20261105_v63.py`。
原資料と report は凍結している。再検証は別作業コピーで行う。
