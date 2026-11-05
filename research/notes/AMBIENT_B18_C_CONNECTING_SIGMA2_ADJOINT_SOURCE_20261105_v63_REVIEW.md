# 残存 sigma² の三行 adjoint：独立 SOURCE 検証

独立検証は PASS した。報告 SHA256 は
`4870a2da98e92ab951af0e63fd4a10f58701f113e9af96fe9ffac91a468bb04a`。
原 SOURCE 報告は `81abe7ed97baf663cc98c5c9d50d795c27f3900667cff6dab6f92d5c2511b852`。
既存定数からの source adjoint と形式的な次数判定のみを確認し、新しい sigma jet、
source3/4係数、接続写像の値は一つも評価していない。

固定24単位逆の左右積、m0 n0=I、五つの実 sigma0 の cycle 恒等式を再構成した。
前段で残った C sigma²_rr に対する L=C n0 は

  L8=−43303/59400,
  L12=−149887/237600,
  L13=30713/59400

のみである。全168候補 path は零または実 B-sector の z0,z4,z5 Hadamard であり、
非零 path は108本である。原保存フィールド名の「168 nonzero」はこの候補全体を
指すものとして限定して読む。r³/r² で割る行7,15,23は L に入らない。

n=Esel(mEsel)⁻¹ と literal head section を固定すると、正次数の sigma 係数は
選択24行に支持される。m sigma=0 の一次・二次の係数を取ると一意に

  sigma1=−n0 m1 sigma0,
  sigma2=−n0(m2 sigma0+m1 sigma1)

となる。これは普通の Taylor 係数の式である。
n1=−n0 m1 n0 および n2 の逆行列展開を代入しても同じ式になり、n1 を零と仮定した
式ではない。m0 head が零という仮定も置いていない。

prefix の評価点は元の z=0 である。y=z+c r² とすると、その点では y=c r²。
最初の k 個の y を保持する prefix では、元 z_j は j<k で零、それ以外で−c_j r²となる。
差を c_i r² で割ると

  [r²]H_i(f)=[z_i r²]f−c_i[z_i²]f−Σ_(j>i)c_j[z_i z_j]f

を得る。この式を全495個の単項式で独立に確認した。
元の全次数四以上で z_i を含む項は除算後 r 次数三以上にしか寄与しないので、
この m2 部分に必要な関数次数は三以下である。adapted y=0 上の式と取り違えない。

実 sigma0 による五つの B 関数 contraction も再構成した。
一つは厳密に零、一つは nuA*bA67 のみ、残り三つは delta18/60/67、delta63、bB63 の
固定結合である。これは関数の型付き恒等式であり、その三次係数の値は未計算である。

V4 sigma0 の primitive frontier では、元 source alias から第一七 bAi=zi を確認した。
row63 の bA node17868は raw f63=node7790への実 alias で、全七 telescope により
z=0 上で gA=node505938 に一致する。こうして全56列を contraction し、未知四次関数は
delta18,delta60,delta67,delta72,delta63,gA,bB63 の七つとなる。
各 r³ell,r²ell²、計十四個は **V4 sigma0 部分だけ**の上界である。
tauAB の既知の−nuB2*bB2も保持した。

V3 sigma1、混合 V2 sigma1、L m2 sigma0、L m1 sigma1 を合わせた coherent contraction は
未完了である。特に sigma1_r に残り得る curve 行の source4は、この L の支持だけから
落とせない。最高次全体の相殺、既存 cache の全方向被覆、値 worker の実行許可、
元 H0・K torsion・nil 三乗の結論はいずれも主張しない。

再現 checker:

```sh
.venv/bin/python code/ambient_B18_C_connecting_sigma2_adjoint_source_review_20261105_v63.py
```

正値の事前計測を定数検査より前に保存した。次数四以下の literal 有理数式に対する
演算前上界657 bit、実際65 bit、所有数上界56,148、実行0.202秒だった。
60秒・2 GiB・150,000 slots・2048 bit・RAW32 MiB以内で、原資料は変更していない。
凍結出力先への再実行は拒否する。
