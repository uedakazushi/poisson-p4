# 真正 N/ellN の有限三次積と全次数 lift の区別

2026-11-05。これは次の計算に使う source 方針で、積の値や元 nil 三乗の消滅を主張しない。N=Nil(D0)、Tnil=H1(元 radical 比較⊗D k) とする。既証の正確列は

    0→N/ellN→Tnil→Cprime[ell]→0

で、次元は14、16、2である。真正 connecting の kernel として十四個の bar-cycle を選べれば、それぞれは何らかの元 N の元へ持ち上がる。全次数 lift の明示式をまだ求めていなくても、この存在は正確列によって保証される。

受理済み product44/cube44 は、元係数環 R0 上で実 N の cycle representative の積を与える多項式式である。従って、真正 kernel の bar-cycle にこの式を ell 零化して適用すれば、未知の ell 倍の lift 補正は結果に影響しない。結果の Tnil class は、実際の N/ellN の積の自然な像となる。源の代表を変えた場合も、その class が変わらないことは actual lift の存在と元 N の積から従う。

ここから、**真正 kernel の有限三次積を調べるためだけなら、全十四個の全次数 lift を先に展開する必要はない**。ただし、bar-cycle と Tnil の実 source projection、出力の境界判定はなお必要である。全 HF55321 と filtered 十六基底だけでは全作用表・quotient reduction は揃わないため、これらを零埋めしてはならない。

Tnil の全十六 classes 自体に代数構造があるとは主張しない。元 product44 の boundary descent は元 H1 と D0 の同型を使い、ell 零化前の全780成分 homotopy を与えていない。従って、持上がらない Tor class にまで積が well-defined とする近道は使えない。計算対象を真正十四次元 kernel に限る。

非零の三次積が得られれば、無条件で元 N³ は非零である。全三次積が零なら N³⊂ellN までは従うが、消滅はなお未証明である。追加仮定 ell H0(D0)=0（実 T2 全長17と同値）の下では H0∩ellN=0 と N³⊂H0 から、元 N³=0 まで従う。標数零では十四基底に関する高々560個の対称三次係数を使える。

次は genuine edge 十値の source と有限 target を確定した後、cube44 の全次数源から必要な有限係数・Tnil projection の依存を抽出する。新しい physical source、全 mixed jet、全作用表の計算を自動で始めず、既存 cache と重なる functionals を先に調べる。この source 方針の独立監査は次の引継ぎ課題である。
