# ONE549 native frame と零 envelope：独立 SOURCE 検証

独立検証は PASS した。報告 SHA256 は
`b02f0d2118aab1499477d6c05603e8e0a7f8738b83cc702795fd96224a234e18`。
frame source `6fe41836…` とpositive source `62d66bc4…` を固定した。
符号付きnative応答のjobは実行していない。

元case619の549番は source768の `298/4+649` であり、
563/609による追加基底補正は零である。full a1の全32項を元sourceと照合した。
native649はKoszul source2、input mask1であり、
Sの二項は row0のG3とrow3の−G0である。
G3の原点値は零、G0の原点値はtauA*tauB/2。
mask1=tauAを掛けるとnuA*tauB/2となり、nuAの原点値が零なので、
**full48成分のS0列9全体**が零になる。
隣接する648はmask0で、S0列8のindex15係数−1/2が非零である。
両者を取り違えず、weighted covectorだけが消すという弱い条件にも置き換えなかった。

J18はcanonical B17にr²ellを加えた18次元の同じR0 quotientである。
全660個の非斉次normal formと全1,089積を独立に読み直し、
元raw33から指定monomial quotientへの写像を確認した。
このためzを含む未知N二次係数がpureplane r²ellへ戻ることはない。
全monic対積がliteral支持だけで零になることも確認した。

canonical C68とC28の既存G8/nu、元Nのpoint・全必要first・rr/reを、
元B11/C16/first cacheと照合した。Nはsource12×physical12であり、
転置はこの十二行列だけに行い、mu4は転置しない。
fullGamma56/D56を保持し、native出力の下八行零だけを使う。
他のsourceに対してbottom成分を捨てる根拠にはしない。

Dの二次係数では未知L2はfull S0列9を掛けるため零である。
Gammaのpure r²ellでは未知L3も同じS0を掛け、未知nonplane L2はz因子を持つため寄与しない。
既知Lrr/Lre内のN0*nu二次等の全trace carryと、ordinary rrの一回だけのfirst積を残した。
N2やN3そのものを零とは置いていない。

独立計算はproducerの正D/U算術と異なり、Boolean semiringで各経路の存在だけを伝播した。
18個のintrinsic S支持と68個のcovector*L支持を、元係数の厳密な零・非零から再構成し、
保存されたpositive cacheの全支持と照合した。
Gamma三次の二channel×五項と、D二次の全項を合わせた**68個すべて**で、
dotの二支持の交叉が空だった。
したがって両native合計の正envelopeは厳密に `(D,U)=(1,0)` であり、
符号の相殺や新しい応答値を計算しなくてもnative部分の零を証明できる。

Gamma pureplaneに入るa1はpivot648のr/ell係数二項だけである。
648のrr/re・両channelの四応答値は、親担当の独立92検証10844d87の零を固定して用いた。
648はS0非零なので、その個別L2S0項を省略して導いた零ではない。

このレビューの結論はnative649部分と既知648のplane couplingに限る。
OLD298/4と全55 FIRST carryの別source監査を合わせる前に、
whole ONE549や全549/617 upperの零を主張しない。
新しいN2/N3、physical/gauge jet、HF・長さ・元H0の値も求めていない。

独立SOURCE所有上界76,461、実行2.406秒、55秒・2 GB・RAW32 MiB内で終了した。
再現コードは `code/ambient_B18_native_one549_source_positive_review_20261105_v63.py`。
