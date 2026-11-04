# 共通必要曲線の因子・還元入力・L成分の独立監査

2026-11-04。対象は `RAW_H_B_COMMON_CURVE_FACTOR`、`RAW_H_B_CURVE_COMPONENTS_INPUT`、`RAW_H_B_COMPONENT_COMMON_ROOT_L` の三凍結packetである。新しいfactorization、Euclid/GCD、点、素数、GB、著者workerは実行しない。上流の実minorとpair monic/normの独立reviewへhashを結ぶ。

保存norm pencilの最終49項Gとfactor packetのsource Gを全係数で比較し、Q*L^4への全掛け戻しを独立QQ辞書演算で確認する。Qはaモニック2次、L=a+(9e²+3e−6)/50。Qのa-discriminantは(32481e²+30546e+8881)/810000で、このe二次式のdiscriminantは非零である。従って平方因子を持たずQQ(e)の平方ではないため、QはQQ(e)[a]で既約である。QとLが必要曲線の二成分であることと、実bad locusの成分であることを区別する。

実gは22 minor generatorのindex15、kはindex20であることを確認する。Hはmonic d4で、a^5Hを係数多項式に戻し、元f,gに対する全係数所属恒等式も再確認する。反転するのはaだけである。Q/Lそれぞれについて、保存された全六quotient/remainderを掛け戻し、a次数の余り条件も確認する。新しい除算は行わない。

L成分では係数を未約分のQQ[e]分子・分母の対として扱い、等式は分母の交差掛けで判定する。保存されたEuclidの全段で、LCと逆元の積1、quotient/remainder恒等式、余り次数減少を確認する。各traceの二source Bézout、最終三source Bézout、全三sourceの共通因子による割り切れを独立に検算する。

その結果、Lのgenericな共通d因子は正確に

`(d+20/(3e−2))²`

である。元のFULLU productをLへ直接代入した結果と保存Uを比較し、全U support除去の掛け戻しも確認する。さらに実Bの式を代入すると

`B|L=(3e−2)*((3e−2)d+20)`

となる。従って分母例外を除いたL上の三本H,g,kの共通根はB=0にあり、FULLUの外である。

serializedされた全field係数の非定数分母を独立に収集し、保存13多項式・次数合計上界65と一致することを確認する。その零点を単位扱いして捨てない。L上でaはeの多項式であり、各例外e値に対してaは一意である。FULLUではa≠0なので、a^5Hのd4最高係数は非零である。よって例外に残り得る実候補は有限点であり、空性や実在は未判定である。

このreviewはQ成分の共通根解析を含まない。Qの入力還元と既約性だけを確認し、新しいQ exact quotient packetは別gateとする。特にQを除去済みとせず、全bad locusが有限だとも主張しない。旧norm曲線を除いたh/G,k/G由来の有限残余も保持する。

Gの「必要曲線」という上流の意味は元norm pencilの証明を前提とし、今回新たに全norm pencilを再計算したとはしない。全FULLUでのB rank3、H7、cubic57、quartic239、尽くしへの昇格はない。専用REVIEW reportに全input/evidence hashと実測値を保存する。
