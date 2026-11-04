# CA-critical の第三embedded閉包：独立幾何・Ass監査

2026-11-04。新しい係数・jet・source・点・CASは計算しない。
独立監査済みの全proper旗模型、V1の実rank20・monic critical曲線、同じR0 conductor、
全接触完成環、Pa/Qbの全fiberおよびfull D theoremをhashで固定し、適用の論理を監査する。

P0の**全**proper fiberが一台長2なので、proper mapはその近傍で有限となる。
有限algebraのfiberのvector-space次元は上半連続で、全fiber長≤2の開を取れる。
有限性・fiber長は同変な性質なので、その最大の条件開、または群translateの合併を使うと
G-stable開になる。P0が20次元critical軌道の閉包にあるため、この開はその軌道に会う。
軌道の推移性により全軌道が開に入る。CA選択liftのcorank1は局所fiber長≥2を与え、
全fiber長≤2と合わせて追加の旗やendpointを排除する。Ac8の完成像から正確にdual numbersである。

finite smooth sourceを初めからnormalizationとはしない。既受理の完成像
`R+I*tau`で非零Iを取り、`tau=(f*tau)/f`から完成分数体の一致を得る。
全fiber一台によりsource完成の別factorはなく、faithfully flat completionは有限写像の
一般次数を保つ。従って既約source/targetの一般次数が1であり、有限な正規sourceが正規化となる。

CB側は選択点の微分rankだけで全fiber1とは結論しない。実canonical conductor商は
同じphysical R0の `R0/I_B0`。accepted CA-critical sourceのmonic曲線と同変G19切断が、
`Z_C`の実formal branchと `V(I_A0)` を結ぶ。`I_A0+I_B0`がm-primaryなので、
このbranchの一般点ではopposite conductorが単元である。
もしZ_CがCB conductorに含まれれば、そのcoherent idealの完成・同branchへの制限も零となり矛盾する。
conductor外の開はG-stableで軌道に会うので全軌道を含む。正規化はconductor外で同型となり、
正則sourceからCBはsmooth29、**全fiber**は被約一点と従う。

独立のmonomial arithmeticで、全actual完成環

\[
 k[[g_{20},l_7,v,x]]/(x l_7,x^2,xv^2)
\]

の非零class xvについて `Ann(xv)=(l7,v,x)` を確認する。また
`x(l7,x,v²)=(x)∩(l7,x²,v²)` を生成元単位で確認する。第二primary商は長4、
nil moduleの長2と混同しない。この計算は有限jetからの全ideal推定ではない。

accepted simultaneous G20作用とambient sliceの形式IFTにより、軌道は群factor×slice原点である。
従って上のassociated primeは、未同定のformal-fiber primeではなく、実代数的軌道idealの完成である。
faithfully flat completionのAss contraction、続いてlocalizationのAss公式により、
閉包Z_Cのideal自体がglobal actual Dのassociated primeとなる。
Z_CはN28のproper20次元部分でありembedded。R20とは、`R∩N`が19次元なのに
Z_Cの20次元軌道がN内にあることからも区別できる。非零scalar不変の同次閉包なので射影次元は19。

Pa、Qbの実G×非空unit-family像も既約20次元閉包にconstructible denseである。
同じ閉包なら各dense constructible集合の含む非空開が交わる。
しかしcriticalのwhole fiber patternはCA一点長2・CB一点被約、PaはCA二被約点・CB一被約点、
QbはCA一被約点・CB二被約点であり矛盾する。special limitで一点長2へ退化できることは
このgeneric-open比較を妨げない。

以上より、actual CA∩CBは少なくとも三つの相異なるglobal embedded associated closure
（各cone20/projective19）を持つ。全Ass list、原点全Ass、閉包境界、Poisson成分追加や尽くしは主張しない。
専用checker `code/ambient_CA_critical_third_embedded_closure_review_20261104_v63.py` と同名REVIEWに保存する。
