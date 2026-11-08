# Concurrent m二normal-double全源とsmooth cubic scrollのT37帰属：独立検証

2026-11-08、別作者 quartic_review。PRIMARY E。**PASS**。
producerをimportせず、既存の独立標準ライブラリFraction／Laurent多項式算術で、
保存した同じsource・canonical・normal・graph・実frameだけを再構成した。
新しいD/operator/M1/anchor/profile/jet又はCAS計算を追加していない。

## 全sourceと実正逆

既独立の十七whole pencil／84複素wallの尽くしは固定入力として用いる。
今回の三残源は、旧full150 P0、元complete K1、**FULL9 affine P2**から全150係数を再構成した。
元J2の非零九次minorとnullity零により、保存particularは全P2解を尽くす。
全350元式とsigned Euler-Pfaffianを直接展開した。

h21のa=−4/11では、局所化前の元row34は−m²/10である。
m=0の全Dはuで割れるので、被約既約quinticのsourceは空。
h3のa=5/2、−1/2では、全八係数の全値で全350式が零となる。
六shearと(r,m)への定数可逆写像は、全八係数の明示逆を直接再構成した。
有限shearの正逆、Euler補正、全九P2 particularを同時に戻す恒等式も確認した。

P0保存unipotentは全εで可逆で、正規化後の作用は
r→r−5εm/2、m→mである。m≠0ではε=2r/(5m)と
diag(1,1,1,m,m)がcanonical r=0,m=1へ送る。
ここだけmを反転し、その補集合m=0は全Dのu因子として保持した。
二canonical間のRAW写像
(a,b,c,u,v)→(−a,b−a,c/3,u,v)、tensor scalar 1/3も全150係数で一致した。
したがって被約既約Dを持つ両whole sourceは一つの実軌道へ戻り、
その像閉包はcanonical軌道閉包と等しい。

## 同じ有限graphとcompleted T37への実恒等式

保存済みD／Euler／実群行列について、rank101／24／48の非零minor、
全核の全行消滅とfree-coordinate単位行列を検査した。
これは既存行列の証明書再検算であり、新Dを計算したものではない。
同じNは群の全51列を消すfunctionalで−59/24を与え、characterは−2。
実48群列、同じU101、Nの全150直和行列の行列式も非零だった。

座標重み(0,0,1,−2,−1)とtensor scale −1から全150／350重みを確認した。
保存graphの各layerは同じU内にあり、degree kの全係数はcharacter −2k。
graph degreeは三、元350式の全多項式はdegree零から六まで恒等零。
t=15552eの整数e一tensorも全150係数で一致した。
先験的な有限性／formal IFTの論理も読解したが、以下のmembershipには
全多項式族とCstar作用だけで十分である。

θ²=−118の商体をFractionで実装した。完成済みT37 JSONのall-one源を、
七つのhigher correctionを省かずに読み込んだ。
保存FとF逆について左右積を検査し、detF=−1/167510398464を得た。
λ=θ/50976、保存した全一次場Rを使い、

```
F_*T = λ P(e=1) + E∧R
```

の全150係数とその逆写像を直接確認した。
さらに全signed quinticについて
λ²D_P(FX)=det(F)D_T(X)の全126係数が一致した。
数値近似／有限体／normal blockだけの一致ではなく、
θの両embeddingに有効な商環内の恒等式である。

e≠0は同じCstar軌道に入り、e一が閉PGL不変成分T37に入るので、
実多項式族の閉性からe零もT37に入る。
全sourceの明示逆とm=0補集合、既十七pencil／全wall排除を戻し、
**reduced concurrent R・rank二plane・m二normal-doubleの被約既約quintic全点はT37**と結論する。
親成分の局所一意性をこのmembershipの必要仮定としていない。

## Whole A25からsmooth cubic scrollの全点へ

先行の `CUBIC_SCROLL_A25_WHOLE_RESIDUE_SEPARATE_REVIEW.md` は、
全normal-root還元を戻したA25についてnode全空、cuspは一実軌道又は空を独立認証した。
今回、その軌道の実在と名前を与えるcanonical Aを直接検査した。

Q=9au−2v²、C=16cu²+9(a−4b)uv−2v³について
486uF=−2C²−Q³を再展開した。
u=1ではF=−256/243[c+9(a−4b)v/16−v³/8]²−3/2[a−2v²/9]³であり、
括弧内は(a,c)のregular座標である。
全minors行列[[u,2v,9b],[v,9a,4c]]は可逆な尺度変更と変数順序で
standard smooth S(1,2)に一致する。
そのdense chartへの全五gradientのpullbackは零、intrinsic affine bracketは
{b,v}=4(9b−v²)/27≠0である。追加のSing planeを仮定で消していない。

Fはcのquadraticで、先頭係数は−256u³/243、定数係数はuで割れないのでprimitive。
discriminantのv=0 restrictionは−512a³u⁵/81で非平方であり、
全discriminantが平方ならこの非零restrictionも平方になることに反する。
従ってFはC上絶対既約、tensorは実rank四、smooth scroll上一般rank二かつcuspである。
canonical-root等の既排除を戻せばgenuine A25 cusp seedであり、上の実T37帰属を持つ。

全cusp sourceが一実軌道へ戻る先行定理と閉PGL不変性から、
**smooth cubic scrollを非正規支持に持ち、その上で一般rank二となる被約既約quintic全点もT37**。
clean、R被約、normalization-flat、又はgeneric四full-P2 fiberを追加仮定していない。
rank零scroll、cubic cone、他導手次数や全分類の完了は主張しない。新成分増分は零。

## 再現と失敗記録

実行前gateしたV2 jobは180秒／2GiB。
`CLASSIFICATION_PLANE_M2_CONCURRENT_T37_SEPARATE_REVIEW_V2_20261108.json` は
**1.397922秒、95,896KiB、ALL PASS**。
inline source、全input hash、既150直和det、全source判定を保存した。
V1は最後の保存quintic readbackで、producerのH,X,Y,Z,W labelをx0,…,x4と
取り違えたparserにより停止した。全150正逆まで成功していたが、V1をPASSと呼ばない。
V1 source／jobと全stdout／stderrを別failure receiptに不変保存し、
V2はそのlabel parserとartifact identityだけを訂正して**同じ全入力・同じ点**を再生した。
producerや旧証明書のbytesは変更していない。
