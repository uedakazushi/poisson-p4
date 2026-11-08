# concurrent m二の一実軌道：次数三の全IFT切片が先験的に尽きる

2026-11-08、PRIMARY E。同じ強制canonical点一つだけを用いる。後続のliteral graph算術の前に、全formal supportの有限性と成功／失敗の含意を証明する。

## 入力と真正切片

whole A8二源の実normalizerと相互RAW同定から、対象の全irreducible-quintic点は一つの実PGL軌道閉包へ入る。そのcanonicalはperiod−1/2、r=0,m=1、整数scale54である。新しい点を選ばない。

保存fullD350×150のrankは101、Euler25のrank24、Euler/GL/scalar全51列のrank48、全DMは零である。Dの保存nonzero101minorの座標列をU、その行をRとする。各座標は実トーラスの固有vectorである。全kernel49内のweight−2成分から、群48列空間に入らないhomogeneous Nを選ぶ。これは保存48minor inverseによる射影残差と、全51群列を殺すlinear functional φ、φ(N)≠0で証明する。新しい150minorは不要である。

D(U)はrank101、D(N)=D(G)=0、NはG外なので

```
C150=G48 ⊕ U101 ⊕ C N.
```

従ってS=P0+U+CNは真の102次元線形sliceである。Euler quotientやchosen tangent subspaceからfull sliceを推測しない。

## 全formal graphの先験次数

normalized actual actionは

```
s^{-1} diag(1,1,s,s^{-2},s^{-1})_*
```

であり、全150点係数、Dの全350×150非零entry、全51実群列が保存weightと一致する。Nのcharacterは−2、全input character範囲は[−6,4]である。

selected101式のimplicit微分D[R,U]はinvertibleなので、形式陰関数定理はU座標をtNに対して一意のformal graph f_i(t)へ解く。f_i(0)=f_i′(0)=0。slice及びselected式はトーラスで保たれるので一意性から

```
f_i(s^{-2}t)=s^{w_i} f_i(t).
```

従ってt^j係数はw_i=−2jの座標だけに支えを持ち、j≤3である。**全formal graphそのものがdegree三のliteral polynomial**である。j2/j3を一回ずつ解き、以後のjet次数を追加する必要も根拠も無い。

元全350式はtensor係数にquadraticなので、この全graphを代入した式のt次数は高々六である。全350行・次数零〜六のliteral係数を保持して検査すれば、未確認higher formal termも未確認higher equationも存在しない。

## 全ambient supportへの逆橋

48本の実Euler/GL/scalar方向をそれぞれの連結代数群作用の局所parameterで取り、その積と102sliceから全raw150空間への作用写像を作る。微分は上の真の直和からinvertibleなので、これは原点近傍でétaleである。Poisson idealはこれらの実作用及びEuler gaugeに不変であるから、全Poisson schemeの局所completionは48smooth group directionsと全slice completionの積となる。closed orbit又はLuna sliceは仮定しない。

全350係数が零ならslice completionはC[[t]]であり、raw schemeはsmooth unique49、Euler24とprojective scale一を除くprojective schemeはsmooth unique24となる。全graphは実多項式Poisson曲線である。正規化トーラスt→s^{-2}tにより全t≠0は一actual PGL軌道となり、群48＋Nの微分下界からそのprojective orbitは24dimで親にdenseである。元の23軌道と全A8閉像はこの閉24親に含まれる。親の名前をこの有限graphだけから決めない。

逆に元350残式の一つが非零なら、weightにより各残式は一つのt冪の定数倍である。全sliceのreduced supportは原点だけになる。étale積から全raw局所reduced dimensionは48、元23実軌道の下界と一致し、その軌道閉包は全ambient23成分となる。厚みやnil指数は求めず、scheme smoothとは主張しない。既129との相違は別義務である。

本jobは同じ点・旧101/48minor・全kernelのreadbackだけを使い、適応的jet、第二点、別profileを追加しない。120秒／2GiBを上限とする。
