# v63 CM：正確な平方自由gcd・元接続係数・四生成上界

2026年11月6日（日本時間）。直前main `4ff9f6d4ad41f6bb0dd81b113f8d43fe902b9702`、Cloud開始SHA `4c3f623fb698f9adc6077015b25e2325b11072af`。対象と証拠を以下に限定する。新規性の独立調査は未完了。

## 二tailの標数零最大公約因子H全94項と平方自由性

元N0全4500項とN2全2317項に対し、固定四十素数の一回CRTとprimitive化によりH全94項を得た。二つの整数商は1795項と724項で、両整数割り算と元三変数の全係数掛け戻しが零である。独立実装はCRT・割り算・producer readerを用いず、逆順の整数積をdense三次元tensorの全10659/5070係数へ展開し、元4500/2317項と完全一致させた。Hを法32003でmonic化すると既存の厳密gcd3899へ一致し、両商のleading係数は生き残る。既存の独立modular商互いに素から、HはQQ[a,e,d]で正確なgcdとなる。Hの各変数次数は6,9,4、全次数11。別の同じ固定特殊化e=0,d=1で、a次数六を保ち、独立なderivative Euclidによりfとfprimeが互いに素と確認した。Hのa-leading=17578125*d³*(80-25d+66ed)とa0のe0係数からa-contentは単元なので、H自身も平方自由である。Hは六因子base環でも非単元で、選択二tail idealの次数一codimension部分はこの被約な主因子H、残る(A0,A2)はproperならheight二以上である。

- Hの既約性は未決。指定の固定六次modular像は可約であり、これは有理数上の可約性を意味しない。
- 六因子base環での選択二tail idealの記述に限る。HとD12の関係、実bad Fitting、全FULLU57/239、補集合消尽は未決。

証拠：

- [`reports/S61_TWO_TAIL_CRT_LIFT_SOURCE_20261105_v63/report.json`](../../reports/S61_TWO_TAIL_CRT_LIFT_SOURCE_20261105_v63/report.json)
- [`reports/S61_TWO_TAIL_CRT_LIFT_WORKER_CUSTODY_20261105_v63/report.json`](../../reports/S61_TWO_TAIL_CRT_LIFT_WORKER_CUSTODY_20261105_v63/report.json)
- [`reports/S61_TWO_TAIL_CRT_LIFT_ROOT_SOURCE_REVIEW_20261105_v63/report.json`](../../reports/S61_TWO_TAIL_CRT_LIFT_ROOT_SOURCE_REVIEW_20261105_v63/report.json)
- [`reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63/report.json`](../../reports/S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63/report.json)
- [`reports/S61_TWO_TAIL_CRT_ROOT_DENSE_REVIEW_20261105_v63/report.json`](../../reports/S61_TWO_TAIL_CRT_ROOT_DENSE_REVIEW_20261105_v63/report.json)
- [`reports/S61_TWO_TAIL_COMMON_FACTOR_SOURCE_20261105_v63/report.json`](../../reports/S61_TWO_TAIL_COMMON_FACTOR_SOURCE_20261105_v63/report.json)
- [`reports/S61_TWO_TAIL_H_ROOT_SQUAREFREE_20261105_v63/report.json`](../../reports/S61_TWO_TAIL_H_ROOT_SQUAREFREE_20261105_v63/report.json)

## 元接続写像二headの四高位係数とKbar基底parameterの確定

同じ元接続写像のHEAD29/30を既存Phi三次以下cacheから一回計算し、別の疎radix5多項式代数・逆順積・geometric inverseで全8081 node hash、safe delta63全B63、八bBのtrusted plane10、全2673 post積と八出力を独立確認した。canonical高位は(p2,b2)=(13621/40344,-1573/1920)、(p3,b3)=(-32129/161376,336618121/20528640)。HEAD30の汎用式に現れるdelta63(z1*ell)は元C2slot21の既証零だけで除き、全十四nonzz二次も独立に零と確認した。未知Phiのsocle係数を真の零とせず、zero-gradientによる全deltaの安全性を用いる。既証のhonest Kbar七基底I4+rh3+U+Sに必要な四parameterは全て実値になった。

- 新しいPhi・IFT・gauge・高次source solveは無い。head0/head4の四highは未計算だが、このhonest七基底には不要。
- 実Kbar作用表、元Kの全形式lift、元H0長・ell指数・nil三乗は別課題。

証拠：

- [`reports/AMBIENT_B18_C_CONNECTING_HEADS2930_CACHE_SOURCE_V2_20261105_v63/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_HEADS2930_CACHE_SOURCE_V2_20261105_v63/report.json)
- [`reports/AMBIENT_B18_C_CONNECTING_HEADS2930_CACHE_SOURCE_CODE_20261105_v63_REVIEW/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_HEADS2930_CACHE_SOURCE_CODE_20261105_v63_REVIEW/report.json)
- [`reports/AMBIENT_B18_C_HEADS2930_ROOT_SOURCE_REVIEW_20261105_v63/report.json`](../../reports/AMBIENT_B18_C_HEADS2930_ROOT_SOURCE_REVIEW_20261105_v63/report.json)
- [`reports/AMBIENT_B18_C_CONNECTING_HEADS2930_CACHED_VALUES_20261105_v63/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_HEADS2930_CACHED_VALUES_20261105_v63/report.json)
- [`reports/AMBIENT_B18_C_CONNECTING_HEADS2930_CACHED_VALUES_20261105_v63_REVIEW/report.json`](../../reports/AMBIENT_B18_C_CONNECTING_HEADS2930_CACHED_VALUES_20261105_v63_REVIEW/report.json)

## 元比較kernelのR0加群生成数は四以下

honest Kbar=K/ellKの七基底を(v,q1,q2,g,w,U,S)、v=z5h4,q1=rv,q2=r²h4,g=rq2,w=rh3とする。rh1は実際にITold=span(v,q1,q2,g)へ入り、rU=w-a*rh1-b*q2である。従って四元U,S,v,q2とそのr倍が七基底全てを生成する。独立証明は列(U,S,v,q2,rU,rv,rq2)の七次minorを取り、rUのw係数一以外を全て任意の未知係数のまま残した。selector行順[5,6,0,2,4,1,3]の行列式は-1。元KのD有限性とell-Nakayamaから、これら四bar類の任意の実K liftが元KをR0加群として生成する。

- 四という上界であり最小性は未証明。D上の最小生成数七、nil ideal全体の最小五生成とは別の主張。
- 元Kへの形式liftは存在だけを示し、全係数を展開していない。元H0やnil指数の結論は無い。

証拠：

- [`reports/AMBIENT_B18_C_K_R0_FOUR_GENERATOR_UPPER_SOURCE_20261105_v63/report.json`](../../reports/AMBIENT_B18_C_K_R0_FOUR_GENERATOR_UPPER_SOURCE_20261105_v63/report.json)
- [`reports/AMBIENT_B18_K_FOUR_ROOT_MINOR_REVIEW_20261105_v63/report.json`](../../reports/AMBIENT_B18_K_FOUR_ROOT_MINOR_REVIEW_20261105_v63/report.json)

## 検証と未完部分

CRTは四十素数全て受理、10.21秒、field積2638783719、bigint積632130、最大owned52055、RAW1137878 byteであった。candidate後の条件付きbit上界は664/555で8192以内。独立dense法は236786積、実中間293bit、source300/276、failure複写込み56966所有、約0.23秒。原二tailは互いに素ではないが、その整数商は互いに素である。H既約性の固定modular検査はnegativeで、同じ像のsquarefree性は別のderivative GCDで独立証明した。別素数・再特殊化・factor familyは追加していない。

接続二headの一回workerは3.67秒、実62bit/演算前125bit、owned67465、RAW744557 byte。独立法は0.553秒、tracked72459であり、全cacheと全termwise出力が一致した。SOURCE上界root51/103・post118/237、hard384/1153を維持した。V1で見つかったshared odd root38648の63→plane10射影忘れは、signed実行前に履歴とcodeを保存し、別V2にした。汎用HEAD30式のdelta_z1ellの扱いも、既存零をpinした補助証明として凍結資料を改変せず保存した。

次のT2 case554はnative四列のweighted S零次を保ち、非零outer644とOLD/FIRSTも全て残す。nativeの別C72 algorithmはcoarse339/1018、収縮後のfine envelopeとOLD等を合計するwhole予算・worker・codepeerは準備中で、signed554値は未計算。全605/617やT2 HF・全長は未決。

一般五点の次のrow17/col79はfull F11を含むcleared minorの特定五座標係数を目標とする。SOURCE compilerが157050>150000のfinal ownership guardで停止し、全bridge/依存計算のprefixを保存した。supportやsigned値は未実行。次はhash-bound prefixから不要parse文書を退役させる別finalization SOURCEを設計し、完了分を再計算しない。

HとD12の比較は旧axis determinantを四次商で一回計算し、一次余り[27788,575,7510,2387]を得た。独立検証は未完了。Kbarの七基底の具体的組込みと初期r作用を既存関係から整理するが、高次作用を零補完しない。元H0・nil三乗・正確な指数、T2長17〜20、全FULLU57/239・全補集合・一般五点尽くし・全分類上界は未決。元nil指数3〜6、下界58（rank4は52）、NEW式28次と五凍結原稿を維持し、保存後も実研究を継続する。

実Kbarの組込み11x7と全11次QQ基底は `AMBIENT_B18_C_KBAR_ACTUAL_SEVEN_BASIS_SOURCE_20261105_v63` および独立review `a5ec34be…` で確認した。全28edge恒等式は残る四highを形式のまま保ち、双方121逆行列積を点検した。typed56x7 cycleはfactored source表現である。QQ線形補空間を元R0 UCTの分裂とは扱わない。

追記：四生成上界は、既存初期関係から得た非零lambdaSにより三生成へ改善された。最小R0生成数は二〜三。独立な商消去は `AMBIENT_B18_K_THREE_ROOT_QUOTIENT_REVIEW_20261105_v63`（aab52381…）にあり、primary `e4d6acd7…` のpolynomial minorを再利用していない。rank r=4でJordanは331/322、深い一scalarだけが未計算である。元ell-DVRやH0への結論は無い。

D12補足：独立Newton/Faddeev算法の全十二段とCayley-Hamilton零行列を記録し、四次商の余り(27788,575,7510,2387)が一致した。H∤D12を六因子baseでも証明した。もしV(H)∩D(U6D12)が空なら、H平方自由よりU6D12∈sqrt(H)=(H)、coprime unitを消してH|D12となり矛盾する。従ってあるH成分の二次元開集合がこの部分に残る。Jacobiの実rankや実bad Fittingの証明ではない。証拠は `S61_TWO_TAIL_H_D12_ROOT_NEWTON_REVIEW_20261105_v63` と `S61_TWO_TAIL_H_D12_NONEMPTY_GAP_20261105_v63`。case554 primaryは両upper零、独立値確認待ち。
