# T37一全rational plane枝：全因子readbackと二H+cubicへの閉必要障害

2026-11-08、PRIMARY E、`/root/quartic_sources`。root明示GO・fresh gate後、凍結した一つの全u∈Gm源だけを一度実行した。真T37の新cell、kernel、operator、point、GB、profile、jetを追加していない。結果はこのuniform入口の陰性であり、T37全成分の非所属を主張しない。

## 1. 全parameterのliteral因子証明書

`TWO_LINES_M3_T37_RATIONAL_PLANE_ENTRANCE.md` の全六a-mapを使う。全元Lagrangian六式の代入、真completed tensorからの全126 signed Euler-Pfaffianと旧literal fullDの一致、whole-cell restricted quinticへの引戻し全21係数、最後のfactor積を再構成した。

SymPyのQQ(u)厳密因数分解は候補発見だけに使い、各factorをQQ[u][A,B,C]のprimitive形へ直した。全21積係数は別Fraction/Laurent積で掛け戻した。nonmonomial parameterは反転していない。45秒/768MiBの同一有限契約に対し、実測0.339738969秒、最大RSS49,700KiB、全PASS。

証明書は `reports/CLASSIFICATION_TWO_LINES_M3_T37_RATIONAL_PLANE_ENTRANCE_20261108.json`、SHA256 `2805d5d081f16161aa5cd69f7ddeb3207d2f69445046cf2ed2095824d5a79bf1`。実行inline全bytes、全入力SHA、実plane、全polynomial、primitive内容、全係数掛戻しを保存した。有限契約cefa668eの「未実行」作成時bytesは変更しない。

正確な結果は

```
D|plane=Q*G/(31104u⁴),
Q=3uA²+2B²−12uBC+72u²C²,
G=216u³A³−(15552u⁵+1)A²B+(46656u⁶+12u)A²C
   +144u²AB²−864u³ABC+5184u⁴AC²−5184u⁴B³+144u²C³.
```

両factorのu係数内容は1で、全uに対し零多項式になることはない。scalarはGm上のunitである。Qのsymmetric matrix determinantは **324u³**。従ってこのcell内の全u≠0でQはsmooth irreducible conicであり、conic退化parameterを一つも捨てていない。u=0は元cell sourceに存在しない。Grassmannianの他cell境界は、この式へu=0を代入しただけでは記述できない。

歴史的chatの同じ数値は未確認としていたが、今回初めて全literal入力から厳密に確認した。旧byteは書き換えない。

## 2. 手で導いたwhole閉必要条件

ell=B−6uC、kappa=1+7776u⁵と置く。上の全八項Gを係数ごとに整理すると

```
G=Q*L + (2*kappa/(3u))*ell³,
L=72u²A−(5184u⁴+1/(3u))*B+(15552u⁵+4)*C.          (1)
```

である。例えばA²BとA²CはQ*Lからそれぞれ−15552u⁵−1、46656u⁶+12uを得る。Aを含まない残四係数は順に

```
B³: 2*kappa/(3u),
B²C: −12*kappa,
BC²: 72u*kappa,
C³: −144u²*kappa
```

であり、(1)のcubeを正確に再構成する。これは新CAS計算ではなく、凍結factorの手係数恒等式である。kappa=0も保持しており、その場合G=QLでD=Q²Lとなる。

このcube条件を退化に対して用いるには、存在するLを不用意に有界と仮定しない。五変数の斉次二次Q≠0について、R3/(Q*R1)は一定rankのvector bundleである。Q*R1は五次元で、Q≠0の全域でmultiplicationがinjectiveだからである。projective cubic G及びprojective linear ellに対し、商における[G]と[ell³]が一次従属というrank条件は閉である。ellを忘れる射影はproperなので、ellの存在も閉条件となる。Q*Gへのprojective multiplicationもproperである。

この閉条件はell³がQ倍となる零商も含む。それを排除した一般点の条件として使用しない。一方、目標のQ=h*xではh,xが相異なるため、Qが一つのcube ell³を割ることは不可能である。従って目標で商[ell³]は必ず非零で、cube条件は正確に

```
C_target = h*x*L0 + c*ell0³.                    (2)
```

へ戻る。既約cubicなのでc≠0。h,xがcubic因子ではないことから、h=0とx=0への制限は両方とも非零のlinear cubeになる。

## 3. 第一二targetへの適用と正確な停止範囲

対象Dは二つのsimple Hと一つの既約cubicの積である。degree二factor Qはその二Hの積以外へ移れず、Gはその既約cubicへ移る。factorの並べ替えによる例外は無い。従ってproper factor flagを通して前節の閉必要条件を適用できる。

第一targetは C=h*y²−x*(x−h)² であり、C|x0=h*y² はcubeでない。第二targetは C=h*y²−(x+h)³ で、C|x0=h*(y−h)*(y+h) もcubeでない。どちらも(2)に反する。この議論はfactorやPoisson tensorへ施すPGL変更に不変である。磁気六係数はactual Dを変えない。

従って、**今回の全rational-plane normal-cone端点族のPGL飽和閉包は、第一二targetのmarked二H+cubicを供給しない**。uの全零／unit境界をprojective factor incidenceで扱っており、固定uの数点の比較ではない。現cellのu=0を架空のsource点に足してもいない。

これは、まずnormal-cone端点を取り、その全族のPGL閉包を考える入口の障害である。元T37 tensorの高次normal項を同時に残す非一様u(e)/period/moving-frame曲線の障害ではない。T37全component非所属も、新成分性も結論しない。このuniform routeはここで停止し、別cell・別parameter点のprofile列へ拡張しない。全lambda源へのpositive actual橋は依然未決である。
