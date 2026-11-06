# 独立三超平面＋既約二次曲面：全源と四有限階数の独立確認

2026-11-07。課題 `E_reduced_nonnormal`。
対象は `THREE_HYPERPLANES_QUADRIC_INDEPENDENT.md` と、四つの零対数括弧源の算術証明である。
主原稿・制御ファイルを変更せず、手計算の確認と再現した厳密計算を区別して記す。

## 1. 非零対数括弧の七配置：手計算で受理

同ノートの記号で、`a≠0` の全源は

    P=aD1∧D2+bQ∂z∧∂w,  ab≠0,
    [D1,D2]=0,  Di(Q)=tr_vertical(Di) Q

である。Pfaffian条件から係数を消去してこの表示を得る逆方向も成立する。
`h,x,y` の joint character は相異なる `0,e1,e2` なので、radicalへの非零restrictionも独立である。
この事実、非退化商のtrace公式、逆二次形式上のcharacter和の条件を使えば、
同ノートのrank五の排除とrank四・三の七配置が得られる。
行列の半単純性は必要ない。rank四の三次元部分空間のGram rankは少なくとも二、
rank三の二次元部分空間のGram rankは少なくとも一である点も確認した。

七配置の全線形係数を比較した。`Q=xy+zw,h²+zw` は記載された二係数の対角場だけが残る。
残り五配置は、記載された `D1⁰+pN,D2⁰+qN` が全解を覆う。
特に次の非共鳴項を落としていない。

* `Q=hz+xy`：追加の∂w係数は `ph+rx+sz` と `qh+ty+sz`。
  `w'=w−rx−ty−sz` で `ph,qh` だけになる。
* `Q=xy+z²`：追加の∂w係数は `b0 h+bx x+by y+bz z` と
  `b0 h−bx x−by y+cz z`。
  `w'=w+2b0 h−2bx x+2by y` で `bz z,cz z` だけになる。
* `Q=x²+yz`：追加係数は `rx+py+sz` と `uh−rx+qy−sz`。
  `w'=w+uh−rx−sz/2` で `py,qy` だけになる。

`Q=xz+yw,hz+w²` では相似保存そのものが、表示されたN以外の追加係数を消す。
これらのshearはh,x,yとQを保ち、除外する係数値もない。

同ノートの五組 `(Qε,Nε)` について、全εで
`Nε(Qε)=0`、vertical divergence零、`[Di⁰,Nε]=0` を手計算で確認した。
`Di⁰` は `Qε∂z∧∂w` も保存する。
従って全p,qを保持する表示された `Pε` は真正の二次Poisson族であり、
ε≠0で二つの対角配置に移る。
対角配置の重み極限ではそれぞれ `bxy∂z∧∂w`、`bh²∂z∧∂w` だけが消え、
退化五次式は非零定数倍の `hxyzw` となる。

連結PGLの軌道極限は元を含む各成分に留まるので、SNC点に関する既存MPS定理が適用できる。
五つの実変形の特殊化には、既知40成分の**有限閉和**を用いる。
この順序により非零対数括弧の全支持の被覆を受理する。
特定の一つのMPS成分への包含や局所唯一性は、この議論には必要ない。

## 2. 零対数括弧の四全源：独立算術

再現入口は次の三ファイルである。

* `code/classification_three_hyperplanes_quadric_independent_review.py`
* `research/classification/jobs/three_hyperplanes_quadric_independent_review.json`
* `reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_REVIEW_20261107.json`

有限job gateを通した後、2.58秒で全項目がPASSした。
checkerは標準ライブラリのFractionによる疎多項式演算とAST literal parserを使い、
producer、SymPy、producerのJacobi・階数補助関数をimportしない。
算術の小関数だけを、既に独立なdependent-source reviewerから再使用した。

四配置それぞれについて、次を独立に再構成した。

1. `hxyQ` の三偏微分に対する全18列の二次syzygy行列。
   非零16次minorと二つの独立な厳密核から、QQ上の階数16・核次元2を確認した。
2. 二つのsyzygyのEuler成分を引いたanchor U,Vと、その基底規格化。
   `det(U,V)=xyQ` および `[U,V]=cU+dV` を厳密確認した。
3. 六つの磁気単項式 `h²,hx,hy,x²,xy,y²` を全て変数として保持した源。
   全10 strict Jacobiが恒等的に零、Pfaffianが
   `−s²xyQ`（sは明記された非零整数倍率）であることを確認した。
4. 全150 bivector係数に関する350本の射影Jacobi微分D。
   Euler方向25列、GL5方向25列、scalar一列、実際の六磁気微分を全て組み立てたF。
   producerのD・Euler・orbit/scale・全実族の四matrix hashが全て一致した。
5. 保存されたliteral normalと全実族にDを掛けると厳密に零になること。
   modular非零minorはQQ上の下界にだけ用い、厳密核でDの上界を与えた。
   Fの上界はFraction消去による全独立核関係を別途保存した。

結果は次の通りである。

| 配置(m,t) | rank_QQ D | rank_QQ F | 全核次元 | 真の余剰接方向 | 射影実族下界 | 全射影接空間上界 |
|---|---:|---:|---:|---:|---:|---:|
| (0,3) | 102 | 48 | 48 | 0 | 23 | 23 |
| (1,2) | 101 | 48 | 49 | 1 | 23 | 24 |
| (2,1) | 100 | 48 | 50 | 2 | 23 | 25 |
| (3,0) | 99 | 48 | 51 | 3 | 23 | 26 |

射影値への換算ではEulerの実効24次元と非零bivectorのscalar一方向を引く。
全六係数源とPGLの積は既約なので、(0,3)では実像の下界と全接空間上界が一致し、
その像閉包は23次元の全成分となる。有限代表点は全Poisson空間の滑らかな点である。
他の三行の成分性は、この階数だけからは従わない。

## 3. 射程

本reviewは、非零対数括弧の全源・実変形の手計算と、四つの零対数括弧源の普遍式・有限階数を確認する。
四配置の完全性を与えるfree curveの幾何、既知成分との相違、
残る一・二・三法線の二次障害projectionは、それぞれ別の証明依存である。
normalが全核を補うことはここで確認したが、その障害の算術を確認したとは主張しない。
全Poisson空間の尽くしは、このreviewの結論ではない。

## 4. (1,2)の残余方向：固定plane-cone経路の一次診断

四源の有限代表点のうち(1,2)を固定する。Fを全実族接空間（rank48）、Rを保存された真の一法線とする。
baseを(h,x,y)、fiberを(z,w)とする。
各Casimir座標H=h,x,yについて、H行を零、残るbase二変数とfiberの四mixed括弧を
任意のbase二次式、zw括弧をbase二次式＋fiber一次×base一次とする空間をS_Hとする。
dim S_H=36である。さらに三H全てを失う変形も許す全空間

    S_all = S_h+S_x+S_y

は、六mixed括弧の各base二次式36係数とzw括弧12係数からなり、dim48である。
各単独S_Hへの非包含だけでは、その線形和への非包含は従わないため、和も別途確認した。

| S | dim S | dim(F+S) | dim(F+S+CR) | rank(DS) | dim ker(DS) | dim(F∩S) |
|---|---:|---:|---:|---:|---:|---:|
| S_h | 36 | 68 | 69 | 20 | 16 | 16 |
| S_x | 36 | 68 | 69 | 20 | 16 | 16 |
| S_y | 36 | 68 | 69 | 20 | 16 | 16 |
| S_all | 48 | 72 | 73 | 24 | 24 | 24 |

全てQQ上の厳密計算であり、各行に対して
`λF=λS=0, λR=1` となる有理係数の疎linear functionalを保存した。
単なるmodular階数の観測ではない。
全150係数の射影JacobiはEuler-vertical追加に不変であり、
Fも全Euler方向を含むので、x/y-Casimir chartのSを元の全150座標へ直接埋めて比較できる。
Pだけ又はRだけを変更した不整合な微分は使用していない。

証明記録は順に以下である。各reportは実行時runnerの全bytesと依存hashを保持する。
先行reportのbytesを上書きせず、同じ有限runnerを順次拡張した。

* `CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_PLANE_DIAGNOSTIC_20261107.json`
* `CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_REMAINING_PLANE_CHARTS_20261107.json`
* `CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_INDEPENDENT_FULL_CONE_DIAGNOSTIC_20261107.json`

従って、実族・PGL・scalar・Euler方向を調整しても、Rの非零一次classは
この固定plane-cone ansatzで表せない。
この診断だけから、一次法線が零の高接触arcの不存在や、残余成分の不存在を推論しない。
別の全ambient変形又は支持障害の経路が必要である。

## 5. 残三族の実MPS変形：全係数の独立再構成

その後に得られた `M1_MPS` と `M23_MPS` の実変形は、上の接空間診断とは独立の証明である。
`code/classification_three_hyperplanes_quadric_mps_review.py` で、三つの全族を標準Fractionだけから再構成した。
有限job `jobs/three_hyperplanes_quadric_mps_review.json` のgate後、約0.034秒でPASSした。
証明書は `reports/CLASSIFICATION_THREE_HYPERPLANES_QUADRIC_MPS_REVIEW_20261107.json`。

t,α,β,γを全て変数のまま保ち、次を確認した。

* 全ten strict Jacobiが恒等零、係数は座標二次・tに多項式であり、保存された実Pfaffianと一致する。
* Y=T(t)xとする座標変更はdetT=t²。
  逆行列を使わず `T P Tᵗ=tD1(Y)∧D2(Y)+bQ(Y)∂Z∧∂W` を全係数で確認した。
  既知の非零対数括弧源そのものであり、t(1+tγ)≠0でQの階数は順に4,3,3となる。
* t=0の全源・三磁気項・anchor bracketが一致する。
  その二anchorと第2節で独立確認した元の全source anchorとの定数GL2変換を実際に解き、
  全係数一致と行列式非零を確認した。
* 六つのfiber translation coboundaryを微分から再構成した。
  その像は次元3で、三磁気項と合わせた六次行列式は厳密に8,2,−1である。

従って任意の六磁気係数が、これらの極限と射影的fiber translationで得られる。
既知非零対数括弧の全支持定理とMPS40の有限閉和の閉性を使えば、
m=1,2,3の全sourceがMPSに含まれる。
以前の法線の二次・三次診断から実変形を推論してはいない。
第2節のm=0新23次元成分と合わせて、独立三超平面＋既約二次曲面の全支持が被覆された。
