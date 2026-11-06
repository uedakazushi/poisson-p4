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
