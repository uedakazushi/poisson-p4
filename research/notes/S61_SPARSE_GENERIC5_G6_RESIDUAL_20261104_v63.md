# G6の実五変数帰結と第四法線二次関係への還元

2026-11-04。D42やD84の既知残余normを再計算せず、低次数の実帰結G6を今回の疎五変数源へ移す。全45未知数と法線二次項を保持する。残余の空性・generic五点尽くしは未証明である。

## 既存の分離と残る目標

既存polarizationは全アフィン代数をE5×Bへ分離し、D84は同じknown-five idempotent eをnormalだけで表す。D42、D84の非零性や原点での単元性だけではBは消えない。raw無限遠には孤立axisと3次元以上の族があり、raw purityや特殊化後のt飽和空性を一般空性へ移すこともできない。残る実際の有限目標は、全Bを保持したnormal idealのirrelevant saturationが単位になるというQQ係数恒等式である。

## 小さい実G6と125線形cofactor

保存された全61変数のG6は9257項、weight6、q次数1または2で、元の70本F4のq一次係数結合G6=ΣK_i(q)F4_iである。今回のhistorical β16の非零九係数を用いる。primitive gaugeに置き換えたりLee係数を消したりしない。

先行support計測では1698項だけが生き、G6のpull後support上界1947、cofactor右辺の相殺前support上界2851、合計4798項であった。674個のq³位置も右辺に一時的に現れる。今回の一回のQQ計算でG6を全係数構成し、125個のK_i係数を元134式のpullへ掛け、全q³相殺と全q次数≤2の恒等式を検証する。

G6を既存q45 forward回路へ代入し、その五汎関数を以前と同じP,C,W0,W1,dの有限回路で作る。既存のadj(B5)/Delta5でlambda_G6を求め、R_G6=G6−Σlambda_G6,j Fselected_jとする。既知五点で消えることと五汎関数の単射性からplane制限は零である。全42 linear-normal factorと126 cross成分を保存し、normal二次項も落とさない。

元の125-cofactor恒等式も、既存129残差と選択五式へ正確に運ぶ。cofactorは未知qに線形だが、その逆元は一つも取らない。G6は新しいアフィン条件ではなく、低いq次数をもつ既存ideal帰結である。

## FULL135の別の同変性gate

RAW134の同変性だけからG6の行空間を移してはならない。`S61_HOMOGENIZED_EQUIVARIANCE_20261101` と元のweight6 reviewの別の証明を使う。70最高二次式のq一次syzygy空間Sは36次元、Schouten部分Bは35次元でSp4安定。Sp4の一次元characterは自明だから、選んだkに対しgk−kはBへ入る。

同変な低Jacobi特解Q0(β)を使うと、全Bianchi恒等式[q,J4]=−[Q0,J4]−[L,J5]によりG6の変換差はβ²F4+βF5型になる。保存gaugeへのqのβ²平行移動も同型の差だけを生む。これは最高q項だけでなく、G6全多項式の変換についてのq非依存係数の恒等式である。従ってFULL135のQQ(β)-行空間が保存される。

元の厳密証明書と独立reviewは、以前の同じβでq最高二次係数のrank135をQQで証明している。最高q次数のrankも可逆affine q変換とparameter行空間同一視で保存される。Sp4×weightedGmの疎源dominanceにより、今回のK上でも135行は独立である。この議論は古い固定135minorが今回のg=identity源で非零だという主張ではない。

## 新しい有限代数目標

plane制限rank5を引くと、FULL135からの残差空間は130次元。RAW134だけのcross rankが既に126なので、FULL135でもcross写像のrankは126である。従ってそのkernelは**4次元**となり、K上で独立な四本のpure-normal多項式(次数≤2)が存在する。元134だけの場合の三本から一方向増える。ただし、今回その四本の具体的な基底を計算したわけではなく、G6残差自体をpure-normalだとも呼ばない。

固定126minorを選ばない完全な有限recipeもある。Xを実130×126 cross行列とすると、Kはformally realであり、Cauchy–Binetによってdet(XᵗX)は非零minorの平方和なので非零。L=(XᵗX)⁻¹Xᵗ、P=I130−XLと定められる。Pはrank4、PX=0で、pure-normal部分U(n)への積P U(n)はこの四次元関係空間を与える。これは存在と有限演算recipeであり、126×126 Gram逆行列やその係数を計算済みとする主張ではない。

同じLでz_h n_kをdegree≤2のnormal多項式P_hk(n)へ戻せる。四次元のpure-normal関係、n_iP_hj−n_jP_hiのcross compatibility、選択五式のn_i²で分母を払ったquarticから作るnormal ideal Jは、D(n_i)の全chartで実残余Bと一致する。従ってBの空性はJ:(n0,...,n41)^∞=(1)、同値にあるdで(n)^d⊂Jを示す有限所属証明へ帰着する。G6は同じidealに属するが、この低次数表示の第四関係を直接利用できる。

この所属証明、Gram係数の展開、normal Macaulay消去、G28をannihilatorとする主張は今回行わない。特にG28が既知proper firstjet像を定めることだけから、全未知BがG28=0へ写るとは推論できない。

再現：`code/s61_sparse_generic5_G6_residual_20261104_v63.py --output work/NEW-output`。一回60秒・2GB・32MiB、全新明示parameter supportは合計150000以下。新しい点・素数・GB・degree探索は使わない。

## 今回の確定結果

一回の実行は1.3769秒で完了した。実G6のpullは1940項で、q次数は1,2のみ。125-cofactor右辺の674個の可能なq³位置はすべて厳密に相殺し、全q次数≤2の係数もG6と一致した。右辺replayの最大supportは2008、保存scalar atomは合計2150項、保守的な新explicit support合計は6098である。

新しいfactored回路は15289 node、RAW411120 byte、gzip約85KiB。G6の直接出力、五汎関数、lambda_G6、R_G6、42法線cofactor、126cross成分と、旧129残差へ運んだcofactor恒等式を保存した。全normal二次項と同じ45未知数を保持する。新しい126逆行列、四本のpure-normal基底、残余unitは計算していない。

## 五点planeの一様な小さい最高次数模型

トーリック三点のtrace零と非共線性により、幾何的係数拡大上で五点は
(1,0,0),(0,1,0),(-1,-1,0),(0,0,-1/2),(0,0,1)
へ移る。この点で消える五次元の二次以下空間には、既存五汎関数proofの基底

```
1+sigma-2sigma²-a²+ab-b²,
a-a²+2ab, b+2ab-b², a*sigma, b*sigma
```

を使える。最高二次式の商では

```
a²=2ab, b²=2ab, sigma²=−3ab/2, a*sigma=b*sigma=0.
```

a²b=2ab²=4a²b、ab²=2a²bより3a²b=3ab²=0、従ってa³=b³=0である。またsigma倍の第三式からsigma³=0、他のmixed cubicも零である。最高二次式五本は独立で、三次単項式十本はすべてその一次倍で生成される。従って最高次数商のHilbert関数は(1,3,1)、長さ5である。

元の五式の商は5個の単項式で生成され、五つの相異なる点へ全射するので、実際にこの被約五点のidealを生成する。rank5/rank10と有限生成は係数拡大から降下する。これは特定の古いquinticや固定minorの輸送ではなく、今回の既知五点のアフィン配置そのものからの証明である。

この小さいplane模型とcross rank126を組み合わせれば、今回の実AはK[n42]上で高々5生成となる。cross rewriteはnを含む項のz次数を減らし、最高三次の全射で残る純z三次を次数≤2へ落とし、五最高二次式で二次単項式を一つに減らせる。実残余Bのnormal閉埋込みと上記の低次数Jによる全chart表示に使える。**K上の全Aが有限長だとは結論しない**。

## 次の実証明の一単位

次段で追う有限targetは、上記Jに対する具体的な(n)^d所属certificateである。G6の第四pure-normal関係を使い、すべての法線を残したまま、多項式cofactorの恒等式として検証する必要がある。normal原点でのJacobian可逆性、特殊化後の空性、G28の存在、Gramの可逆性はいずれもその代わりにならない。

これにより「どの残余をどのidealで消すべきか」と「G6が追加する低次数方向」が実係数と有限operatorで固定された。大きいnormal Macaulay行列やGram展開を自動で開始することはせず、次に使うcertificate空間のsupportと次数を別予算で測る。
