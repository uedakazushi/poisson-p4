# Raw compact Yの純3次元性の反証と孤立無限遠のsheaf分離

2026-11-03。新normal-direction graphのraw compact Yに対して、「純3次元あるいはpure3のCohen–Macaulayを証明すれば、既存のaffine次元≤2と合わせて未知解を消せる」という方針を調べる。この**変更しないraw Yの純3次元性は偽**である。既知の孤立無限遠を実際のsheaf idealとして分離する正確なcolon recipeも指定する。colonの展開基底や、残余全部のw飽和は未計算である。

## 1. 元134式だけに属する孤立点

`S61_MINIMAL_AXIS_INFINITY_20261101` は、非零最小冪零modular型の16次元既約incidence成分を構成し、W16を優越しgenerically finiteであることを証明した。一般fixed betaでこの成分の点は無限遠にある被約孤立点である。sampleのbeta微分16minorは569で、full-homogenized45minorは

    -124102513209883356037120000000000000000/59049

である。重要なのは選んだ45行が**追加G6を含まない元134式だけ**に属することである。今回のproducerは元134の係数から最高q部分、44projective q列、homogenizer列を直接再構築し、このQQ determinantを新たに計算した。beta minorも保存行列から再計算した。旧source全hashを確認する。

この点のrelative tangentは零で、beta方向が全射なのでtotal incidenceはそこでregular dimension16、source closureが唯一の成分である。full45minorが非零な開はこの成分の非空開であり、componentと悪い閉集合のproper projectionを使ってbetaを縮めれば、その全generic fiberは被約孤立点となる。

normal=0の射影無限遠はtinyplane degree3 spanから空なので、この無限遠成分はnormal-direction graph Yにも現れる。一方、全135のraw infinityには既存の3次元binary/flag族がある。それは134のYにも含まれる。したがってYには一般betaでも0次元の成分と3次元以上の成分が共存し、純3次元ではない。単に最大次元やinitial Hilbert値を見る問題ではない。

## 2. 保存minorからglobal Deltaを定義する

元134式を全てq次数2へwで斉次化する。全46座標に関するJacobianから、保存45行と、sampleで固定したq座標一列を除く残り45列を選ぶ。これによりglobal homogeneous section Delta（次数45）が有限determinant recipeで指定される。そのsample値は上の非零数である。

Y上でDelta≠0ならprojective relative Jacobian rank45なので、その点は被約孤立点である。体K=C(W16)上ではD(Delta)∩Yは有限な被約0次元成分の和で、clopenである。そこからw=0だけを選んだ部分Y_infもfinite etaleかつclopenである。これは単に「supportが消える」という主張ではない。

## 3. 全未知アフィンを保持するexact sheaf ideals

Yのideal sheafをIとする。射影のirrelevant saturationはsheaf化に含め、

    I_inf = (I+(w)) : Delta^infinity,
    I_rest = I : I_inf

と置く。第一式はD(Delta)内の単純孤立点のうちw=0にあるものだけを閉じたschemeとして取る。上でこの部分が既にclopenと分かったので、余分なembedded closureを加えない。第二式はそのidempotent idealのannihilatorであり、**single colon**で補成分を与える。

従ってsheaf idealsとして

    I = I_inf intersect I_rest,
    I_inf + I_rest = O_ambient,
    Y = Y_inf disjointunion Y_rest

である。D(w)ではI+(w)=OとなりI_inf=Oだから、I_rest=Iとなる。従って未知アフィンの全点、正次元成分、非被約構造は正確に保持される。Deltaで単純な未知アフィン点を直接saturateして消す方法ではない。

このrecipeは少なくとも保存axis成分の非空generic部分を取り、generic betaを縮めれば対応するfinite etale成分を丸ごと分離する。選んだDeltaが零になる他の孤立成分まで全て取ったとは主張しない。全45minorsを用いるFitting版に替えれば全単純孤立無限遠を同じ方法で扱えるが、ここでは巨大minor集合を展開しない。

## 4. 範囲と次段階

これは実際のscheme分解の有限ideal-operation recipeであり、colonのGroebner basisやfinite etale代数の次数は計算していない。Y_restがpure3であることも未証明で、他のw-torsion成分や交点でのembedded構造が残る。raw Yを変更しないCM探索や、その全体の空性検査は行わない。

一般アフィン五点尽くしに向けては、正確なw飽和やdirection-chartのB_i nilpotenceに関するQQ identityがなお必要である。modp nilpotenceだけを一般QQへ移すことはしない。

保存: `code/s61_affine_raw_purity_split_20261103_v63.py`、`reports/S61_AFFINE_RAW_PURITY_SPLIT_20261103_v63/`。旧20261101証明書は変更していない。
