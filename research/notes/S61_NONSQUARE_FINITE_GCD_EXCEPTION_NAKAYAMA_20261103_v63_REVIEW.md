# 有限gcd例外cover上のraw QQ全射性：独立監査

2026-11-03。author packet `S61_NONSQUARE_FINITE_GCD_EXCEPTION_NAKAYAMA_20261103_v63`の
**有限閉cover Sigmaだけ**を対象に、算術とfinite-model証明を独立に確認した。
全三parameter open U上のQQ rank24・元のproduct rank41やambient removalへ広げない。

## 入力と実raw sourceの独立再構成

凍結archive全memberのsizeとSHA256をstreamで照合し、新しいdirectoryへread-only inputのcopyとして展開した。
author sourceやwork fileは変更していない。
四つの実cleared quadricの804個のparameter係数から、P=a*dへの変換と各rowのa-shift(0,2,2,1)を独立に再構成した。
六つのlight linearを掛け、全24×56 cubic係数を生成すると、literal raw sourceの**全4824項**がQQ上で一致する。
元source SHAはf9590655c56c4c7790d43bc780950ff9e7bf4a3339a4a0879e2901369a3b78aa。
後続のSchur residualではない。
全QQ分母はp=32003の単元である。

さらに同じ四quadricから、W=0での実Q1=A*B/6075、Q2-combination=40*x0*B/6075と
3Sc−7Vc=20Pを確認した。
Q0のsplit表示4Q0=U*R−5V*T、Bのradical係数−20a*DB、DB消去後のdual normからN4を再計算した。
実Q3をB=W=0へ制限した四変数二次式から、全4×4係数行列H3も再構成した。
保存済みancestor係数と全て一致し、P^6*(20P+3e+20)^2*Jを除いた行列の指定3次minorが
P*(20P+3e+20)*既知QQ単元*M9となることを確認した。
これでNとMは単なる名前の多項式ではなく、実rank2 gcd例外の必要条件へ結び付く。

## 有限V-modelにescapeがない理由

V=Z_(32003)、a=−P*(20P+3e+20)/180として

    S=V[P,e]/(N,M)

を取る。保存されたQQ sparse cofactorを、独立のFraction多項式積で掛け合わせ、

    R(P)=A(P,e)N(P,e)+B(P,e)M(P,e)

を文字通り確認した。
N,M,R,A,Bの全係数はVへ入り、Nのe^4係数72442188とRのP^28先頭係数はVの単元である。
従ってSにはe次数4とP次数28のmonic関係があり、Sは有限V-moduleである。
例えばe^i P^j、i<4,j<28で生成できる。
**flatnessやfree rank28は仮定も結論もしていない。**

FULL localizerは

    U=a*P*(500a+90e²+3e−24)
       *((50a+18e²−9e−2)P+(60e−40)a)
       *(100a+18e²−9e−2)*(5P−6e−4).

特殊代数全体でのU単元証明は、特定のopenに制限した証明ではない。
有限V-module S/USのmod pが零ならNakayamaによりS/US=0、従ってUはS自身で単元。
この論証はSの全componentとnilを保持する。
特にP、a、J、および20P+3e+20=−180a/Pも単元なので、上のminor正規化は元coreの既知単元しか除いていない。
非proper openでのmodp unitからのQQ promotionを使っていない。

## 特殊代数の独立証明

独自のdictによるFp多項式演算とnormal formで、次を確認した。

- literal N,M,Uの全係数はQQ式をp-unit分母で還元した値と一致。
- 5個のmonic reduced GBについて、各generatorの(N,M)へのsource liftと、N,MのGBへの還元が成立。
- 全10 S-pairが零へ還元する。
- standard monomial全体は正確に28個。
- N*Abar+M*Bbar+U*Cbar=1がFp[P,e]内で文字通り成立する。
- 完全normal formからP,eの全28×28 multiplication matrixを再構成し、全1568係数が保存値と一致。

GBを新しく計算したり、saturation APIの結果を前提にしたりしていない。
Buchberger判定と両ideal包含を確認したため、28次元特殊代数のnilpotentも落としていない。

## raw moduleの独立native検証

authorのFLINT matrix/rank codeは使用せず、新しい標準C++整数剰余実装を用いた。
独立に再構成したP,e作用からaの作用を作り、全literal raw係数を評価した。
28次元regular representationで**全672×1568=1,053,696係数**を構築した。
P,eの可換性、N=M=0、Uと保存cofactorによる逆の左右両積も検査した。

その実行列から指定672列を取り、保存inverseとの**左右両積**を独自の64bit整数dot productで確認した。
検査した恒等行列係数は合計903,168個。
新しいrank solverやGauss計算を必要とせず、全行射影のright inverseが実sourceに存在することを確認した。
native部分は約0.2秒。

## NakayamaによるQQ結論と限定

raw cokernel Kは有限S-moduleなので有限V-moduleである。
上の正則表現上の全射性からK/pK=0。
NakayamaによりK=0、QQへtensorして、**閉じた有限cover Sigmaの全座標代数上でraw24×56が全射**となる。
従ってその全複素点を含み、特殊fiberで見える一部分だけの結果ではない。

この証明が扱うのは、実rank2 gcd例外を含むこの有限coverである。
全U、他の高rank proportional場合、pair gate・five-quadric gate、全QQ product41、quarticやambient removalは別課題。
元product41への翻訳には、別に証明された17方向kernelと実座標bridgeが必要である。
CRTや有限体の全open結果を無条件にQQへ移していない。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/s61_nonsquare_finite_sigma_review_20261103_v63.py --output <new-output-directory>

review codeと独立native source、operator/raw-modp/逆確認receiptを保存する。
archive inputは検査後にcopyを削除し、ELFも保存しない。
新しいGBなし、120秒・2.2GB、全監査は約5秒。
