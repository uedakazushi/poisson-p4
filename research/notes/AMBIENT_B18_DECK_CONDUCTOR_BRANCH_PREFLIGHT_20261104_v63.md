# 二つの deck/conductor 実 scheme 写像と、長さ12に不足する source gate

2026-11-04。source/theorem の調査。新 jet、GB、beta 点、高次 core 係数は計算しない。

## 1. 実 root sign による二つの全次数写像

二つのモニック関係の unit Schur 後の正則九環 S を使う。
R0 の二消去座標の解は tauA²,tauB² のみで決まるため、
sigmaA:tauA→-tauA と sigmaB:tauB→-tauB は R0 を固定する実 S の involution。
Npair の像は L=6*(tauB-beta(X,tauA))=0 である。

N 正規化 B_N 上で実12出力差 Fj は零。
sigmaA へ移した源では Fj=-2*bAj*tauA、sigmaB へ移した源では Fj=2*bBj*tauB。
一方

    L|sigmaA(N)=tauA*UA,
    UA=6*(beta(X,tauA)-beta(X,-tauA))/tauA,
    UA(0)=-48,
    L|sigmaB(N)=-12*tauB.

UA は odd Hadamard の正則商で、原点の -48 は実 Npair tangent の beta_tauA=-4 から得る。
従って UA は単位。tauA,tauB は B_N 内で非零なので、まずこの整域内で消去すれば

    Gj|sigmaA(N)=-2*bAj/UA,
    Gj|sigmaB(N)=-bBj/6.

これは conductor を割った環で tau を逆にしているのではない。
整域内で証明してから conductor quotient へ制限するので、そこで tau=0 でも有効である。
各側の全12個 b は元8個の I0 ideal を生成するため、厳密に

    J|sigmaA(N) = IA*B_N,
    J|sigmaB(N) = IB*B_N.

従って二つの**実 scheme** Spec(B_N/IA B_N), Spec(B_N/IB B_N) から core への写像を得る。
既に定義した正則トレース根を restrict し、conductor 上で元の0/0表示を再評価しない。

## 2. 各 source は三変数三関数にできる

実 conductor pullback の12x8一次行列は両者 rank5。
rows18,19,60,28,29 と、Nsource の X108,X109,X113,X118,X119 の五 minor は
CA が16、CB が-1/64。
同じ Nsource28/20切断と純CA48/CB48を用いてこの五関数を unit IFT で消せば、
残る三 source パラメータは X112,X117,p_N、残る三関数は b67,b72,b63。
これらは具体的な元 DAG 関数の有限入力であり、未知の branch 方程式を代入してはいない。
三関数が二つで生成される、CIである、何本の曲線を持つ等はまだ主張しない。

## 3. 初期 divisor の incidence は source branch の証明ではない

固定正規化 Phi へ四つの root-sign hyperplane を代入した全 binary quartics を保存した。
row63 section の根の incidence は次の通り。

| binary 因子 | 指数 | その方向を含む conjugate hyperplane |
|---|---:|---|
| u-w | 3 | L0 |
| 9u-44w | 2 | L0,LA |
| 19u-124w | 2 | 無し |
| 43u-148w | 1 | LA |
| 611u²-4792uw+7856w² | 2 | LB |

最後の因子は二つの node preimages を表すが、射影四次曲線上では**一つの節点**。
これを二つの相異なる actual core branch と数えてはいけない。
R、Pa、Qb との対応は有力な候補であって、表の接方向だけで高次 lift や generic 長さは証明しない。
長さ3の main-Npair 曲線と長さ2の critical-CA 曲線を必要な source として挙げられるが、
この二つの実写像・長さは未構成であり、長さ12の下界はまだない。

## 4. 既存三 strata の長さと、正しい輸送の gate

既存 source-bound な全環から、R の normalizer-pair meeting は ordinary double、
Pa の余分な CA pair は smooth20、Qb の余分な CB pair は Z4=k[u,v]/(u²,v²)。
generic 長さはそれぞれ2,1,4である。
Qb の元 source metadata は二 Artin tangent 座標を physical f=row18 と rho=row108 に取る。
両者は現在の torus で weight0、physical w140 は weight1と直接確認した。

これらを現在の core の multiplicity 下界へ運ぶには、以下の source 図式を明記すべきである。
P0 近傍の proper normalization chart が既知の punctured Pa/Qb flag pair を含み、
valuative criterion でその arc が P0 の唯一の幾何 source 点へ延びる。
target が Gm-fixed、flag fiber が有限なので connected Gm は各 flag を固定する。
fixed-locus は fiber product と可換で、共通 G19 の equivariant unit section を除いた後に
現在の zero-source tensor に入る。Qb の有限 Artin algebra の tangent generators は weight0 なので、
その全 algebra が weight0。従って正 weight w を零にしても長さ4を失わない。
Pa の smooth余分枝は reduced curve になり、R の長さ2は既存 zero-slice 認識を用いる。

この図式と現在の local source chart の完全性が確認できれば、三つの異なる曲線が core へ入り、
associativity から e(core)>=2+1+4=7 を得る。Pa/Qb の反対 flag count と R の非N性から相異なる。
各曲線の原点 multiplicity は少なくとも1で十分であり、接線を推測する必要はない。
本 packet はこの輸送 application を別受理 gate として残し、literal な原点 arc 係数は構築しない。
まして不足する5を接錐の根から足して e>=12 と宣言しない。

## 5. CM 完了へ必要な追加証明

実三関数 CI curve C3 は CM・e12。full core はその四 graph ideal H による商。
同じ一次元で e(full core)=12 と source から証明できれば、H は有限 support の submodule。
CM一次元環には非零有限長 submodule がないので H=0となり、全 graph 冗長性と CM が従う。
現在不足するのは actual source branches と generic 長さの下界であり、二次/三次接方向の数ではない。
代わりに source14/17 polynomial model の determinantal 完全性を証明する方法も考えられるが、
それもまだ入力構築前の別提案である。
