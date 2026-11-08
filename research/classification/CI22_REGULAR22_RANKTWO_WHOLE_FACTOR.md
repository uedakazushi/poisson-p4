# regular22_1・rank二 CI22 の全source：実二重quadric因子とE被覆

2026-11-08、PRIMARY E-FINITE-SURFACE。
**全rank四点を既知 H+2rank3Q 支持へ戻し、実既約被約quintic枝を空とする。**
新Jacobi operator、source profile、点、minor、GB又はjetは作らない。
既存の全十源・元350／126配列と独立受理済み逆場合分けを用いる。

## 1. 全sourceを失わない二つのmapと零補集合

surfaceは

\[
 S=V(2x_0x_1+2x_2x_3+x_4^2,\ q),\qquad
 q=x_1^2+x_3^2+x_4^2.
\]

`CLASSIFICATION_NONNORMAL_CI22_TANGENT_20261107.json` の model `regular22_1` にある
全百五十tangent／intrinsic conormal sourceを Euler商した十係数を p0,…,p9 とする。
全支持の逆場合分けは
`NONNORMAL_CI22_RANK_TWO_SEPARATE_REVIEW.md` §4で独立受理済み。
二乗式 (p2−p8)²=0 から、複素閉点の被約支持で p2=p8=t と置ける。

| case | 全parameter map | 残す自由係数 |
|---|---|---|
| I | p2=p8=t≠0、p4=−t/2、p5=m、他は零 | t,m |
| II | p1=s≠0、p4=2s、p0=r、p5=m、他は零 | r,s,m |
| III | t=p1=0 | 下記全五源A又はB |

case IIのzero tをgeneric kernelへ代入して作ったmapではない。
元式からp3=p6=p7=p9=0等を戻したwhole sourceである。
case IIIは正確に (p0,p4)*(p7,p9)=0 なので、
Aの自由座標(p0,p3,p4,p5,p6)、Bの自由座標(p3,p5,p6,p7,p9)の二線形sourceが覆う。
両全五源では全元350 Jacobiと全126 actualDが恒等零であることは既受理済み。
従ってrank四点はI又はIIに限られる。
scheme上で二乗式を一次式へ書き換えたとは主張しない。

## 2. 二つのwhole mapのactual Pfaffian

元の全126 quadratic rowsへmapを代入し、同時に元百五十係数から作った
signed Euler Pfaffianと比較した。全係数で厳密に

\[
             F_{\rm I}=0,\qquad
             F_{\rm II}=s^2x_4q^2.                    \tag{1}
\]

を得た。両mapの元350 Jacobiも全て零である。
追加compatibilityを自由と仮定したnecessary supersourceからの推測ではない。
r,m,s,tは全複素parameterを保持し、r又はmのgeneric値を選んでいない。

case IIの式は特に短く手で確認出来る。非零entryは

\[
\begin{aligned}
 P_{01}&=-mx_3^2,\\
 P_{02}&=rq-sx_1^2+2sx_2x_3-sx_3^2-mx_1x_2,\\
 P_{03}&=sq+mx_1x_3,\\
 P_{12}&=s(-x_1^2+x_3^2+x_4^2)-2mx_1x_3,\\
 P_{23}&=2sx_1x_3+m(x_3^2-x_1^2),\\
 P_{24}&=x_4(2sx_1+mx_3).
\end{aligned}
\]

他の上三角四entryは零。この疎性から実五次は

\[
 F=P_{24}(x_1P_{03}-x_3P_{01})
                  +x_4(P_{01}P_{23}+P_{03}P_{12}),
\]

となり、代入すればrは最初から現れず、mの全項が相殺してs²x4q²になる。
case Iの全十括弧と全係数零のreadbackも同名reportに保存した。
単にrankやnormal symbolを比較したのではない。

## 3. 全点の結論

IとIIIは一般rank二以下であり、指定rank四枝では空。
IIではs≠0なのでactualDは正確に

\[
                     D=H_{x_4}+2Q_q.
\]

qはrank三の既約quadric、Hはそのvertex線{x1=x3=x4=0}を含む。
base conic上ではx4=0は相異なる二点を切るのでsecant配置である。
既受理の全 H+2rank3Q 被覆はこの配置・全磁気・全cofactor／zero補集合を含む。
従って全regular22_1 rank四Poisson点は、名前付き既知成分の有限閉和に入る。
全PGL normal formを戻しても、実因子型と既知閉和所属は保たれる。

特に、Sを非正規支持に持つ**実被約既約quinticのrank二regular22_1枝は全空**。
normal CI22、nil32、singularM1、coneの別model又はprojected Veroneseへ、この結果を流用しない。
新しいcomponentの増分は零であり、全Eや全分類の完了も主張しない。

## 4. 検証と保存

`jobs/ci22_regular22_ranktwo_whole_factor.json` のgate後、
保存された元350／126配列だけを二mapへ代入した。
全百五十sourceを同じmapから再構成したPfaffianも(1)と一致した。
全計算 **0.108408秒、RSS51,312KiB**。
inline source、全二map、全括弧、全残差、factorizationと入力SHAは
`CLASSIFICATION_CI22_REGULAR22_RANKTWO_WHOLE_FACTOR_20261108.json` に固定した。
旧producer・reviewer・report・currentjob・主稿・global台帳・Gitは変更していない。
別作者の全逆coverを入力とし、今回の小さい代入がcoverの独立証明だとはしない。
