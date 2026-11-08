# nil22長translationの最後のrank二core：有理全軌道入口と停止したT19 frame

2026-11-08、PRIMARY E。純手の全source／実逆／不成立frame診断を保存する。
新operator、point評価、rank、minor、profile、GB、jet又はCASは無い。
既知親又は新componenthoodを主張しない。旧証明書とGit/global台帳は不変。

## 1. 正確な残sourceと有理代表

独立受理 `TWO_QUADRICS_NIL22_LONG_FIVE_WHOLE_SOURCES.md` §5の
\(\eta=1/12,\xi^2=-1/432,\kappa_0=-1/6\) を扱う。
二符号は実normalizerで交換され、全s非零は既明示weighted scalar mapで同じactual軌道、
s零はこの同じ全Poisson源の閉境界である。
従ってこのwhole枝のnamed parentには全s非零の一軌道の陽性包含で十分。
陰性chart判定から閉parent非包含を推論しない。

\(j^2=-3,\xi=j/36,s=j\) を選ぶ。
これは対象のs単元を一つ選んだ全軌道代表であり、全sourceを選択sliceへ縮めてはいない。
actual map \(z'=jz,w'=jw\)、h,x,y不変、tensor scalar一を掛ける。
逆はz,wをjで割るだけ、determinant \(j^2=-3\ne0\)。
全sの明示正常化を先に戻せば、元全20sourceへの実逆も保持される。

新座標の \(R=hx-y^2/2,S=hz-yw\) を用いると全六括弧は
\[
\begin{aligned}
 P_{xy}&=R/2+w^2/144-hw/12,\\
 P_{xz}&=2xw/3-zw/18-yz/3-3y^2-w^2/12+yw,\\
 P_{xw}&=S/6+w^2/36-3hy,\\
 P_{yz}&=-S/2-w^2/12-3hy+hw,\\
 P_{yw}&=-3h^2,\\
 P_{zw}&=6R+w^2/12-3hw.
\end{aligned}                                                \tag{1}
\]
全係数はQQ、h行零、cofactorは元と同じ \(y\partial_x+h\partial_y+w\partial_z\)。
sourceの \(n_{17}\) はこのmapで−3となるので、(1)をn17一の同じ正規係数表と呼ばない。
この有理化そのものの新算術実行はしていない。

## 2. 因子／pencilの必要幾何

元 \(U=R-w^2/24,W=S+jw^2/6\) の二factorは
\[
 q_+=U+(j/6)W=R+(j/6)S-w^2/8,
 \qquad q_-=U-(j/6)W=R-(j/6)S+w^2/24.
\]
q−は \(h(x-jz/6)-(y-jw/6)^2/2\) でrank三、q＋はrank四。
全pencilはKronecker M1＋regular二次Jordan型で、rank三memberは一つのdouble root。
共通radicalは零である。

T20は全completed sourceで二quadricに共通radicalを持ち、proper incidenceの閉性が必要比較になる。
T17はpencilにrank二以下のmember CDを持つので、そのproper incidenceを本pencilが満たすかが必要比較になる。
これらの幾何条件だけを他成分全体の非所属へ外挿しない。
T18のmoving-H曲線が別coreを覆うことも、本coreへの帰属そのものにはならない。

## 3. 二つの非閉鎖T19 frameを停止する

T19のh-Casimir sourceは
\(q_4=AB+2tDH/(r-s),q_3=CD-uB^2/(r-s)\) を持つ。
真全completionを保持した必要factor frameとして
\[
 D=h,H=h+eW,B=h+eY,C=h+2eY+e^2c,
 A=-h+e(Y-W)+e^2a,
 \quad t=(r-s)/2,u=r-s
\]
を試みた純手展開では、q4／q3はorder e二で所望rank四／三を持つ。
ただしr=er0,s=−er0の必要leading balance後にも、actual inverse及び許容Euler gaugeの
\(P_{aY}\) に \(r_0hW/e\) が残る。従って**このframeは実有限Poisson曲線を与えない**。

H−Bをorder e二、Aの独立項をorder e三にする第二valuationでは、
rank四Pfaffianのbalanceをr=order e二,s=order eと出来るが、
canonical modular場のlimitは二本の長二block、すなわちN二乗零へ落ちる。
本対象のJordan(3,2)とは異なるため、このrouteも現在のsourceには閉じない。
二つの不成立を受け、同じframeで係数又はorderを反復して探すことを停止する。

これはT19の閉包への**非所属証明ではない**。
既知parentの別actual whole family、又は全軌道から導かれる一回の最弱fullambient判定が次の独立義務である。
新anchor/operatorは親の数学的契約／gate／GOより前には実行しない。
