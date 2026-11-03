# 正軸における新成分の被約滑らかさ：多重度論法の独立監査

2026-11-03。新25次元成分 `R4_POSITIVE_WEIGHT6_25` を C とする。
既に独立監査された実五次 initial と、C に沿う全 Poisson scheme の一般長さ5を
組み合わせると、**C の被約多様体は全ての非零正weight6軸点で滑らか**である。
その点を通る25次元成分は C のみである。ただし24次元以下の別成分、埋め込み素、
全 Poisson 環の形式的厚みを尽くしたとは主張しない。全 Poisson 環は同じ点で非CMのままである。

これは新しい物理係数計算ではなく、凍結済みの独立計算を用いる数学的帰結の監査である。
作者の凍結packet `R4_NEW25_POSITIVE_AXIS_SMOOTH_REDUCED_COMPONENT_20261103_v63` の
証明・コード・receiptを入力として固定し、その主張と限定を独立に承認する。
付属コードは依存の全保存hash・既承認の数値と範囲・有限個の多重度不等式を検査する。
環論の論証をコードが自動証明したという意味ではない。

## 1. 同じ実全局所環での入力

正軸 eta=1 の実canonical tensor を P+ とする。全126方向・全350式から作った完成環には、
実際の GL と独立scalar の23方向を用いる形式積

    O_hat_(Poisson affine cone,P+) = C[[g1,...,g23]] completed-tensor A

がある。A は八変数 z0,...,z7 の商として得られる**全次数**形式横断環である。
これは補助的な固定pencil sourceでも、五次で打ち切った環でもない。
全射影化では自由scalar一本を除いて22個の滑らか方向を戻す。
`Q_POSITIVE_WEIGHT6_SMALL_MODEL_..._REVIEW` がこの実orbit23と95式の陰関数消去を検査している。

`Q_POSITIVE_WEIGHT6_INITIAL5` と独立reviewは、真の initial ideal に含まれる斉次ideal J5について

    dim C[z8]/J5 = 3,  degree C[z8]/J5 = 5,
    radical J5 = (z1,z2,z4,z6,z7) intersect (z1,z3,z4,z5,z6,z7)

を厳密QQで証明している。C の実次元26と上の23方向の積から dim A=3 である。
従って J5 と真の initial ideal は同次元で、Hilbert多項式の上界により e(A) <= 5 を得る。
五次idealと全initial idealの等号を仮定していない。

独立に承認された `Q_NEW_FAMILY_LENGTH5_COMPLETION` とそのroot reviewは、C に沿う
全 Poisson scheme の一般長さが正確に5であることを証明する。
その根拠は実閉fiber長5、一般Jacobian98、一般augmented99、有限自由rank5の論法であり、
閉点での非被約厚みを一般へ外挿したものではない。

## 2. 完備化で枝が分かれても多重度比較は有効

C は連結群 GL5 とscalarに不変な被約既約成分である。上の実群作用による形式積は
C のidealにも適合する。したがってその被約完成局所環は

    O_hat_(C,P+) = C[[g23]] completed-tensor B

と書ける。B は A の商である。有限型複素多様体の優秀性からこの完成環は被約かつ
等次元26であり、B は被約かつ等次元3である。**代数的Cが既約であることだけから
完備化後も整域だとは仮定しない。**

B の極小素を q_i とする。全 Poisson scheme の一般長さ5は各 q_i でも保存される。
実際、優秀環の完備化はregularかつflatであり、元のCの一般点上にある各極小素での
closed fiberは体である。長さ5の組成列をflat base changeすると各単純商は一つの
剰余体になり、局所長さは5のままである。23形式変数の積を除く操作もこの長さを変えない。

Hilbert--Samuel多重度の結合公式を全Aに適用し、Cから来る等次元極小素の寄与だけを残すと

    e(A) >= sum_i length(A_(q_i)) e(A/q_i)
         = 5 sum_i e(A/q_i)
         = 5 e(B) >= 5.

B が被約なので最後から二つ目の等号では各一般長さが1である。
一方 e(A)<=5 なので

    e(A)=5,  e(B)=1.

また他の三次元極小素は正の寄与を持つため存在できない。
多重度公式は低次元極小素と埋め込み素を数えないので、それらは排除されない。

## 3. 被約等次元環の多重度1から正則性

ここで用いる標準事実を、解析的整域性を先取りせず説明する。
完備な被約等次元三次元C局所環 B で e(B)=1 とする。無限剰余体上でmの最小reduction
(t1,t2,t3)を選ぶと R0=C[[t1,t2,t3]] は B に有限に入り、e_(t)(B)=e_m(B)=1 である。
各極小商はR0上有限かつ同次元3なのでR0の非零元は各極小素に入らない。
被約性からBはR0-torsion-freeとなる。R0上の一般rankは多重度と等しく1である。
従って B は Frac(R0) 内の有限整拡大として埋まり、R0の整閉性により B=R0 となる。
よって B は正則であり、解析的枝も一つである。

これを適用すると C は P+ で滑らかである。射影Cの次元は25である。
非零正軸は独立監査された有限対角変換で一つのprojective軌道をなすので、
同じ結論は eta!=0 の全点に移る。eta=0 の Q_(8/9) での滑らかさは結論しない。

## 4. 接平面と正しい限定

B の接錐は三次元線形空間である。Aへの商写像からそれは V(J5)=L3 union P2 に含まれる。
既約三次元の線形空間は二次元P2に含まれず、唯一の三次元成分L3と等しい。従って

    T_(C transverse,P+) = L3 = {z1=z2=z4=z6=z7=0},
    free coordinates = (z0,z3,z5).

全 Poisson 完成局所環の埋め込み次元8と次元3は変わらない。
CMなら e(A)>=edim(A)-dim(A)+1=6 となるため、e(A)=5 は全Aの非CM性とも整合する。
滑らかな被約成分と非CMな全環を混同しない。

このpacketは、正軸における新成分の被約滑らかさ、唯一の最大次元枝、実多重度5を証明する。
全局所環が新成分だけに支えられること、残るP2の実枝への持ち上げの成否、全initial生成、
全成分にわたる正規性、他のQaへの所属は未解決である。

review code: `code/q_positive_weight6_reduced_component_smooth_review_20261103_v63.py`。
receipt: `reports/Q_POSITIVE_WEIGHT6_REDUCED_COMPONENT_SMOOTH_20261103_v63_REVIEW.json`。
