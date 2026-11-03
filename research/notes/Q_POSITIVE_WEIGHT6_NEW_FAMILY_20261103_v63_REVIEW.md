# 実三変数Poisson族・新しい25次元成分・下界58の独立監査

2026-11-03。Q_POSITIVE_WEIGHT6_NEW_FAMILY packetを独立に検証する。
補助modular N²条件を課さない三変数族は**全次数の実Poisson族**であり、
そのPGL5飽和像閉包は、全projective Poisson空間の新しい25次元既約成分となる。
従前の57成分とは異なるため、成分数の下界は**58**に上がる。
網羅性、一次ジェット像の優越性、一般fiberの五点尽くし、一般被約性は結論しない。

## 独立な全多項式検証

reviewerは新producerの構成関数を呼ばず、旧209式証明書の線形99×3mapと
三つの明示二次補正から全99入力を再構成する。
そこから全affine90係数を作り、五変数へ二次同次化する。
divergence25係数を独立に取り、P−Euler∧div(P)/5でcanonical150係数を再構成する。
全150係数が保存された新族と一致し、各係数のパラメータ次数は2以下である。

別の四変数微分・積engineで全140 affine JacobiをQQ[n0,n1,n2]上で確認する。
さらにreviewerの定数物理二次積表とC++/FLINT engineで、
全350個の**strict cubic canonical Jacobi**をパラメータ次数0–4の全係数で検証する。
族の次数が2以下なので、これは切断jetではなく全多項式恒等式の検査である。
新producerの350式はEulerを使うprojective quartic係数だが、reviewerは独立にstrict350も検査する。
全25 divergenceも零である。
唯一の旧補助cubic(4/3)n0²n2はN²条件であり、実族には課さない。

## rankから実既約成分まで

原点および(1,1,1)で、全GL25微分・scalar・三つの族微分を独立に作る。
保存150×29行列と全係数で一致し、QQ rank26と指定26次小行列式を確認する。
したがってGL5×Gm×A3の実像閉包Yは既約coneで、affine次元≥26。
その射影像は次元≥25である。

原点像Pは正軸eta1であり、既に独立検証した全canonical形式模型の実projective局所次元は25以下。
Yを含む全Poisson既約成分CはPを通る。
閉点における既約成分の局所次元はその大域次元なので、dim C≤25。
従って射影YもCも正確に次元25で、一致する。
単に接空間が大きいことや部分切断の次元から成分を主張していない。

## 従前57成分との差異

既存の独立監査済み51個のrank4成分カタログでは、射影次元≥25のものは
T01,T39,T40(各25)、CA,CB(各28)だけである。
次元上界はその保存カタログとlegacy定理の監査を入力とし、51族の全計算を新たに再実行した意味ではない。
正軸原点Pは、独立検証済み全閉包排除定理によりこの五つ全てに属さない。
閉成分がYを含めばPも含むため矛盾する。
その他のrank4成分は次元≤24なので、25次元Yを含めない。
rank≤2は閉条件で、族はrank4の点Pを含むため、残る既知六rank2成分にも含まれない。
よってこの全P4成分は従前57個のいずれとも異なる。

## 一般退化因子と、被約性の未解決

reviewerはcanonical150からEuler contractionと五つの4×4Pfaffianを直接作り、
全quinticがh³QとなることをQQ三変数上で確認する。
Qの5×5detは(1,1,1)で非零。
h=0制限の4×4detは**QQ[n0,n1,n2]上で恒等零**、同点でrank3なので、
一般にQは滑らかでHはその接超平面である。
一点のrank3だけから一般上界3を推論していない。

全strict350をcanonical126で微分した実ambient Jacobianも確認する。
QQ点(1,1,1),(1,2,3),(2,3,5)でrankは全て98、affine接空間28・projective接空間27である。
これは群＋族のrank26とは別の計算である。
この三点の接空間超過から、一般非被約性・nilpotent構造・他成分との交差を結論しない。
rank100による滑らかな開の証明は得ていない。

review code `code/q_positive_weight6_new_family_review_20261103_v63.py`、
receipt `reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63_REVIEW/`。
全入力と独立物理engine・出力のhashを保存する。
