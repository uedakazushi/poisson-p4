# 新25次元成分の一次像15・一般像fiber1の独立監査

2026-11-03。全7変数のraw jet行列と関数体階数を独立に再構成し、
新25次元成分のprimitive beta像が正確に15次元であることを検証する。
W16を優越せず、一般betaの五点問題への反例とはならない。

## 全sourceとgeneric上界

reviewerは保存実affine三変数族から、一般評価点p4を残した値Cとlinear jet Lを取り直す。
三つの族方向・四つの評価点方向・GL4全16方向・perspective全4方向・scalarをすべて保つ。
perspective式はproducerの閉公式をコピーせず、二次vector field
V(t)=−t_i*t のpushforward微分からconstant/linear係数を抽出する。
全30×28多項式行列が保存された行列と一致する。

非零有理定数による20回のSchur消去を、全7変数の多項式環で各行・列indexとともに再計算する。
これらは全係数環で可逆なので、generic rankは20+rank(残差10×8)である。
残差について**別backendのSymPy QQ(n0,n1,n2,p0,p1,p2,p3)**で正確なrrefを計算し、
rank5・残り五行全entry零を検査する。
さらに独立に書き出した入力で、凍結済みFLINT有理関数体engineも再実行しrank5を確認する。
従ってraw image次元25は全係数体での上界・下界であり、標本rank25の外挿ではない。
別の25次有理minorは、その下界の具体的証人として確認する。

## nuisance10方向とprimitive像

非退化なCのGL4値方向はrank6。
perspective方向は値を変えず、linear jetに独立なrank4方向を与える。
同じ行列で両rankと値Pfaffian−16/81を確認する。
GL4とperspectiveを全てsourceに含めているため、局所的な値正規化とperspective除去は
このraw像を六値方向・四trace方向・primitive像の積へ分ける。
Sp4と残るscale方向はprimitive像の側に保持し、さらに商を取っていない。
よってprimitive像次元は正確に25−6−4=15である。

実sourceが既約だから像閉包も既約であり、W16内の既約hypersurfaceとなる。
ただしその展開多項式、defining prime、他の既知hypersurfaceとの一致は未計算である。

既に独立確認した全新成分のprojective次元は25、rank4 S61切断の次元は25−9=16。
同じ全群と評価点を保つこの像が全成分の一次像を与えるので、
その15次元像上の一般S61 fiber次元は16−15=1となる。
fiberの被約性・既約性や全特殊fiber構造は結論しない。

review code `code/q_new_family_firstjet_review_20261103_v63.py`、
receipt `reports/Q_NEW_FAMILY_FIRSTJET_IMAGE_20261103_v63_REVIEW/`。
一般W16の全fiberを既知五点が尽くすかどうかは、依然別の未解決問題である。
