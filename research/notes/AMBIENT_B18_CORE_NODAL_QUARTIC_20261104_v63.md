# 二つの実 core 初期式の節点四次曲線と次の exact membership 目標

2026-11-04。固定済みの二つの QQ 二次式 f,g だけを扱う。
新 source jet・core 三次係数・GB・beta 点は追加しない。

## 1. pencil と唯一の節点

座標は (x112,x117,tauA,tauB)。f は core row67、g は row72-row60 に由来する。
元 `CORE4_QUADRATIC_INITIAL` の全係数から対称行列 A,B を作る。
rank(A)=4、rank(B)=3、gcd(f,g)=1 であり

    det(A+tB)=(25*t^2+5*t+3)/691200.

有限二根は単純である。B の唯一の vertex は QQ 点

    v=(7/5,2/5,7/24,1),   f(v)=0.

従って C=Proj QQ[x112,x117,tauA,tauB]/(f,g) はこの点で特異。
他の特異点があれば、その点 q で勾配が従属なのである pencil 行列が q を kernel に持つ。
有限の単純 pencil 根では corank1 で adjugate は非零定数倍 q*q^T。
det の微分は q^T B q の非零定数倍となるが、q が C に属せばこれは零で、単純性に反する。
無限側は B の kernel の一点 v だけである。従って他の特異点はない。

tauB=1 patch で v を中心にする。f の一次項は 2*(Av) として

    Av=(37/960,-13/240,-7/40,3/160).

特に tauA の微分 -7/20 は単位。f を消した滑らかな曲面内で、g の局所二次接錐は
tauB の変分を零、tauA=(37*x112-52*x117)/168 とした g の二次式である。
添付の2x2 determinant は非零。従って v は幾何学的 ordinary node、delta=1。

## 2. 既約性と QQ 上の明示正規化

単なる pencil 判別式だけから既約性を結論しない。以下で全曲線を具体的に覆う。
b=tauB=0 を v の線形補空間とする。g の制限 Q は rank3 の滑らかな平面 conic。
その QQ 点

    w0=(8,1,7/6,0)

を直接確認する。これは固定 quadric の補助点であり、新しい beta/source の点ではない。
T(w)=v^T A w とすると T(w0)=1/20。

Q 上から C への写像は

    Phi(w)=2*T(w)*w-f(w)*v.

Q と T=0 に f=0 が同時に生じないことを、二つの binary quadrics の非零 resultant で確認する。
従って Phi は basepoint を持たず、Q 上の全域射となる。
任意の C の点 p!=v を p=w+lambda*v と書けば、g(w)=0。
T(w)=0 なら f(w)=0 が必要で、上の resultant に反する。
よって lambda=-f(w)/(2*T(w)) が一意であり、Phi は C\{v} を全て覆い、射影 v への逆を持つ。
したがって C は幾何学的既約である。

gcd(f,g)=1 から完全交叉は純一次元・CM。
v 以外の一般点が滑らかなので generically reduced、従って全体が被約となる。
Q は滑らかな conic で QQ 点を持つので Q≅P1_QQ。
w0 からの線による二次パラメータ表示を Phi に代入し、四つの binary quartics を literal に保存した。
それらの gcd は1、係数行列 rank4、f(Phi)=g(Phi)=0 を全係数で検算した。
写像は有限双有理であり、C の正規化は P1_QQ、次数は4。
二つの node preimages は T(w)=0 の相異なる二点。
完全交叉(2,2)の算術 genus は1、唯一の node の delta は1、幾何 genus は0。
二つの接線が QQ 上で分裂するとは主張しない。

## 3. 実 core への限定した帰結

二つの実 core 関数 k1,k2 は先頭式 f,g を持ち、残る五関数は三次以上である。
f,g が正則列なので、filtered Koszul の次数帰納により

    gr(QQ[[core4]]/(k1,k2))=QQ[core4]/(f,g).

右辺は整域であるから実二関数の商 S2 も整域で、二次元完全交叉・CM・multiplicity4。
これは五つの残り実関数が零という意味ではない。
全 core 商は S2 をその五関数で割ったものであり、ここから socle や深さは決めない。

## 4. 次の計算を係数数の増大へ戻さないための入力

同次多項式 h が (f,g) に属するかは、保存した四つの binary quartics Phi を代入して
P1 上の binary form が恒等零かで厳密に判定できる。
これは被約整域の正規化への単射であり、点サンプルではない。
node の二つの preimages を同一視する gluing も保存した二次式 T=0 で追える。

実 S2 の場合も、C の node を除く generic point を一つの関数体 chart とする。
原点 blowup の strict transform では、f,g の Jacobian rank2 の unit minor によって
二つの法線座標を actual k1,k2 の unit Hensel 系から解ける。
従って正規化 P1 の関数体 K=QQ(P1) と radial parameter rho に対する

    S2 -> K[[rho]]

の具体的な代数的 operator を構成できる。先頭写像は Phi による整域への単射なので、
この写像も単射。残る実五関数が S2 内で零かを調べる次の exact membership 目標は、
この一つの generic chart でその五代数的関数が恒等零かを確かめることである。
もし非零なら S2 の整域性から全 core の次元は高々1へ下がり、その最初の非零 binary section の
零点と node の二 preimages が局所 branch を調べる具体的有限候補になる。
ただし有限 jet が零だっただけではこの恒等零を証明できない。
node の二 chart の conductor/gluing を保ち、full core の H0 を決めるには別の比較が必要である。

この generic strict-transform operator は次の提案であって、本 packet では構築・実行していない。
高次 source 求解や core 三次係数への自動引き上げも行わない。

再現:

    .venv/bin/python code/ambient_B18_core_nodal_quartic_20261104_v63.py --output work/B18-core-nodal-replay
