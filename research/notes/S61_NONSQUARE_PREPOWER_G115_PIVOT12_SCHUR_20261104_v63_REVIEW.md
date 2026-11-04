# G115第十二pivot：独立source・整域・typed Schur監査

2026-11-04。新scalar candidate／D5展開なし。N6は別のQQ peerが1325 sourceと十五積から
全9047係数と非零を検証済み。ここでは残余成分との積恒等式、chart、source逆を監査する。

## 4＋4＋3 parent

parentの103 typed nodeも独立確認する。pair-firstにした行と列の交換数はともに16で偶数。
A4=[Hbar0,Hbar1,anchor22,B*tail2] のdetは−T、adjugate最後の行は−κ。
選択pure六列のpair部分はrank一なので、他四pure列をanchorと交換したcofactorをvとすると
pair後のpure blockはQ=A−u v/k、k=κ3となる。旧chartではkも積単元の補題から単元である。
全16自由多項式恒等式 Adj(kA−uv)=k² AdjReduced を別checkerで検証し、
detQ=−cprod*Saux/kから Q^-1=−AdjReduced/(cprod*Saux) を得る。
最終operatorにkの逆は必要ない。一般非可換2×2 block逆と後続3block・Schurの両側逆を
独立検証する。同じ11行・11列のpivot順を入れ替えただけなので、残余4×11は既受理の
KAPPA11 Schur complementと一致する。全source／target順と逆演算を照合した。

## 固定六minorとの実source恒等式

略称G115/G85/G20はcurrent15の11/7/2というprefix後の行で、裸の135 generatorではない。
実heavy係数はrho_i=c_i*a*F*C²倍である。保存actual表から三つ全てについて
g*N10,j−h*N6,j=0を厳密QQ多項式で確認する。従ってB3消去後の行は、D(g)で
(c10/g) Ksmall、Ksmall=g*current11/c10−h*current7/c6となる。

pairの選択六pure列はrank一であるため、anchorが非零の稠密開上で、pair消去による
この六列への効果はliteral current2を使う一行のpure消去と同じである。P4四行と
current2を置く5blockのdetをD5、Ksmallとtarget19を付けた6blockのdetをD6とすると、
bordered determinant恒等式はD6/D5を与える。current2を前へ動かす交換は四回で偶数。
従って実old4×11の[1,10]成分rに対し

    g*D5*r = c10*D6 = c10*a²*N6.

gとD5は非零であり、最初は稠密D(gD5)でこの比を導く。両辺自体はold localized domain
FULLU[Delta407^-1,Saux^-1,T^-1]の元なので、その整域全体の恒等式へ戻る。
anchor、g、D5を最終chartで別途反転する仮定は不要である。

## 新3×10とscope

D(N6)で上式右辺が単元だからg、D5、rは全て単元となり、
r^-1=gD5/(c10*a²*N6)を使える。新しい逆nodeはN6と既知FULLU単元c10*a²だけ。
D5はliteral5detの有限recipeであり、展開した値を持たない。

全101 typed node、source／targetのpermutation、Schur左右演算とその逆を照合する。
old15/22逆を合成し、FULLU∩D(Delta407*Saux*T*N6)上の同cokernel3×10を得る。
生存sourceはcurrent15 [9,12,14]、targetは元32 [1,3,4,6,7,8,9,14,17,18]。
old25→normalized57→canonical63と実eta/Cdelta cleanupも同じ演算子に明示され、
raw七heavy係数を零にする。残る七tailと全chart補集合を保持する。

rank3なら元actual cubic rank57が従う十分条件は維持されるが、rank3、全FULLU57、
quartic239や尽くしは未証明。新しいscalar source、jet、点、GB、rankを計算していない。

再現: `python code/s61_nonsquare_prePOWER_G115_pivot12_schur_review_20261104_v63.py`。
