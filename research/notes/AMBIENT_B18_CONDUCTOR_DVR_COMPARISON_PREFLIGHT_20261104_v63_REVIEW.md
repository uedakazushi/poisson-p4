# conductor DVR 比較 preflight の独立レビュー

2026-11-04。一次 packet `AMBIENT_B18_CONDUCTOR_DVR_COMPARISON_PREFLIGHT_20261104_v63` を監査した。
作者 builder を実行せず、固定された実座標と QQ 係数を独自に読み、以下の算術と数学を確認する。
実 DVR 作用・完全な K 表示・torsion・socle・B18 深さは未計算のまま受理する。

## 1. 同じ物理 R0 上の parameter

CA の F9 最後の二座標は、実 positive odd 係数 j46,j95 そのものである。
raw physical q46,q95 でも任意の source 変数でもないことを元ノードで確認した。
ell=5*(j95-j46) とし、既存二族の物理 R0 Jacobian を独自 Fraction 行列積で輸送すると、
CB の zero7 plane 上では

    ell=(35*rB-85*sB)/3+O(m²)

となる。保存済み二つの実 cubic の全10係数を独立に評価し

    gA3(1,1)=125/24,
    gB3(17,7)=-181500,
    gB3(1,7/17)=-181500/4913

を確認した。ell が零の tangent line への三次制限が非零なので、両 CI 曲線上で
ell の切断長は正確に3。曲線環は CM 一次元で、ell は parameter かつ非零因子である。
Weierstrass による monic 三次表示から D=QQ[[ell]] 上自由 rank3 を得る。
CB の高次座標加速度は、低次項零の g の三次係数を変えない。

## 2. 非 canonical な DVR 分解

既受理の完全列を

    0 -> T -> R0/(IA*IB) -> R0/(IA intersect IB) -> 0,
    T=A3^6 direct_sum omega_A3

とする。A3=QQ[[r,s]]/(r,s)^2 は同じ CA 座標である。
scheme union は二つの D-free rank3 環への単射で D-torsion-free、有限 D-module なので
自由 rank6。有限商 A3 は generic rank を変えない。

A3 の basis(1,r,s) 上で ell は -5*r+5*s。その乗法行列と双対行列を独自に作り、
どちらも rank1・square zero と確認した。従って三次元 block はそれぞれ
D/(ell^2) direct_sum D/(ell)。T は七 block である。
free quotient への D-module 拡張は分裂するので

    R0/(IA*IB) = D^6 direct_sum (D/ell^2)^7 direct_sum (D/ell)^7.

special fiber は20。環の積表や canonical な20 basis を計算したわけではない。
この有限 kernel の長さ21を、任意の比較 module の torsion 長の上界に使ってはいけない。

## 3. raw tensor と actual cokernel

個別像 Ai=R0 direct_sum Ii*taui が同じ R0 上に source-bound されているので、
raw tensor の自由四正規化への像は

    R0 direct_sum IA*tauA direct_sum IB*tauB direct_sum IA*IB*tauA*tauB.

従って raw cokernel は三つの quotient の直和であり、D 型は

    D^12 direct_sum (D/ell^2)^7 direct_sum (D/ell)^7.

最小生成数は26。raw kernel は IA tensor IB→IA*IB の kernel、すなわち
Tor2(R0/IA,R0/IB)。既受理の全 Tor 分解から21個の三次元 block、
(D/ell^2)^21 direct_sum (D/ell)^21、長さ63・D生成数42となる。
その actual D0 での像の零性は従わない。

実出力差 F=delta+bA*tauA-bB*tauB は raw 像の**部分代数**の元である。
そのため raw 像で割った加群上に乗法を誘導する。
独自に自由四代数の二つの乗法行列から三座標へ射影し

    [[delta,0,-nuB*bB],
     [0,delta,nuA*bA],
     [-bB,bA,delta]]

を導いた。三行は IA,IB,IA*IB を法とする。
source annihilator を保つことは bA in IA,bB in IB から従い、
二つの普遍行列の commutator の全 monomial が各行の ideal に属することも独立に確認した。
raw cokernel 自体は ring quotient とは扱っていない。

十二個全ての delta,bA,bB の元 source ノードと、68本の CB zero relation・612積の対応を確認した。
右完全性により actual C は Craw/(sum Fj*Craw)。
実 D-basis と作用が bind できた**後なら**、26自由生成、14 torsion relation、
12*26 action columns による26/326以下の表示を得る。
この実 basis や action の計算は今回の入力には無く、実行済みとは認めない。

## 4. kernel で必要な完全性 gate

K,C は IA*IB で零化され、R0 上有限なので D 上有限。
しかし ell=0 の matrix kernel をそのまま K/ellK と同一視できない。
P:R0^640→R0^9、Q:R0^48→R0^4、U:R0^9→R0^4 とし、UP=QB の実 lift があれば

    L=ker[U,-Q],
    K=L/(columns(P,B)+columns(0,Syz(Q)))

となる。L→K の射は (x,y) を [x] へ送り、その kernel を上記二種類の列が尽くすことから証明できる。
完全な Syz(Q) と L の完全な表示が必要であり、現在の二次初期式から補ってはいけない。
28 intrinsic 列には零 lift、他612列には CB の全68多項式と C0 九 basis の
有限 target Hadamard 差を使う recipe は正しい。ただし recipe の存在と実全 lift の構築は区別する。
完全な K cokernel 表示が得られた後に限り、IA*IB で tensor しても右完全性から同じ K を得る。

D 上の有限 perfect complex で H1=K,H0=C、他の homology が零のものが完成すれば、
D の global dimension1 の普遍係数完全列から

    H2(derived fiber)=K[ell],
    0 -> K/ellK -> H1(derived fiber) -> C[ell] -> 0

を得る。K に D-torsion があることと K[ell] が非零であることは同値。
また有限 DVR module では mu_D(K)>rank_D(K) と同値である。
source presentation 二つを単に ell=0 にした matrix は、この完全 complex の代わりにはならない。
この gate と C[ell] 補正を残した一次 preflight は数学的に正しい。
E4 の H0 も未決なので、K torsion 判定だけから full D0 の深さを断定しない。

再現入口:

    .venv/bin/python code/ambient_B18_conductor_DVR_comparison_review_20261104_v63.py --output work/B18-DVR-review

新 source 求解・高次 jet・GB・Smith・作者 worker は一切実行しない。
