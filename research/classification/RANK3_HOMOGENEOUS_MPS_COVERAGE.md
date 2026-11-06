# 階数三二重二次曲面の頂点非包含配置を MPS 成分で覆う

2026-11-07。分類義務 `E_quadric3`。退化因子が実際に \(H+2\mathcal Q\)、
\(\operatorname{rank}\mathcal Q=3\)、H が二次曲面の頂点を含まない場合、
**全ての点が既知 MPS40 成分の和に属する**。
接触配置や完全分類は本稿の結論に含めない。

全逆形式の表示は `RANK_THREE_QUADRIC_SUPPORT.md` の全90係数計算に依存する。
以下では定数二形式の退化枝、共鳴パラメータ、実変形の特殊点を全て保持する。

## 1. 全点を一変数の純二次族の閉包へ還元する

H の補集合で \(Q=u^2+v^2+z^2\)、残る座標を w とする。
全ての閉逆形式は

\[
 \alpha=Qd\beta-dQ\wedge\beta,\qquad \deg\beta\le1,
 \qquad \operatorname{Pf}(\alpha)=kQ^2,\quad k\ne0
\]

である。\(\Omega=d\beta\) は必ず非退化である。
rank0 は Pfaffian 零を与える。rank2 の場合には
\(\iota_X\mathrm{vol}=\Omega\wedge\beta\) と置くと
\(X=Lx+a\in K=\ker\Omega\)、\(\operatorname{tr}L=0\)、\(XQ=-kQ\) となる。
Hessian の根基 R は一次元で L 不変、三次元商上の trace は \(-3k/2\)。
従って \(L|_R=3k/2\ne0\) なので \(R\subset K\)。商上の L の階数は高々1だが、
非退化二次形式に関する随伴を用いると \(\bar L+\bar L^\dagger=-kI_3\)。
左辺の階数は高々2、右辺は3で矛盾する。
符号と全点性の独立確認は `RANK3_HOMOGENEOUS_DEGENERATE_FORM_REVIEW.md` に記す。

\(\iota_Z\Omega=\beta\) とすると
\(\mathcal L_Z\Omega=\Omega\)、\(ZQ=\kappa Q\)、\(\kappa\ne1\)。
Hessian を H と置き \(S=\Omega^{-1}H\) とすれば
\([L,S]=(1-\kappa)S\)。全ての \(\operatorname{tr}S^m\) は零なので、
rank3 の S は Jordan 型 (4) の regular nilpotent である。
線形座標を選び

\[
 \Omega=da\wedge dd+db\wedge dc,\qquad Q=b^2+ac
\]

とできる。全解は

\[
 L=L_\kappa+xS+yS^3,\quad v=t\partial_d,\qquad
 L_\kappa=\operatorname{diag}\left(
 \frac{3\kappa-2}{2},\frac\kappa2,
 \frac{2-\kappa}{2},\frac{4-3\kappa}{2}\right).
\]

実際、regular nilpotent の symplectic 中心化代数は \(\langle S,S^3\rangle\)。
\(\exp(pS+qS^3)\) は Ω,Q を保ち、\(\kappa\ne1\) なので x,y を消せる。
さらに \(\kappa\ne4/3\) なら根基上の固有値が非零なので平行移動で t を消せる。
全源 \((\mathbb A^1\setminus\{1\})\times\mathbb A^3\) は既約である。
従って \(\kappa=4/3,t\ne0\) も、対角一変数族の PGL 飽和閉包に含まれる。
この閉包論法と全正規形の独立確認は `RANK3_HOMOGENEOUS_SYMPLECTIC_REVIEW.md` に記す。

\(A=(3\kappa-2)/2\)、\(B=\kappa/2\)、\(C=(2-\kappa)/2\)、\(D=(4-3\kappa)/2\)
と書く。h を Casimir とする対角族 \(P_\kappa\) の六係数は

\[
 (-Aa^2,\ 2Aab,\ (1-\kappa)b^2+Bac,\ b^2+Dac,\ Bbc+Dad,\ Cc^2-2Dbd)
 \tag{1}
\]

である。順序は (ab,ac,ad,bc,bd,cd)。全パラメータで Jacobi は零、
\(\operatorname{Pf}P_\kappa=(1-\kappa)(b^2+ac)^2\)。
以下でこの**族全体**を一つの既知成分に入れる。

## 2. κ=2 から SNC 点への実 Poisson pencil

\(P=2P_2\) とする。非零係数が

\[
 R_{ha}=-bh/9,\quad R_{ab}=b^2/3,\quad
 R_{ac}=ad,\quad R_{bc}=bd
\]

である二ベクトル R は、\([P,R]=[R,R]=0\) を射影二ベクトルとして満たす。
従って \(P+\tau R\) は全 \(\tau\) で実際の大域射影 Poisson 族である。
退化五次式は

\[
 -\frac h9\left(36a^2c^2-18\tau a^2d^2+72ab^2c+36b^4
 +12\tau b^3d-2\tau b^2c^2+\tau^2b^2d^2\right).
\]

\(X=-(b/9)\partial_a\) として \(R'=R-E\wedge X\) と置く。
これは同じ射影二ベクトルを表し、h を Casimir とする。
\(P+18R'\) に次の可逆線形座標を用いる：

\[
 A=a/12+b/4-c/3+d,\quad B=a+b,\quad C=-a+b,
 \quad D=-a/12+b/4+c/3+d.
\]

新しい六係数は厳密に

\[
 (8AB,\ 4AC-B^2,\ BC/6,\ 8BC,\ 4BD-C^2,\ 8CD). \tag{2}
\]

その Pfaffian は \((12AC+B^2)(12BD+C^2)/3\) である。
座標重み \((h,A,B,C,D)=(0,0,1,1,0)\) に対し、
係数単項式の重みから二つの微分方向の重みを引く規約を用いる。
対数項は重み0、三つの補正項はそれぞれ1,2,1である。
従って (2) の PGL 軌道閉包は、退化五次式 \(48hABCD\) の対数 Poisson 点を含む。
五つの座標超平面は単純正規交叉をなす。

MPS Theorem 1.4 の既確認の命題
「SNC log symplectic locus と交わる全 Poisson 成分は既知40個」
を適用する。\(P+18R\) を含む任意の成分 C は connected PGL5 の作用で保たれ、
その SNC 軌道極限も C に属する。従って C は既知 MPS 成分である。
既知成分へ特殊化することだけから一般点の所属を推論していない。

元の座標の重み \((h,a,b,c,d)=(0,2,1,0,-1)\) では P の全項が重み1、
R の全項が重み−1である。この作用を \(\rho_\lambda\) と書けば

\[
 \lambda\rho_{\lambda^{-1}}(P+18R)=P+18\lambda^2R\longrightarrow P.
\]

よって同じ閉成分 C は \([P]=[P_2]\) を含む。

## 3. 接空間の挟み撃ちにより族の全パラメータへ移す

\(P_2\) における**全150係数**の射影 Jacobi 微分行列は350×150であり、
標数零で階数は少なくとも102である。従って射影 Poisson スキームの接次元は
\(125-102=23\) 以下である。
ここでは h-Casimir や純二次という制限を接空間に課していない。

`MPS40_CROSSWALK_20261031.md` とその一次資料照合済み証明書によれば、
MPS40 の各成分次元は23,24,25のいずれかで、全て23以上である。
上の C が \(P_2\) を含むことから局所次元は23以上。
局所次元と接次元が23で一致するので、\(P_2\) の局所環は regular であり、
この点を通る全空間の既約成分は C 一つだけである。

族 (1) の像閉包は既約で \(P_2\) を含むため、全て C に含まれる。
第1節の全源の閉包還元から、共鳴値と全 affine 定数項も含め、
指定した rank3 配置の全点が C に含まれる。
滑らかな二次曲面の成分との交点を根拠に族全体を延ばした議論ではない。

## 4. 検証と射程

新計算は active task `E_quadric3` に束縛した有限 job として実行した。
全150接空間で実 pencil を見つける計算は約2.82秒、
有理座標変換と Pfaffian 確認は約0.20秒で完了した。
計算入力・checker 全文・全項恒等式・独立非零 minor は
`CLASSIFICATION_RANK3_HOMOGENEOUS_MPS_REVIEW_20261107.json` に保存した。
独立査読は `RANK3_HOMOGENEOUS_MPS_REVIEW.md` を参照。

独立な Fraction 疎多項式 checker は、全 κ の Jacobi、全 τ の pencil、
gauge、座標行列の行列式 −4/3 と変換後の全六式、二つの重み極限を確認した。
独立に再構成した350×150微分行列の102次小行列式は、
素数1000003を法として410150≠0である。この有限体計算は標数零での階数下界に用いる。
MPS40の全次元も再集計し、23次元29個、24次元8個、25次元3個と確認した。
独立検証は約0.13秒で完了した。

本命題は既知成分数を増やさず、残る rank3 接触配置と一般ファイバー G を解決しない。
