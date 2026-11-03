# B18におけるCA・N正規化sourceの全分岐scheme：独立監査

2026-11-03。対象は正規化射のsource上の最大Jacobian minor idealである。実targetの交叉Dやodd idealとは区別する。

親の有限正規化近傍を前提に、固定旗の20/19係数基底を物理的な二ベクトルから作り直す。Lie作用はproducerの微分関数を呼ばず、150係数上の添字置換と単項式微分として構成する。全16 parabolic生成元の全基底への作用、中心、非零orbit10 minor、実旗9方向を別に検証する。連結parabolicの作用と旗chartにより、source上で19個の滑らかな形式因子を除ける。

全150×9商行列と全150 Schur成分を再構成した。8次blockの行列式は全断面で1なので、残るSchur列は最大minor idealそのものを生成する。QQ上でproducerと異なるlex順序を使い、全Schur式と予測idealの両包含を確認する。小さい単項式primary交叉はLCM算法でも照合する。

CAのcone完成環は `C[[g19,w,u,v]]/(wu,wv,uv,u³)`、Nは `C[[g19,w,v]]/(wv,v²)` である。CAの被約分岐schemeは二軸のnodeでjunctionに三次nil、Nは滑らかな軸に平方零のembedded junctionを持つ。scalarを除くprojective版では滑らかな変数は18個となる。

CAのcone次元20・埋め込み次元22・depth19・多重度2、Nのcone次元20・埋め込み次元21・depth19・多重度1は、小単項式環と滑らかな因子から従う。これらの数値をCA/Nのtargetや実Dへ転用しない。

検証コードは `code/ambient_B18_CA_N_critical_review_20261103_v63.py`、結果は同名scopeのREVIEW JSONに保存する。
