# 実 proper-source E4 の十四変数十七多項式と明示 flag 因子

2026-11-04。一回の bounded builder。source 48 の暗黙係数関数を展開せず、
proper-source 線形束を用いる。これは E4 の全 scheme 模型であり D0 と同一視しない。

## 1. 完全な係数入力

前提 `PROPER_POLYNOMIAL_MODEL_PREFLIGHT` の二つの19次元線形空間、実 P0、
source gauge15 minor -2/3、同じ物理 cut15 minor32を保持する。
CA は中心+九 free fiber 基底。CB の相対旗五変数は

    ell=x+a*h+b*z+c*w,   k=h+d*z+e*w.

frame の determinant は1、inverse の parameter degree は2。
変換された全十CA列（中心一列+free九列）の zero36係数を求め、
固定CB十九空間への定数左逆 readback19 と完全 annihilator17 を適用した。

    transformed36 = CB36*readback19 + injection17*equations17.

この全360係数の恒等式、および divergence4による150係数への全再構成を検算する。
CB19像に入ることと17式の零は全 ideal として同値であり、点や tangent だけの判定ではない。
最終入力は17x10係数行列。各entryは五旗変数の QQ 多項式で、全17式は九CA変数に線形。
全十四変数での sparse 17式も保存する。

二つの実 unit15 section の形式逆関数定理から、source gauge incidence と
current physical cut15 の normalized tensor incidence は全次数で形式同型。
affine tensor scalar は GL4中心に既に含まれ、追加の projectivization はしない。
GL4全16方向を quotient したのでもない。

## 2. 五 pivot の事前検査

展開した17式の原点17x14 Jacobian を直接求め、前提と同じ rank5・minor8/3を確認する。
選択列が CAfiber と相対flag のどちらか、および選択五式の全五 pivot 変数に関する次数を記録する。
九変数へ移す五 unit IFT は可能だが、五つ全てを fiber の線形消去と見なしてはいけない。
joint degree が1を越える場合、rational12式を作ったとは主張しない。

## 3. Npair の六つの polynomial 旗方程式

CA の binary f を physical reference で

    f=A*x²+B*x*h+C*h²

と書く。A(0)=-1は単位、B,Cは独立 free fiber 座標である。
実 Npair は

    I_Npair=(C, B-A*a, b,c,d,e)

で与えられる。これは以前に検算した raw25 成分の単位変換に由来する。
C=0、B=A*a とすると f=A*x*(x+a*h)。CBの marked線を ell=x+a*h に取り、
Euler vertical を除くと正確に CB fixedflag 表示になる。
H は二つの binary 変数を混ぜても z/w を含む全十六三次空間内に残る。
この形では A を逆にすることすら不要で、全17多項式への substitution が恒等的に零である。

I_Npair は六つの独立 coordinate 条件。
選択五多項式の線形classesと、pivotに選ばれなかった一つの flag coordinate が
六 conormal classes を張ることを actual QQ matrix で確認する。
その flag の名前と非零6 minor は report に保存する。
全17式は I_Npair に属するので、Nakayamaから完成局所環で

    I_Npair=(選択五多項式, そのflag).

従って五 unit IFT 後の残る十二式は、**実際の未消去 flag coordinate**を全次数で共通因子に持つ。
Bruhat rankや二次共通因子から推論したのではない。
この packet はその十二式の有理表示や全係数除法をまだ構築しない。
次の計算はこの小さい polynomial 模型の row/flag構造を使って行えるが、
GB・local standard basis・高次 jet を自動追加しない。

## 4. 検算と範囲

全 matrix-entry 多項式の term数、次数、共通分母、分子bit数、出力bytesを report に記録する。
42840項の主係数 ceiling、全出力150000項・32MiB、60秒・2GBを上限とする。
source proper模型と同じ単位 gauge に基づく E4 表示が今回の結論。
graph 冗長性、full core のCM、socle、D0の深さは未決定のまま。

再現:

    .venv/bin/python code/ambient_B18_proper_polynomial_model_20261104_v63.py --output work/B18-proper-polynomial-replay
