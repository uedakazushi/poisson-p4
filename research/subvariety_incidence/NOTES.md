# Poisson部分多様体incidence：現行研究ノート

2026-11-10 / v11。再開時にはRESEARCH_POLICY.mdを先に読む。原稿にある証明の全文は重複させない。

## 完全等方Log点の既知の構造

q=u0²+2u1w0+2zw1、i_(P0)Ω=dq∧i_ξ(du0∧du1∧dz)、ξ=(1,1,-2)。固定Q,σ延長ファイバーはSpec C[a,b,c]/(ac,c²)。cの埋込み接方向はfull ambientでも二次までPoissonに持ち上がらず、五係数Λの値は-9。全150係数の一次Jacobi階数95、射影接空間30次元。さらに二つの独立消滅汎関数が30次元接空間上に9*t21*t29と-4*t25*t27を与える。Log(1,1,1,2)パラメータ源の28次minorは-57177414000000で、局所次元28を確定。証明はpapers/quadric_ambient_quadratic_obstruction.texとpapers/quadric_log_component_intersection.tex。

## 階数四成分への明示同定（検算済み）

H=u1=0のaffine chartでx=u0/u1、y=-w0/u1、z=z_hom/u1、w=w1/u1とする。P_a=P0+aK1の括弧は
{x,y}=a(y-x²/2-zw)、{x,z}={x,w}=0、{y,z}=-6z²、{y,w}=2zw-2x²+4y、{z,w}=6z。
modular場はM_a=-a∂x+(-a*x+10z)∂y+10∂w。a≠0でJordan(3,1,1)、a=0でJordan(2,2,1)。

a≠0では**アフィン線形**変換X=-x/a+10z/a²、Y=y/a²、Z=z、W=w+10x/a（Jacobian=-a^-3）で長鎖源の
g=6Z、U=(60Z/a²)∂W、d=0、K=-ZW/a²+50Z²/a⁴、F=-ZW/a²-10Z²/a⁴
に六括弧全部が一致する。残余二次式Q=hY-X²/2-ZW/a²+50Z²/a⁴の半Hessian行列式は-1/(32a⁴)≠0。したがってa≠0の全PoissonテンソルはR4_MODULAR_RANK1_LONG_NONCONSTANT_25の実sourceに属し、閉包でpも属する。**同じPfaffian因子型だけから帰属を推測したのではない**。本帰属を確認する独立SymPyコードと証明書はcheck_rank4_long_nonconstant_membership.py、rank4_long_nonconstant_membership.json。通常の日本語の証明はpapers/rank_four_long_component_membership.tex。

## 幾何学的意味と未解決箇所

既存カタログのrank2 Log1112（28次元）と長鎖非定数rank4（25次元）の双方がpを通るという二つの**成分性定理を入力**として、局所環は非等次元・非Cohen--Macaulay。これはambientのnilradical非零と同値ではない。すでに完全なrank4帰属の一例は分かったが、pを通る他の成分、全接錐の分解、ac障害、ambientの非被約性、特殊Q3の全層別は未解決。次は二つの既知成分のscheme論的交叉と、全ambientの追加二次障害を照合する。既知のその他のQ3・cubic・del Pezzo・P3の有効な成果は現行原稿に保持。新規性の一次文献調査は未完了。GitHub main正本、Actions不使用、PDF非追跡、Libraryはpush失敗時のみ。
