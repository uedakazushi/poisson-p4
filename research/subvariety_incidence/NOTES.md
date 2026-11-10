# Poisson部分多様体incidence：現行研究ノート

2026-11-10 / v10。再開時はRESEARCH_POLICY.mdから読み、原稿の証明を全文重複させない。

## 完全等方Log点

q=u0²+2u1w0+2zw1、i_(P0)Ω=dq∧i_ξ(du0∧du1∧dz)、ξ=(1,1,-2)。固定Q,σの延長ファイバーはSpec C[a,b,c]/(ac,c²)。embedded cは全ambientでも二次Poisson変形へ延長できず、五係数Λ(E∧Jac)はその障害に-9、全150方向の線形化に0。ambientの一次Jacobi rank95、射影Zariski接空間30次元。これは非滑らかさを証明するがambient非被約性を証明しない。元証明はpapers/quadric_ambient_quadratic_obstruction.tex。

## 二つの二次接錐式と対数型の下界

新しい消滅汎関数は
Λ1(F)=[u0 z²w1]F0124−[u0³z]F0124+[u0u1zw0]F0124、
Λ2(F)=[u0u1²w0]F0134−(2/3)[u0u1²z]F0124。
軌道20、q変形4、ξ変形3、K1,K2,K3の30方向にEulerゲージ24とP0尺度を足すとrank55。30座標t0,...,t29上の二次式は**9*t21*t29, -4*t25*t27**で、互いに素だから接錐次元≤28。

pの対数二形式はω0=q*u0*u1*z*dlog(u0/u1)∧dlog(q³/(u1⁴z²))。Log(1,1,1,2)パラメータ源は次元28。その350×28微分行列の行番号(2,3,4,6,7,8,9,10,11,16,19,20,21,25,26,27,28,29,30,38,39,42,43,45,46,52,53,55)からなる28次minorは正確に**-57177414000000**。局所次元は28。P0+aK1はrank4でPfaffianがq*u1²*zの非零倍、したがってrank2対数型成分と少なくとも一つのrank4成分がpで交叉する。

## 残る研究上の問題

rank4成分ID、全接錐の既約分解、acのambient障害、ambient局所環の非被約性、特殊Q³全配置は未解決。次は退化因子q+2H(u1)+H(z)と既存147カタログの実sourceを比較する。証明はpapers/quadric_log_component_intersection.tex、計算はcode/subvariety_incidence/check_local_dimension28.py・local_dimension28.json。一般smooth Q3・del Pezzo・cubic等の有効な結果は現行主稿に保持。新規性・外部査読は未確認。GitHub main正本、PDF非追跡、Actions不使用、Libraryは実push失敗時の緊急退避に限る。
