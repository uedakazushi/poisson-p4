# Q_(8/9) に到るT38の明示四次体フレーム

2026-11-03。既知の全QaのT38所属定理を再計数せず、Q_(8/9)における**具体的な有限代数的sourceと固有座標**を与える。
全150係数の実canonical Poisson曲線を保存する。T38のsmoothing edgeは今回は全て零であり、
正weight6軸のT38所属や残るP2の実枝を証明したものではない。

## 四次体とperiod

K=QQ[theta]/(theta^4-50theta²+305)を使う。この多項式は5-EisensteinなのでKは次数4の体である。
実表示theta=sqrt(25-8sqrt(5))は一つの埋め込みを選ぶだけで、全恒等式はK内で確認する。

    rho=(33-theta²)/16,
    z=(theta-7)*(rho-1)/4,
    x=-(3rho+1)z-(2rho+1),
    y=-(2rho+1)z-(3rho+1).

rho²-rho-1=0、2rho(z²+1)+7z=0となる。d=(0,1,z,x,y)は全て相異なる。
全10差について逆元とQQへの非零normを保存する。

    alpha=x²+rho*z*(x+y)+rho*x*y,
    beta=-rho*(z+x+y)-x,
    gamma_ij=d_i*d_j*(d_j-d_i)*
      [alpha+beta(d_i+d_j)+rho(d_i²+d_i*d_j+d_j²)+d_i*d_j].

厳密なK計算で gamma14=gamma23=gamma34=0、gamma13=gamma12/2、gamma24=-gamma12/2を得る。
p=gamma12は非零で、Pf(gamma)=p²/4も非零である。
**periodはgamma自身ではなくそのPfaffian dual**なので、保存されたT38 period basis
(0,1/2,0,0,-1/2,1)のp倍と正確に一致する。固有値の標識置換は用いない。

## 実際の全150曲線と固有フレーム

s=3(rho-1)/2、mu=-4/(9s²(s+3))とし、Nの四つのsubdiagonalを
(1,s,(s+3)/2,-9/(s+3))とする。N_t=N+t*diag(d)、Wk=N_t^k-tr(N_t^k)I/5として

    Pi_t=mu*(alpha*t² W1 wedge W2 + beta*t W1 wedge W3
             +rho W1 wedge W4 + W2 wedge W3)

を構成した。行列の線形場をwedgeする規約は既存の実物理150係数規約と同じである。
全係数はK[t]内の有限多項式で、全25divergenceと全350strict Jacobi係数を物理変数で直接検査した。
t=0における全150係数はa=8rho/[9(rho²-1)]=8/9の物理的Qaと一致する。

S0は単位下三角で、i>jに対して

    S0_ij=-link_(i-1)*S0_(i-1,j)/(d_i-d_j)

と定める。S_ij(t)=S0_ij*t^(j-i)とすればN_t*S=S*t*diag(d)である。
U=S-Iは五乗零なので、逆行列はI-U+U²-U³+U⁴である。
両逆積と全四つのtracefree冪の対角化をK(t)上で確認した。

この座標でPi_t/(mu*t^5)はT38の三項log代表

    p*x1*x2 d1 wedge d2 + (p/2)*x1*x3 d1 wedge d3
      -(p/2)*x2*x4 d2 wedge d4

のcanonical代表と全150係数で一致する。raw log代表とcanonical代表の間のEuler差も、
P_can=P_raw-(Euler wedge div(P_raw))/5により直接検査する。
従ってt!=0では、凍結されたT38全smoothing sourceのedge=0部分の実際のPGL変換である。
単にT01の曲線をT38と呼び替えたものではなく、T38の明示period条件と全tensor同一視が根拠である。

この有限フレームは、coalescing eigenframeにsmoothing edgesを加える次の構成の入力となる。
そのedgeのvaluation、正軸の極限、P2の持ち上げ、未知成分の有無は今回の主張に含めない。

code: `code/q_T38_Qa8over9_compact_frame_20261103_v63.py`。
report/certificate: `reports/Q_T38_QA8OVER9_COMPACT_FRAME_20261103_v63/`。
