# Poisson部分多様体incidence：現行研究ノート

2026-11-10 / v14。再開時はRESEARCH_POLICY.md。詳しい証明は新接錐原稿と同伴incidence原稿を参照。

## 模型と完全な二次障害

q=u0²+2u1w0+2zw1、i_(P0)Ω=dq∧i_ξ(du0∧du1∧dz)、ξ=(1,1,-2)。固定Q,σファイバーはSpec C[a,b,c]/(ac,c²)。全150二次双ベクトル係数の射影Jacobi350方程式の微分はrank95、Eulerゲージ24と尺度1を除く接空間は30次元。

既存30基底のs_i=t_(20+i), a=t27, b=t28, c=t29で、全二次初期イデアルは
I2=(s0*c-s2*a, s1*c, c*(s2+c), s3*c, s4*c+s6*a, s5*a, s5*c)。軌道20方向とbは二次関係に現れない。全255左核×465接方向対の行列はrank7。新producerは旧producerをimportせず全行列を再構成する。

## 平方零元と被約台

h=a*(s2+c)。g1=s0*c-s2*a、g3=c*(s2+c)とすると
h=c*(s0+a)-g1、h²=a*(s0+a)*g3-h*g1。
hはac単項式を含み、7式のどれもacを含まないため二次関係に属さない。全次数二初期関係を尽くしたので、**真の接錐でもhは非零平方零元**。単に二次近似の非被約性を観察したのではない。

sqrt(I2)=I2+(h)=P28∩P26∩P24、
P28=(a,c)、P26=(c,s2,s6,s5)、P24=(s1,s3,s5,s2+c,s0+a,s4*c+s6*a)。
I2:h=(s1,s2,s3,s5,s6,c)。二次模型のnilradicalは24次元線形空間上の平方零の巡回加群。真の接錐では高次初期式によりhのannihilatorや台がさらに小さくなり得るが、次数二のhは消えない。

## 重複度一と混合acの消去

対数型表示ω0=q*u0*u1*z*dlog(u0/u1)∧dlog(q³/(u1⁴z²))と、源の非零28次minor -57177414000000 を再計算した。局所次元28。唯一の28次元接錐候補はP28で、そこではs1,s5が単元になりc,aが消えるので一般長さ1、Hilbert–Samuel重複度1。低次元nilpotent構造はこの最高次係数に寄与しない。

混合acについて、明示されたR_ac（新原稿第6節）がdF_(P0)(R_ac)=-F(K1,K3)を全350係数で満たす。従ってC[a,c]/((a,c)^3+(c²))上の二次補正は存在する。固定ファイバーのac関係が周囲でも障害だという予想はここで否定される。純c方向のc²障害と全次数の持上げは別問題。

## 最大次元対数型成分の滑らかさと交叉

e(R)=1の結合公式より、唯一の28次元既約成分Zlogの局所整域はe=1。複素数体上有限型でexcellent・formally unmixedなのでNagata重複度一判定から正則。T_p Zlog=V(a,c)で、完備局所環はC[[t0,...,t26,b]]。その滑らかさと全Poissonスキームの特異性・接錐の非被約性は両立する。全25 gl5作用の斉次双ベクトル係数行列はQQ-rank21、スケールP0は像にあるためPGL軌道次元20。軌道はLog1112とlong rank4の共通部分にある。25次元long成分はLogに含まれないので局所交叉次元23〜24（上下限、正確な値は未決定）。証明は現行papers/isotropic_log_tangent_cone_ja.texの局所重複度定理の直後。新producer check_isotropic_orbit20.pyでGL微分を検算。

## 共通23次元族の具体的な構造

q_v=u0²+2u1w0+2zw1+v w1²、xi_t=(u0+t u1,u1,-2z)、A1=u1∂u0−u0∂w0、A2_v=(z+v w1)∂w0−u1∂w1に対し、P(λ,v,t,b)=P(q_v,xi_t)+λK_A1(q_v)+bK_A2_v(q_v)は全パラメータでstrict Jacobi零。h=u1のchartでPf=-3λ z(x²−2y+2zw+v w²)、modular場は−λ∂x+[−λx+(10+b)(z+vw)]∂y+(10+b)∂w。

λ≠0、D=λ²+(10+b)²v≠0なら、X=(−λx+(10+b)(z+vw))/D、Y=y/D、Z=z、W=w+(10+b)x/λは非線形ではなくアフィン線形。modular=∂X+X∂Y、g={Z,W}=6Z。全六括弧は長鎖非定数g源に一致する。λ=0の族は、epsilon≠0でxi_(t,eps)=(u0+t u1,(1+eps)u1,(-2-eps)z)、f=(u0−t u1/eps,u1,z,q_v)とresidue (2+b/2)r1∧r2−(3+eps)r1∧r3−eps r2∧r3のLog1112源に含まれ、epsilon→0で戻る。全10二形式係数の恒等式を確認。

gl5作用25方向と(v,t,b)の3方向の150×28係数行列はrank24、24次minor=10871635968、Euler gauge24方向との合計rank48。尺度を除く23次元の共通像が実際に存在する。数値下界23だけなら埋込み次元30におけるKrull交叉次元不等式でも分かるため、数学的貢献は具体的な両源への接続である。コード check_isotropic_intersection_dim23.py と同名JSON証明書に全入力を保持。

## 保持する成分帰属と注意点

H=u1=0のaffine chartでx=u0/u1、y=-w0/u1、z=z_hom/u1、w=w1/u1。
P_a=P0+aK1はa≠0でアフィン線形変換X=-x/a+10z/a²、Y=y/a²、Z=z、W=w+10x/a（det=-a^-3）により長鎖源g=6Z、U=(60Z/a²)∂W、K=-ZW/a²+50Z²/a⁴へ全六括弧が一致する。残余quadricの半Hessian det=-1/(32a⁴)。modular Jordan型は(3,1,1)→(2,2,1)。したがってR4_MODULAR_RANK1_LONG_NONCONSTANT_25に帰属。既存の成分性定理を入力に、28次元Log1112との交叉は非等次元・非CM、depth≤25。この結論は接錐計算そのものとは別の依存を持つ。

未解決：全初期イデアル、低次元接錐の真の成分、hの正確な支持、元の局所環の被約性、他のrank4成分、全次数変形。次元26/24をambient Poisson成分と呼ばない。新規性・優先権の網羅確認は未完了。先行Q3・cubic・del Pezzo・P3の数学は同伴原稿に保持。GitHub main正本、PDF非追跡、Actions不使用、Libraryは実push失敗時のみ。
