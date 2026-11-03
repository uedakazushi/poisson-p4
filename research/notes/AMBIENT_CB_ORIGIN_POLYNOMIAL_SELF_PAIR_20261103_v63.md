# B18原点を通る別の実族：CBの二つの旗

2026-11-03。先行P_a族はCAの二旗を持ち、CBでは一旗だった。
ここでは係数を変えた**別の族Q_b**を構成し、CBの実二旗を全150係数で検証する。
P_aの凍結済み全CB一旗定理は変更しない。

## 新しい多項式族

同じCA fixed20行列Cで

    Q_b=P0−b*(2C[1]+C[9])=∂y∧V_b,
    5V_b=(5x²+8b*hx)∂h−2b*x²∂x
          +(5hx+3b*xz)∂z+(5h²+5z²−7b*xw)∂w.

Q_bはQQ[b]上の全150係数を持ち、b=0でP0へ戻る。
25 divergenceと350 strict cubic Jacobi係数は全多項式として零である。
constant factor∂yを持つ階数二の族であり、新しい全Poisson成分を主張するものではない。

三つのDarboux一次式はx、K=x+bh−b²z、H=x+2bhで、それぞれ5Vのeigenvalue−2b,3b,8bを持つ。
この比率から候補を選んだが、以下の確認は接方向やtrace必要条件だけではなく、実際のCB sourceである。

## 二つの実CB source

第一座標変換G1は

    h_new=x, x_new=K, y_new=y, z_new=h, w_new=w, det G1=−b².

第二変換G2は

    h_new=H, x_new=x, y_new=y, z_new=z, w_new=w, det G2=2b.

両者はQQ[b,1/b]上で実可逆である。
各G_i*Q_bを外積G_iと実逆座標代入によって物理5×5行列として計算し、
canonical(rawCB(f_i,0,H_i))との全150係数恒等式を確認した。
source cubicは残る二変数のHamiltonian条件から正確に積分した有限式であり、jet打切りはない。

    f1=−b*h*x,
    H1=−3b*h*w*z/2−h²*w+z³/3+z³/(3b²)
       +(h−x)z²/b³+(h−x)²z/b⁴,

    f2=2b*x²,
    H2=−b*w*x*z+z³/3−h*w*x/(2b)+w*x²/(2b)
       +(h−x)²z/(4b²).

これらのh,x,z,wは各変換後の標準CB座標である。
元の座標に戻した実旗は

    first:  L=e_y, W=span(x,K), marked h=x,
    second: L=e_y, W=span(h,x), marked h=H.

b≠0ならxとx+2bhは異なるので、二つの旗は実際に異なる。
両表示は同じテンソルQ_bに属する。

## 全source微分・orbit・同一Q上でのcross-tangent

実20座標map canonical(rawCB(f,mu*f,H))を先に作り、muと三binary係数の積を保持して微分した。
19fiber+ratio1+全flag9の150×29微分は両sourceでrank29。
19次定数blockを除いた10次minorは

    first: 432*b^5/15625,
    second: 786432*b^10/15625.

従って両immersionは全b≠0で成立する。
GL25のrankは19、Q自身のscalarを加えても19、dQ/dbを加えると20。
射影G-saturationの実像次元は19（cone20）である。
選択minorはそれぞれ150994944*b^12/390625と−301989888*b^11/390625。

第二CB sourceの全20fiber/ratio接方向を、G1*G2^-1で第一frameへ実pushbackした。
両branchの同一Qでの接空間の和はrank37、交わりは21次元。
共通fixed19列を除いたcross minorは82944*b^25/15625。
これはwhole imageのembedding dimensionを37と同定するものではない。
intersection schemeのclean性も主張しない。

## 全proper fiberから得る一般点の二branch

既知の全proper CB模型のP0 fiberは一点支持でscheme length2。
propernessからP0を含むtarget openでsource mapをfiniteにし、finite moduleの上半連続性により全fiber length≤2とできる。
Q_bのparameter線のinverse imageはb=0を含むopenである。
そこでb≠0なら構成した二つの異なる旗があるため、全fiberは正確に二つの被約点となる。
両sourceのrank29とcross rank37>29から、reduced CBのcompletionは二つの異なるsmooth29 branchの実imageとなる。

    A=T/(I1∩I2),  T/Ii≅C[[29]].

全fiberを用いているので、このopenでは第三のformal branchを落としていない。
ただし、全b≠0の全旗分類は行っていない。

scheme intersection C'=T/(I1+I2)のembedding dimensionは21なのでdepth C'≤21。
exact sequence0→A→C[[29]]⊕C[[29]]→C'→0よりdepth A≤22<29。
従ってこの新しい19次元射影saturationの一般点で、reduced CBはnonnormal・singular・non-CMである。
射影branchはsmooth28、depth上限21。

固定されたrank2 quadratic-normalizer base R29の上でA=R⊕tIという既存局所表示へ移すなら、
二つの被約source fiberはdiscriminant ν≠0のself-pair部分を与える。
その場所ではAは二つのR-sheetのfiber productで、branch intersectionがR/Iになる。
従ってこの族のcone20像と一般点のintersection edim21は、原点conductorの別の計算へ渡せる。
本packet自体は原点のI、全conductor、CM性・multiplicityを分類しない。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/ambient_CB_origin_polynomial_self_pair_20261103_v63.py --output <new-output-directory>

全150係数、二source、全Jacobianとminor:
`reports/AMBIENT_CB_ORIGIN_POLYNOMIAL_SELF_PAIR_20261103_v63/`。
実行約3秒、制限120秒・2.2GB。GBなし。
