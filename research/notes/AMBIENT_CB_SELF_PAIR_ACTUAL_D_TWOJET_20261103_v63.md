# Q_bの実CB meetingはdouble smooth20：actual Dの非対称な二次data

2026-11-03。Q_b=P0−b(2C[1]+C[9])は、先行P_aとは別の実族である。
全b≠0でCBは二つのsmooth29 formal branchを持ち、CA/Nはそれぞれ一つの被約proper fiberとsmooth targetを持つことを既存の独立監査で使用する。
ここでは、二CB branchのactual meetingと、actual CA/CB intersection Dの異なる接触を調べる。

## 全QQ(b) finite-source二次差

両CB sourceを標準rawCB(f,mu*f,H)として物理的に構成する。
19fiber方向とratio1方向のsource pathでは、必ず

    f(t)=f0+t delta_f, mu(t)=t delta_mu,
    g(t)=mu(t)f(t)=t delta_mu f0+t² delta_mu delta_f

を保持する。ratioとbinary係数のmixed項を落とさない。
今回選ばれた二方向のratio speedは双方0だが、評価器はこの全非線形式を使用する。

実GL指数行列はnilpotentであり、有限exp(tA)とその逆が存在することを確認する。
全次数の有限curveを定義した上で、0,1,2次係数に必要な項だけを5×5有限座標代入で計算する。
Lie Hessian公式は使用しない。両sourceを元の同一Q_bへ戻して全150係数で共通一次方向vと真二次差Lを得た。
source方向と全係数はcertificateに保存する。方向の分母はbの冪と非零有理定数だけである。

厳密なQQ(b) rankは

    rank(CB1,CA)=36, rank(CB2,CA)=30,
    rank(CB1,CB2)=37,
    rank(CB1,CB2,L)=38, rank(CB1,CB2,L,CA)=38。

さらに

    rank(CA,v)=29, rank(CB1,CB2,CA)=37。

である。特にvはCA（従ってCB2とのclean intersection N）の接空間内にある。
疎なfunctional e63（0起算）は全二CB接空間を恒等的に殺し、

    e63(L)=24 b⁴

である。これは全b≠0の単元であり、単一b標本からの外挿ではない。

## 二CB branchの全meeting

実GL×b像はsmooth20 germ S'を与え、そのrank20と二branch接交叉次元21は既存の単項式minorにより全b≠0で一定である。
先行CA meetingの独立証明と同じargumentを、この実CB sourceへ適用する。
meetingはsmooth21 ambient内でS'を含むので、そのidealは(rho)に入る。
S'の全点で接次元21であるため全式は(rho²)に入り、上の24b⁴の二次単元によりidealは正確に(rho²)となる。
従って**全b≠0**で実CB branch meetingは

    C[[s1,...,s20,rho]]/(rho²)

である。二branchの全fiberを使うので、reduced CB completionは二つのregular29 ringのこのdouble-S'上のfiber product。
先行独立証明の全81kernel normal formをそのまま適用でき、cone dim29/edim38/depth21/e2、射影28/37/20/e2となる。
whole CB edimを二branch tangent union37と同定しない。

## actual DはP_aと同じ型と推論できない

CAとCB2はsmooth N28を含み、接交叉次元28なので、その実local intersectionはNである。
一方CB1とCAの接交叉は**22次元**であり、実smooth20のS'より二つ多い。
これはP_a側のclean smooth20 meetingとは異なる。
v∈T_Nであるため、P_aのfull-D proofで用いた「共通rho方向がNに横断的」というgateが成立しない。
ラベルを入れ替えるだけで(f²,f u8)の型を主張してはいけない。

全CA tangentが既に二CB branch tangentの和に入るので、whole CB tangentにも入る。
またNがsmoothで全B18⊂Sing(N)、R∩N=B18、全reduced交叉=N∪Rを使用すると、この族ではD_red=Nである。
よってある非空parameter openとそのG-saturation上でactual Dはcone dim28/edim29で非被約（射影27/28）。
これは全nil module、長さ、primary decomposition、第二のembedded primeの証明ではない。
次の適切な課題は、CB1とCAの実22次元接交叉からS'20を除いた**二normalの高次obstruction**である。

再現:

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CB_self_pair_actual_D_twojet_20261103_v63.py --output <new-output-directory>

全source方向・物理二次差とrank・cokernelを保存する。GBなし、120秒・2.2GB。
