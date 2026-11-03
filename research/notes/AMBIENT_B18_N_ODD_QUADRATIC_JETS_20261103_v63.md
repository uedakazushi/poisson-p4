# actual N odd generatorsの二次項：七一次＋二つの新しいminimal class

2026-11-03。generic B18代表点で、Nの**actual target image**
A=R⊕I*t、B=R[t]/(t²−nu)を使う。
Rは独立受理済みの共通Y28 actual target底であり、全odd係数b_iはBの基底1,tに対して一意である。
ここではactual b_iの二次項を計算し、mu(I)≥9を証明する。
source critical idealの生成数を、そのままactual Iの生成数と呼ぶ操作ではない。

## source critical tangent plane

N19 fixed-fiber基底とflag9の順序は前稿と同じである。
19基底の番号は0始まりで、critical tangent plane上の実Poisson族は

    P(W,V)=P0+W*B2+V*((5/2)*B1+(11/2)*B8+B16).

これは先に全critical idealを計算した9断面の
t0=5V/2,t1=W,t2=11V/2,t6=V、他0という二次元planeである。
全source Jacobianはfiber19列とflag9列を合わせた150×28行列J(W,V)であり、各entryはW,Vにlinearである。

## 選択27 actual base座標を保つvector field

Nの27単元target行をr、source列22（旗E_xh）以外をcomplementとする。
Jcomplement=B0+B1、J22=b0+b1と書く。
全27行を保ち、列22成分1であるvertical source vectorのcomplement係数を次数2まで求める。

    c0=−B0_r^(-1)b0_r,
    c1=−B0_r^(-1)(b1_r+B1_r*c0),
    c2=−B0_r^(-1)B1_r*c1.

全150 actual target微分の次数0,1は零で、二次項は

    S2=B0*c2+B1*c1.

選択27行も二次まで厳密に零となる。
さらにbranch row27のS2も零であることを、全QQ係数で確認する。

## なぜS2がactual odd係数の二次項なのか

sourceのtracefree quadratic coordinateをtとする。
上記vector field Dはselected27 base coordinatesを保ち、Dt(0)=1である。
残るactual base座標の微分は、nu=t²とそのbase derivativeの単元倍から
D(row27)=unit*2t*Dtとなる。
W,V planeではその次数1,2が零なので、t自身のplaneへの制限は次数3以上である。

actual target座標Fi=a_i(R)+b_i(R)*tにDを作用させると

    DFi=Dt*(b_i+2t*(a_i)_nu+2nu*(b_i)_nu)

となる。b_iのsource一次項はこのcritical tangent planeで零であり、tは次数3以上なので、
次数2ではS2_iが**actual b_iの二次制限そのもの**になる。
even部分の微分やtrace shiftを無視した近似ではない。

## 七つの一次classと二つの新しい二次class

前稿で全150 odd係数のactual R28一次行列Aを計算済みであり、rankは7である。
独立rowは18,19,28,29,60,61,67である。
このplaneのR28への微分はrank2で、Aを掛けると全て零になる。
従って任意のm_R Iの二次項はこのplaneへ制限すると零になる。

row67は正確にdb67=−8*dnuを満たす。
chosen lift i_nu=−b67/8の二次制限は−9V²/4である。
ただしこのclassは、一次項0の別のodd generatorを加えると変わるのでcanonicalとは呼ばない。

全rowから上の七basis rowの一次結合をsubtractし、actual R firstjetが零となるgenerator g_iを作る。
その二次項には

    g46|plane,2 = −4*W*V,
    g72|plane,2 = −6*V²

が現れる。WV,V²の係数行列のdeterminantは24で非零である。
ゆえにこれら二classはI/m_R I内で、七つの一次classから独立である。
結論は**mu(I_N)≥9**となる。

この独立性はactual I内の元で検証したものであり、paired-flag arcの不存在から生成数を推測していない。
したがって以前の八生成ansatz二つはいずれも除外されるが、full Iの等号やexact depthは本稿単独では主張しない。

またm_A=m_R⊕I*t、m_A²のodd部分はm_R I*t、even追加部分nu*I²はm_R²に入る。
従ってedim A=dim R+mu(I)≥28+9=37、projectiveでは36以上となる。
これはnormalizer source微分のrankをtarget接空間次元と取り違えたものではない。

* code: `code/ambient_B18_N_odd_quadratic_jets_20261103_v63.py`
* evidence: `reports/AMBIENT_B18_N_ODD_QUADRATIC_JETS_20261103_v63/`

全150 Jacobian、vertical vectorの三項、branch零、全actual一次row結合、全二次generatorを固定する。
full odd idealを決めるNakayama比較が成立するかは、別のproof companionで監査する。
