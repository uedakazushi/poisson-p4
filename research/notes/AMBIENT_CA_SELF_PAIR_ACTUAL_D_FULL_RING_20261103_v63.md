# CA自己旗族のactual交叉D：embedded成分を含む全完成環

2026-11-03。Dは実被約CAと実被約CBの**scheme intersection**である。
P_aの非零aのある非空開とそのG飽和S（cone20、projective19）で、全完成環を決定する。
これはB18原点の全Dを決定する結果ではない。

## 受理済みの幾何的入力

1. CAは二つのsmooth29形式枝を持ち、そのmeetingはC[[s20,rho]]/(rho²)。normalizationの全fiberを使用した全81 kernelのfiber-product記述がある。
2. CBはsmooth29。第一CA枝とCBの交叉はsmooth28のN、第二CA枝とCBの交叉はsmooth20のSである。
3. QQ(a)での全150行の真二次差から、Sの稠密開全体でT_CB⊂T_CA、従ってT_D=T_CB（29次元）と証明した。単一標本からの外挿ではない。

上記は `AMBIENT_CA_SELF_PAIR_MEETING_DOUBLE_20261103_v63_REVIEW` と `AMBIENT_CA_SELF_PAIR_ACTUAL_D_TWOJET_20261103_v63` およびその独立reviewを入力とする。

## 二枝とNの同時座標

第一枝をB1=C[[s20,rho,u8]]、第二枝をB2=C[[s20,rho,v8]]とする。
meetingへの両写像のkernelはそれぞれ(u8,rho²)、(v8,rho²)に直せる。
Nは第一枝内のsmooth hypersurfaceで、そのmeetingとのscheme交叉はSそのものである。
よってNの方程式はmeeting上でunit*rhoとなる。これを第一枝のrhoに取り、同じmeeting classを第二枝へ持ち上げれば、N=(rho)と二枝のfiber-product座標を同時に保てる。

T0=C[[s20,rho,u8,v8,tau]]で、第一枝は(v8,tau)=0、第二枝は(u8,tau−rho²)=0とする。
tauは(0,rho²)へ写る規約であり、先行reviewの逆規約からはtauをrho²−tauへ置き換える。
CAの全81関係は

    u_i*v_j, u_i*tau, v_j*(tau−rho²), tau*(tau−rho²)

である。単なる接錐近似ではなく全fiber-product kernelである。

## actual CB上の両ideal包含

CBの完成環をB=C[[s20,u8,f]]とし、N=(f)、S=(f,u8)とする。
N上でrho,v_i,tauは零なので、Bへの制限を

    rho=f*d, v_i=f*b_i, tau=f*c

と書ける。cはunitである。実際、d tauは二CA接空間の和37を殺すが、CBはその和を越える（rank38）。extra ambient equationsはT_CB⊂T_CAにより一次項を持たず、この判定を変えない。
従ってu_i*tauとtau*(tau−rho²)の制限は、それぞれf*u_iとf²のunit倍である。
actual D idealをQ⊂Bとすれば

    f*(f,u8) ⊂ Q。

逆にN⊂DなのでQ⊂(f)。F=f*k∈Qについて、**Sの稠密開全体**でT_D=T_CBだから、dF|S=k|S*df=0となる。
dfは非零なのでk|S=0、すなわちk∈I_S=(f,u8)。これはminimal81関係だけでなく、全ての追加ambient equationにも適用できる。従って

    Q ⊂ f*I_S。

両者より**全idealの等号**Q=f*(f,u8)を得る。

## 全完成型と不変量

cone完成環は

    D_hat = C[[s20,u8,f]]/(f²,f*u1,...,f*u8)。

reduced ringはC[[s20,u8]]（smooth28のN）、nilradicalはC[[s20]]*fでsquare-zeroである。
dim28、edim29、depth20、Hilbert–Samuel multiplicity1。
射影版ではs19として、dim27、edim28、depth19、多重度1となる。
多重度1だけからreduced/CMと判定してはいけない。

明示的primary decompositionは

    (f²,f*u8) = (f) intersection (f²,u8)。

従ってassociated primesはminimal Nとembedded Sの二つであり、S成分のnil moduleのgeneric長は1である。
Sは実G×A1像の閉包として代数的に定義され、非空開でこのlocal associated primeが得られる。
忠実平坦なcompletionとlocalizationにより、そのcone20（projective19）の閉包はglobal actual Dのembedded associated componentである。
R19は別のminimal componentであり、このSをRと同一視しない。

全Dのprimary decompositionやB18原点へのspecialization、CB自己旗族での別のembedded成分はなお別課題である。
