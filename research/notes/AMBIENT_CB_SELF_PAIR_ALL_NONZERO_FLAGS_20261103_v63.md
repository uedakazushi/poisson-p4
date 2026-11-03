# 新CB self-pair族Q_bの全CA/CB/N旗：全b≠0で判定

2026-11-03。新しい多項式族

    Q_b=P0−b*(2C[1]+C[9])=∂y∧V_b

について、先行packetのgeneric finite-open結論を、別packetで**全b≠0**へ強化する。
以前のP_a=P0−a*(3C[1]+C[9])族と混同しない。
旧P_aの全CA二旗・CB一旗という結論は変更しない。

## 全proper端点を覆う内在的行列

新族の全75×5定covector contraction行列を再計算するとrank5、minorは−384*b^4/625。
全175×5の(Euler∧v∧Q_b)行列はrank4、minor2b、核は厳密に<e_y>。
従って、全b≠0で定covector Casimirは存在せず、projective constant factorは一意である。

CAのf=0端点およびCBのa_CB=0端点は、全係数のまま二つの独立な定covector Casimirを持つ。
両方とも上のrank5から不可能。
残るCA g=0とCB b_CB=0のsource旗は、そのLが必ずe_yである。
よって全ての候補は四変数(h,x,z,w)のmarked ellと二次元Wで尽くされる。
これは小さい旗chartだけの判定ではなく、既知proper模型の全階数二端点を覆う判定である。

## 全不変一次式と必要source恒等式

新族のtracefree場は

    5V_b=(5x²+8b*hx)∂h−2b*x²∂x
          +(5hx+3b*xz)∂z+(5h²+5z²−7b*xw)∂w.

ellにw項があると、Darboux cofactorのw²,hw,zw係数がそのcofactorをxの倍数に制限し、
h²係数がellのw項を零に強制する。
残る三変数ではUFDの因数性とeigenvalueにより、全marked ellは正確に

    x             : lambda=−2b,
    K=x+bh−b²z   : lambda=3b,
    H=x+2bh      : lambda=8b

である。三eigenformのdetは−2b³。
全一般Darboux方程式とeigenbasisを証明書へ保存した。

標準raw式をcanonical化し、全16Hamiltonian係数のまま確認した恒等式は

    CA: V(ell)=f_k*ell/5,       V(k)=−f+f_k*k/5,
    CB: V(ell)=2f_k*ell/5,      V(k)=−f+2f_k*k/5.

以下ではkはWの第二生成元である。

## CBは厳密に二つの旗

ell=KまたはHではx∈Wとなりk=xを取れる。
CBの二恒等式からf=lambda*x²/4+alpha*ell²となり、V(x)=−2b*x²/5との比較でlambda=8bだけが可能。
従ってHは先行の第二CB旗を与え、Kは不可能である。

ell=xでは

    f=alpha*x²−b*x*k,
    V(k)−3b*x*k/5 ∈ <x²>.

x²成分を除く全9×4行列はrank2、kernelはspan(x,h−bz)=span(x,K)。
従って第一CB旗のW=span(x,K)だけが可能となる。
先行packetは両旗に対して実QQ(b) sourceを全150係数で構成している。
marked h=xとmarked h=Hは全b≠0で異なる。
ratioは両方[1:0]で、固定旗19係数行列がinjectiveなのでsource係数の余分な解もない。

従って全CB geometric fiberは厳密にこの二旗だけ。
先行の両29次source微分minor432*b^5/15625、786432*b^10/15625が全b≠0で非零なので、
全scheme fiberも二つの被約点である。

## CAとNはそれぞれ一旗

CAのell=K,H候補では、必要式f=lambda*x²/2+alpha*ell²とV(x)の比較が、
lambda=4b/3を要求する。このeigenvalueは存在しない。

ell=xの場合だけが残り、

    f=alpha*x²−2b*x*k,
    V(k)−8b*x*k/5 ∈ <x²>.

対応する全9×4行列のrankは2、kernelはspan(h,x)。
従ってCAの全旗は、実際の元のfixed-frame旗L=e_y,W=span(h,x),marked h=xだけである。

このCA sourceはC[2]=0なので、N closed19-coefficient sourceにも入る。
同じ旗base上のN sourceはCA sourceの閉subbundleであるから、CA一旗の全判定はNも一旗であることを与える。
元の物理150係数でCA全150×29・N全150×28微分を再構成し、固定fiber block後の9次minorを計算した。

    CA: −262144*b^9/15625,
    N:   1572864*b^9/78125.

いずれも全b≠0で非零。
properness、一点支持fiberの接空間零、finite Nakayama同型により、全b≠0でCAはsmooth cone29、Nはsmooth cone28である。
射影ではそれぞれ28、27次元。

## CBの二branch帰結の全parameter化

全CB fiberが二被約点と分かったため、先行のgeneric two-branch証明は全b≠0へ適用できる。
CB reduced completionは二つの異なるsmooth29 branchの実imageである。
先行の同一Q上でのcross-tangent unionは全b≠0で37、intersection tangent次元は21。
従って実branch intersection C'のedimは21。

    A=T/(I1∩I2),  T/Ii=C[[29]],
    0→A→C[[29]]⊕C[[29]]→C'→0

からdepth A≤22<29となる。
ゆえにreduced CBは全b≠0でnonnormal・singular・non-CM。
射影branchはsmooth28、depth上限21。
whole CBのembedding dimensionを37と同定せず、intersectionのclean性も主張しない。

この結果は原点b=0の全odd ideal・全conductor・full completionやglobal bad locusを分類しない。
新しい全Poisson成分も数えていない。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/ambient_CB_self_pair_all_nonzero_flags_20261103_v63.py --output <new-output-directory>

全QQ(b)行列とminor:
`reports/AMBIENT_CB_SELF_PAIR_ALL_NONZERO_FLAGS_20261103_v63/`。
実行約2秒、制限120秒・2.2GB、GBなし。
