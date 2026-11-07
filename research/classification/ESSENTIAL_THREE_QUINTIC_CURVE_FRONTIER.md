# 本質三変数の既約五次：具体曲線の有限化と文献の確認範囲

2026-11-07、PRIMARY E。全Poisson源はfree plane quintic exponents(2,2)に限られ、
固定Fの全磁気項は実shearで除去でき、一つの23次元軌道になることは既済である。
ここでは曲線代表の具体化を進める。全rational cuspidal quinticの尽くしは、
未確認の引用を残したまま完了と呼ばない。

## 直接読んだ資料

1. Dimca–Sticlaru, *Free divisors and rational cuspidal plane curves*,
   arXiv:1504.01242v4。primary PDF SHA256は
   `dc121721ced3e6add92cdc956dd473f4d3954733a7ecfc4fb8342cf99954a86b`。
   Introduction (1.1)、Theorem2.8(i), pp.6–7、§3、Examples4.1,4.3,4.4、Theorem4.6を読んだ。
   https://arxiv.org/pdf/1504.01242v4
2. Flenner–Zaidenberg, *Rational cuspidal plane curves of type (d,d−3)*,
   arXiv:alg-geom/9709001v2。primary PDF SHA256は
   `c999bf5403f354e34247b405f719945a11e8c8d40f69f9aae320f7be8ae64f87`。
   Introduction Remarks(1), p.2 と Remarks(2), pp.10–11 のquintic部分を読んだ。
   同本文の主定理はd≥6・三尖点以上であり、これを全quintic定理と読み替えない。
   https://arxiv.org/pdf/alg-geom/9709001v2
3. T. K. Moe, *Rational Cuspidal Curves*, arXiv:1511.02691v1、2008年修士論文。
   著者がoverviewで新結果を主張しないと明記するsecondary資料である。
   PDF SHA256は `7abc47523b28e1288850b8f9b500fd914240a8fef3dbe4a7d062fe7ef759d0eb`。
   §6.1、Table6.1,p.78、§6.3の代表式を読んだ。
   https://arxiv.org/pdf/1511.02691v1

Moeは全quinticの射影分類として M. Namba,
*Geometry of projective algebraic curves* (1984), Theorem2.3.10, pp.179–182を参照する。
Namba原本文は今回まだ取得・照合していない。
PDFは文献読解用workに置き、数学証明書のprivate入力にはしていない。
必要な定義式と確認箇所は本稿および公開candidate inputへ記録する。

## 尖点重複度四の全族は一つのfree代表へ閉じる

この節はNambaのリストやMoeのparameter表示に依存しない直接証明である。
Fを既約五次、p=[0:0:1]に重複度四のunibranch cuspを持つとする。
適当な座標で接線をy=0にすると、F=y⁴z+f5(x,y)と書ける。
既約性からf5のx⁵係数は非零である。これを1へscaleし、xをx+λyと変えてx⁴yを消し、
zをbaseの一次式でshearしてxy⁴,y⁵を消すと、全族は

    F_{C,D}=y⁴z+x⁵+C x³y²+D x²y³

となる。C,Dの値をgenericに固定していない。
全二次gradient syzygy W=(A,B,H)をzの次数で展開する。

    A=A0+z A1+a2 z²,
    B=B0+z B1+b2 z²,
    H=H0+z H1+h2 z².

A0,B0,H0はbinary二次、A1,B1,H1はbinary一次である。
WF=0のz³係数からb2=0、z²のx⁴係数からa2=0となり、同じz²行から

    B1=−h2*y/4.

z係数のx⁵,x⁴yを比較するとA1=0であり、残るx³y²係数はC h2=0を与える。

C≠0ではh2=B1=0で、z係数はB0=p xy+q y²、H1=−4p x−4q yを強制する。
z零次のx⁶,x⁵yからA0=k y²となり、x⁴y²,x³y³係数は

    5k+2Cp=0,   3Dp+2Cq=0.

残るH0は一意である。従って全二次syzygy核は正確に一次元で、free (2,2)ではない。

C=0では

    B0=3h2 D*x²/16+p xy+q y²,
    H1=−4p x−4q y,
    A0=k y²,
    5k+9h2 D²/16=0,   Dp=0.

D≠0ならp=0、q,h2の二自由係数だけとなり、全核は

    U_D=(0,−y²/4,3D x²/4+yz),
    V_D=(−9D²y²/80,3D x²/16−yz/4,9D³xy/40+z²)

で張られる。直接

    det(E,U_D,V_D)=−9D² F_{0,D}/64.

従ってこの全枝はfree (2,2)である。
D=0も合わせて残すと、核は一次syzygy (0,y,−4z)のx,y,z倍の三次元で、free (2,2)ではない。
以上は全十八係数のz次数・binary次数比較であり、特殊C,Dの未処理補集合はない。

C=0,D≠0では、C上のDの三乗根を選んでbaseをdiagonalにscaleすると、
F_{0,D}は x⁵+x²y³+y⁴z へ移る。
従って全重複度四cuspidal族のうちrank四Poissonを許すfree (2,2)曲線は、
今回の最初の代表の一つのPGL3軌道だけである。
全frame・六磁気源も既済inverseでその一つのPGL5軌道へ入る。

Moe p.81 のC1C parameter式は、PDF自体に次数四のs²t²を次数五座標に混ぜた印字があり、
そのまま斉次parametrizationとして使えないことを画像で確認した。
本節はその式を黙って修正した証明ではない。

## 他の具体候補と、まだ必要な確認

以下はMoeのTable6.1のconfigurationと、読解済みの具体式を整理したものである。
この表だけをprimaryな全分類の代わりに引用しない。

| configuration | 具体入口・確認状態 |
|---|---|
| 一尖点(4,5) | 上の全二parameter証明でfree条件を尽くし、最初の代表へ認識済み |
| A12 | DS Proposition3.6(ii), Example4.4のk=0。具体式は `z(yz−x²)²−2xy²(yz−x²)+y⁵`。新しい全18核は未実行 |
| E8+A4 | DS Theorem4.6のk=2、`y(yz+x²)²−x⁵`。正の一次対称性を持つ別のmonomial代表もあるのでτ=12だけで採用しない |
| E6+A6 | DS Remark3.3、Example4.1のd=5、`y(yz+x²)²+2x³(yz+x²)−x⁴z`。新しい全18核は未実行 |
| A8+A4 | Moe C5の具体parametrizationとNamba参照を得た。primary分類・全核未照合 |
| E6+A4+A2 | DS Proposition3.4、Example4.3のd=5,a=2,b=1。FZにもunique代表との記述がある。具体斉次式・全核を残す |
| 3A4 | FZ p.2、pp.10–11はprojectively uniqueと明記。Moe C7に具体parametrization。全核を残す |
| A6+3A2 | FZの同箇所はprojectively uniqueと明記。Moe C8に具体parametrization。全核を残す |

一般のrational cuspidal曲線では同じconfiguration内に複数射影軌道又はparameterがあり得る。
全free locusが有限軌道であることは既済だが、その有限性だけで上表の一式ずつが尽くすとは推論しない。
τ=12でもfree (2,2)でない例があるので、全二次gradient核とSaito行列式を必ず確認する。

## 三尖点以上のprimary記述にある小さい誤りと補修

FZのquinticの記述は、三尖点E6+A4+A2、三尖点3A4、四尖点A6+3A2の三つだけとする。
しかしpp.10–11の排除case(5)、configuration A4+A4+A2+A2で
「dualはcubic」とする一文は正しくない。dualの次数は四である。
Moe Remark6.1.1も同じ誤りを明示し、局所dualの計算で補修する。

この排除自体は次で直接補修できる。
各A4 cuspは適当なprojective affine chartで

    x=t²,   y=c4 t⁴+c5 t⁵+… ,   c5≠0

と書ける。dualの局所座標は slope=y'/x'、intercept=y−x*y'/x' である。
c4≠0では、slopeの二乗でinterceptのt⁴項を消すとt⁵項が非零で、
dualは再びA4 cusp、δ=2となる。
c4=0ではdualの最低次数は3,5で、δ=4である。
従って二つの元A4 cuspはdual上でδ少なくとも二ずつを与える。
同じdual点に写る場合は枝間交点のδも加わり、下界は減らない。
dualは既約quarticで算術種数三なのに、全δが少なくとも四となり矛盾する。

この補修で当該caseの非存在は確かめられるが、FZの全Cremona排除・Nambaの全曲線normal formを
すでに独立再証明したとは拡張しない。
特にA6+A4+A2という数値分配は実在候補として増やさず、FZ/Moeの排除対象として扱う。

## 次の有限義務

最初の代表の全component/audit packetは凍結し、旧producer/reportを変更しない。
次に、未処理の具体候補を一つの有限入力へ整理し、全18核でfree条件とparameter補集合を処理する。
primary分類の射程確認を並行し、全E支持の尽くしに必要なrepresentativeと未処理normal formを明示する。
未知の代表点を増やして全分類と呼ぶ、同じ成分の点を増やす、無期限jetを続けることはしない。
