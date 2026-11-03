# actual pair rank-failureのQQ有限性：独立監査

2026-11-03。author packet `S61_NONSQUARE_PAIR_FINITE_BAD_LOCUS_20261103_v63`を受理した。
このreceiptは**有限性だけ**を証明する。pair failureの空性、proper finite Z_(p) model、全raw unitや元product41を含まない。
後続の強い証明が得られても、この凍結scopeは変更しない。

## 実804係数から全pair空間を復元

archiveと全input SHAを照合した。
四quadricの804項からA,Bとpair chartを独立に再構成し、A=x0=W=0を全係数で確認した。
chartの三変数(u,v_light,z)はx2=v_light、x4=z、u=x3+7*x1/6から戻せるので、
部分的なcutではなく実際の三次元pair空間全体を表示する。
B=−10a*P*J*zであり、Bは**自由変数のまま**である。

Q0,Q3を制限しP^4でclearした全12係数がliteral native sourceと一致した。
保存normalized sourceへの変換も、既知a,P,Jと非零QQ定数のrow倍だけであることを確認した。
第0行に対する10項cofactorのliteral QQ dotはa*P^4。
従ってこの行はa*P*J≠0のどの点でも零にならない。

## 自然parameterと全必要方程式

新座標t=a/P、v=3e+20t、H=P+8tの逆は

    P=H−8t, e=(v−20t)/3, a=t(H−8t),
    J=5H−2v−4, Jacobian=P/3.

既知単元以外を除かない実可逆変換である。
全12自然係数と、Q3からh0=(160t−15H−2v+31)/36倍のQ0を引く行操作を再検証した。
後者はdeterminant1のrow operationであり、h0をinvertしない。

u²列はq0=27f/40000、q3=−3E/32000、
f=(v+3)(3v+4)、E=3H(v+4)+5v−12となる。
他の三列とこの列のminorから得る三式について、除いた因子がtの非負冪と非零QQ定数だけであることを確認した。
全式は正確にalpha_i(H,v)*t+beta_i(H,v)の形である。
さらに三つのt消去minorと保存projection polynomialの間のf冪・QQ定数の全因子恒等式を確認した。
fを除くのは次の一般chartだけであり、f=0は別に全て処理する。

## 特殊二枝は空

f=0ではu²列の第0行成分が零になる。
rank<2なら、非零の第0行と全minorから第1行のu²成分も零、すなわちE=0が必要である。
従って(v,H)=(-3,9)または(-4/3,7/3)。

両枝について、残る五列の全10 raw minorを独立に計算した。
既知unit t、P=H−8tとQQ定数だけによる保存normalizationを全係数で照合した。
零minorは零のまま扱い、零定数をinvertしない。
各枝で保存QQ cofactorを独自のFraction sparse積で掛け、正確な和1を確認した。
従って両特殊枝にrank failureは存在しない。

## 一般chartのprojectionは有限

f≠0に限ってf因子を除いた二つの必要polynomial N(H,v),M(H,v)を、実minorから照合した。
独立のFraction polynomial積により

    C(H,v)*N+D(H,v)*M=R(v)

を文字通り確認した。Rは非零でv次数33。
さらにNをHのpolynomialとみた全係数と保存QQ[v] cofactorの和が1であることを確認した。
したがってRによりvは有限個であり、どのvにもNはHの零polynomialとして消えない。
各vでHも有限個となる。resultantの単なる戻り値や標本は証拠として使っていない。

## t-fiberが線にならないこと

C0=-beta_3/f、b3=beta_2/f、alpha_3の保存表示を全QQ[H,v]で照合し、

    U*C0+V*b3+Z*alpha_3=1

を独立のFraction多項式積で確認した。
もし必要三式の解で全alpha_iが零なら、三式から全beta_iも零になる。
f≠0でC0=b3=alpha_3=0となり、上の恒等式に矛盾する。
よって**rank-failureの候補点では**少なくとも一つのalpha_iが非零であり、各(H,v)にtは高々一つ。
alpha_iだけが全空間でunit idealを生成するとは主張していない。
beta_2に現れる非定数factorを点ごとにinvertしてもいない。

これでpair rank-failureの複素点は有限。
有限型QQ schemeであるからfailure locusの次元は0以下となる。
この段階では零点の不存在も、特殊素数でのunitからのQQ promotionも結論しない。

## raw defectへの限定された帰結

独立監査済みのgcd例外RAW-goodと、pair⇒fiveおよびconditional linear-syzygy lemmaを組み合わせると、
FULL U上のactual raw24×56 rank defectはこの有限pair failureに含まれ、次元0以下である。
元product41への17-kernel/coordinate bridgeは別の監査事項として保持する。
full57、quartic、ambient removalも未完了である。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/s61_nonsquare_pair_finite_bad_locus_review_20261103_v63.py --output <new-output-directory>

source chart、全恒等式の監査receipt:
`reports/S61_NONSQUARE_PAIR_FINITE_BAD_LOCUS_20261103_v63_REVIEW/`。
約3秒、制限120秒・2.2GB、新しいGB・素数・heavy solverなし。
