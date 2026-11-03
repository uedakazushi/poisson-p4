# 四線形変数構造とpair gateへの縮約：独立監査

2026-11-03。author packet `S61_NONSQUARE_FOUR_QUADRIC_BINARY_REDUCTION_20261103_v63`を受理した。
全804係数からの実恒等式と、pair gateが補助五二次式のrank5を含意する証明を確認した。
pair gate自体の全域成立や全QQ raw24/product41は、このreviewでは証明しない。

## 実線形式とminor

archiveと全input SHAを照合した。
実Q0,Q1,Q2,Q3から、既知unit a²JだけでQ2primeを作り、

    Q1=A*B/6075+W*C1,
    Q2prime=40*x0*B/6075+W*C2

を再構成した。W除算はliteral polynomial divisionであり、Wをinvertしていない。
C1,C2は実線形式で、

    D=40*x0*C1−A*C2

が成り立つ。従ってD∈(A,x0)である。
実線形行列

    [[A/6075, 40*x0/6075, W], [C1, C2, −B]]

の三つの2-minorは正確に−D/6075,−Q1,−Q2primeとなることを確認した。

## 六entryは全域で正確に四変数しか使わない

新座標(A,X=x0,W,B,S,U)の明示inverseを再構成し、det=1/(200*a*P³*J)を確認した。
全a*P*J≠0で実可逆であり、未知factorは使わない。
authorのlinear coefficient matrix積とは別に、実C1,C2へinverse座標を直接代入した。
両者はS,Uを全く含まず、A,X,W,Bだけの線形式になる。
全保存係数とも一致した。

六entryの新座標係数行列はS,U列が零であり、(A/6075,40X/6075,W,−B)の4次minorは
−40/6075²という非零QQ定数。
従って全coreでlinear span rankは**正確に4**であり、6ではない。
六独立entryを仮定するgeneric determinantal議論へこの行列を流用しない。

## Dのuniform binary制限

全B=W=0におけるDを実係数から計算すると

    D = −8*(3Pe+2P+12a)/(1215P²)*A²
        −256a²J/(243P)*A*X+1280a³J²/243*X².

保存された式と全係数が一致した。
frameの可逆性により、この制限でもAとXは独立座標。
X²係数1280a³J²/243が既知単元なので、Dの制限は全a*P*J≠0で非零である。
最初のA²係数や未知因子はinvertしない。

## pair gateから補助rank5を得る二段階の制限

pair gateはQ0,Q3を**A=x0=W=0、Bは自由**へ制限した二式の線形独立性である。
新しいcoordinate frameと全実式によって、この空間ではQ1,Q2prime,Dが零になることを確認した。
実pair制限の変数(B,S,U)と全2×6係数も証明書に保存したが、そのrank判定を新たに探索していない。

五式の定数係数関係を第一の空間へ制限し、pair gateを仮定するとQ0,Q3の係数が零になる。
残るQ1,Q2prime,Dの関係を、**第二の異なる空間B=W=0**へ制限すると、上の非零DからD係数が零となる。
最後にW=0でB*(c1*A+40*c2*x0)=0となり、B≠0とA,x0独立性によりc1=c2=0。
従ってpair gateは全a*P*J openで補助rank5を含意する。
第一制限にB=0を追加して変数を失っていない。

## raw linear-syzygy lemmaとの接続

既存のconditional lemmaの純代数部分も確認した。
線形係数関係sum Li*Qi=0を(B,W)で割り、Q0,Q3のcoprimalityを仮定する。
次数よりL0,L3は(B,W)の線形部分に入る。
W=0でBをpolynomialとしてcancelした後、pair制限へ落とすと、L0,L3のB係数は零となる。
Aとx0は独立なので、残るL1,L2は

    L1=40k*x0+d1*W, L2=−k*A+d2*W

の形となる。
元の関係でWをpolynomialとしてcancelすると、五二次式の定数係数関係だけが残る。
上のpair⇒fiveから全係数は零、従ってraw24×56はrank24となる。
Wをlocalizeした議論ではなく、polynomial domainでの正当なcancelである。

独立監査済みのglobal比例排除とfinite-Sigma Nakayama証明により、gcd gate失敗点はFULL Uで既にRAW-good。
したがって残る十分条件はpair gateだけで、別の五二次式rank gateを解く必要はない。
これはraw-bad locusがpair failureへ含まれるという縮約であり、pair gateが全域で成立したという結論ではない。
元product41へのkernel/coordinate bridgeやquartic、ambient removalも、このreviewから完了としない。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/s61_nonsquare_four_quadric_binary_review_20261103_v63.py --output <new-output-directory>

実frame、全linear係数、binary D、B自由のpair制限:
`reports/S61_NONSQUARE_FOUR_QUADRIC_BINARY_REDUCTION_20261103_v63_REVIEW/`。
実行約9秒、制限120秒・2.2GB。新しいGB、素数、数値標本は使用しない。
