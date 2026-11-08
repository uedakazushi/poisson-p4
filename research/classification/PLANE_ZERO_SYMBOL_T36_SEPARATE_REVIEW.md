# 零symbol全39源のT36 actualcurve：独立の全係数手監査

2026-11-08。producer plane_collision_review、reviewer exceptional_coverage。
`PLANE_ZERO_SYMBOL_WHOLE_T36_ACTUAL_CURVE.md` と凍結receiptを独立に読解し **PASS**。
対象の全零symbol39-sourceのうち rank四・被約既約actual quinticを持つ全点は、既知閉成分T36に入る。
新成分増分は零。新ambient/Jacobi operator、点、minor、kernel、jet又は別算術runは行っていない。
独立手展開とproducerのexact計算receiptを区別する。

## 1. 真のcompleted sourceと全実逆

入力 `MPS_SMOOTHING_COMPLETION_20261103_v62/T36.json` の全八非零括弧を読み、
特に `Pbc=p0*bc+10*t2*t3*uv/p0` の必要correctionを保持した。
source `(a,b,c,u,v)`、target `(A,B,C,u,v)` に対し

```
p0=epsilon, t0=t1=3/epsilon, t2=t3=1,
A=a+(u+v)/epsilon,
B=b+5(v-u)/epsilon,
C=c+5(v-u)/epsilon.
```

mapはunit upper triangularの実GL5でdet一、逆は三shiftを引く。
左右両合成はidentity。全epsilon!=0で同じ一つのactual source pointと同じmapである。
追加tensor尺度は一、Euler gaugeは零。
period open `3p0²/25` はepsilon!=0上で非零、他のparameter反転は無い。

## 2. 全十raw括弧の独立手展開

Leibniz ruleと上の逆代入を各括弧に適用した結果を、次の全表として得た。
記号Q=u²-uv+v²をPfにだけ使い、括弧では元の全係数を残す。

| pair | epsilon零係数 R | epsilon一次係数 L |
| --- | --- | --- |
| A,B | -B(3u+v)/5 | 0 |
| A,C | C(u+3v)/5 | 0 |
| A,u | uv/5 | 0 |
| A,v | -uv/5 | 0 |
| B,C | B(2u-4v)+C(4u-2v) | BC |
| B,u | 3u²-2uv+v² | 3Bu/5 |
| B,v | 2uv-v² | Bv/5 |
| C,u | -u²+2uv | -Cu/5 |
| C,v | u²-2uv+3v² | -3Cv/5 |
| u,v | 0 | -uv/5 |

例えばABのepsilon^-1分子は
`3u²-v²-2uv+(v-u)(3u+v)=0`。
BCでは元のcorrection `10uv/epsilon` を含め、残分子は
`5(v-u)²+10uv-5u²-5v²=0`。
他八行も表示表の通りで、全負冪は消える。
各括弧の全二次式の等号なので、これは全150係数のraw等号
`M P(M^-1 X) M^t=R+epsilon L` を証明する。
Mとその逆が実逆であるから、この等号を逆変換すれば元全150係数を厳密に復元する。
単にD又はprojective classが一致したという推論ではない。

元truecompleted sourceの全projective Jacobiは既exact theoremを用いる。
epsilon!=0で実PGL像、かつ上表はepsilonの多項式だから、同じ全Jacobi恒等式はepsilon零にも延びる。
新しい350×150 operator又は別のfullJacobi runは不要である。

## 3. signed Pfと全lane仮定

座標順(A,B,C,u,v)のsigned Euler-Pfaffianを手で五項に分けて展開すると

```
F_R=-3A(u²-uv+v²)²+(3/5)B*u⁴-(3/5)C*v⁴.
```

例えばAのcofactorは
`-(3u²-2uv+v²)(u²-2uv+3v²)+(2uv-v²)(-u²+2uv)=-3Q²`。
残るEuler四項と合わせるとB係数は3u⁴/5、C係数は-3v⁴/5となり、
表示式以外の全degree五係数は零となる。
これはproducer receiptの全126係数と一致する。

B/C係数のgcdは一なので、F_RはC[u,v]上のprimitiveなbase三変数一次式。
C(u,v)上で一次式として既約であり、Gauss lemmaからC上絶対既約、従って被約である。
F_R非零なのでprojective rankは四である。

Rのbase-base行はbase-linear×normal-linear、mixed行はnormal-quadratic、normal-normal行は零。
従ってS=(u=v=0)はrank零plane、全normal weight零symbolは零、
Rは元のfull30 weight一方向だけに入り、full9 weight二方向を勝手に制限した一般源ではない。
P2=0の一つの実点として、既whole39 theoremのdomainに入る。
DのS沿い位数は四であるが、whole theoremは全multiplicityを持ち、この位数だけを対象へ仮定しない。

## 4. 一endpointから全actual39-sourceへ

既受理の原 `PLANE_ZERO_SYMBOL_WHOLE_SINGLE_ORBIT.md` とその別作者reviewを適用する。
全39の対象Poisson点について、anchor共通零はfactor/rank条件で除かれ、
conormal kernelは全fiber aff1、nonzero extension e in H1(O(-5))を持つ。
A=O(-1)^3なのでnonreduced binary cubic eはH0(A)=0に反し、三相異根型だけが残る。
その全型はnormal GL2とconstant base GL3のactual block PGL orbitで一つ。
CE transgressionがc1(O(-5))!=0なので全H2CEは零で、任意full9二cocycleは全global一cochain六係数のd像。
元Pi1がbase-linear、normal bracket零なのでfinite shearの高次項は無く、
全Pi2はactual linear shearとその符号反転逆で除ける。

従ってこれはformal bundle同型又はG上界ではなく、全sourceの一つのactual PGL orbitへの正逆である。
今回Rは同domain内の一つの実点だから、Rのorbitはその全唯一orbitに等しい。
epsilon!=0のP_epsilonはknown closed PGL-invariant T36に入り、
全A1曲線は非零projective tensorを与える（epsilon非零はrank四、零はF_R非零）。
閉性からR in T36、PGL不変性からその全orbit、従って対象全39sourceがT36へ入る。

あるknown boundary pointの所属から逆向きのmembershipを推論したものではない。
旧orbit23/tangent29からsmooth/unique componentも主張しない。
rank二・可約Dなどwhole theoremの外のsourceを一括して本命題へ含めていない。

## 5. 凍結exact receiptの監査と計算の区別

producer最終report SHA256 `27589fb06da531a9576cc2645e3edc8c37e933c04c193a3633749bc8a495d7bd`、
実source全文SHA、四input hash、凍結job SHA a08d8a0b…を照合した。
同reportの全150forward/逆・det1・全126Pf・Pf共変性・gcd・zero-symbol supportは全PASS、
実行0.189590秒/RSS49012KiBという**producerのexact算術receipt**である。
本reviewの算術根拠は§§1–3の別作者手展開であり、その実行時間を自分の別runと呼ばない。
rootの指示に従い、重複する別Fraction job/runを作成・実行していない。

初回Poly(0).terms()の末尾guard実装失敗recordは元bytesで保持した。
最後のzero skip以外の入力・曲線は同じであり、失敗を成功receiptに読み替えない。
過去の「全四double-zero leading quarticsが22型」という会話中診断にも依存しない。

新独立note/純receipt/manifestだけを作成し、producer/旧report/job/source/台帳/原稿/Gitは変更しない。
全対象laneのT36 membershipを受理し、新component数は零である。
