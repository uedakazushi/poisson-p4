# B18 actual D の深さを四つの零weight有限加群へ還元する

2026-11-04。既受理の相対三次モデルを使う新しい構造的還元であり、未計算の
actual coefficient を零と置かない。全87 graph 差、特に六本の三次差を残す。
**全 D の depth19/20 または full nilradical 自体は今回決めない。**
結論は、正weightを無限に追う必要をなくし、四つの rank2 ambient module 内の
有限検査へ落とせることである。ungauged core、source jet、GB は再実行しない。

## 1. 共通の実切断と CA の正確な代数

共通 G19 を除いた actual 非被約交叉を `D` とする。元の cone の深さは
`19+depth D`、既受理の上下界から `depth D in {0,1}` である。

共通 physical R10 の正weight base coordinate を `w=u140` と書く。零weight base は

```
R0=C[[x27,x108,x109,x112,x113,x117,x118,x119,x61]]
```

という正則9変数環である。名前は中心化した実際の physical rows を示す。
CA の受理済み actual quadratic normalizer と actual image は

```
BN=R0[[w]] ⊕ R0[[w]]*tau,        tau²=nu in R0,
A =R0[[w]] ⊕ (I0+w J0)R0[[w]]*tau.             (1)
```

ここで tau はweight0、I0は八本のactual零weight odd係数で、J0は二本のactual
正weight odd係数をwで割った係数で生成する。これは別々の抽象99模型をCA/CB間で
同一視する議論ではない。(1) は共通physical baseを使った CA 内の記述である。

次を R0-algebra として normalizer の零weight部分内に定義する。

```
C0=R0 ⊕ I0*tau,
H =R0 ⊕ (I0+J0)*tau.
```

両者とも代数である。実際 `(a+b tau)(c+d tau)=(ac+nu bd)+(ad+bc)tau`、
I0およびI0+J0はR0のidealだから閉じる。また

```
A_0=C0,          A_n=w^n H  (すべての n>=1).     (2)
```

これは一次近似ではない。八本の零weight generator と二本のweight1 generator
に対する equivariant Nakayama と実際の全次数有限二次模型から従う。
従って A は completed graded algebra `C0 ⊕ w H[[w]]` に等しい。
H/C0 はJ0の二本をj1,j2として `j1*tau,j2*tau` の **高々二元**で生成される。

## 2. 全交叉を A の homogeneous quotient として保つ

共通minimal20 targetのCA kernel、CB kernel、全87 graph差という受理済みの
同じ ambient presentation を使う。まず CA kernel だけで割れば (1) の A。
残る CB kernel と全graph差の像が A の閉 homogeneous ideal Q を生成し、

```
D=A/Q.                                         (3)
```

全相対 theorem により Q のhomogeneous generatorsはdegree0,1,2,3から選べる。
CA自身の99関係は既に零で、CBの高々99関係はdegree<=2、87 graph差の内訳は
degree0/1/2/3で4/53/24/6。従って Q のgeneratorは高々186本でよい。
この個数は最小性の主張ではない。**三次generatorは六本のgraph差だけ**である。
lower-weight graph差と全個別関係を捨てて三次だけを検査してはならない。

degree d のgeneratorを `w^d c_(d,alpha)` と一意に書く。
d=0ならcはC0、d>=1ならcはHの元である。ここでwをD内の単元としていない。
分解 (2) の中で coefficient を指定した記法である。

`L_d` をc_(d,alpha)の C0-linear span、d>3ではL_d=0と置き、

```
Jinf=H*L0+H*L1+H*L2+H*L3  subset H.
```

Qのdegree nをw^nで係数化したR0-submoduleをQ[n]と書くと、正確に

```
Q[0]=L0 subset C0,
Q[n]=sum_(0<=d<n) H*L_d + L_n    (n>=1).         (4)
```

証明は homogeneous ideal の次数n成分を `sum Q_generator,d * A_(n-d)` と書き、
n−d=0にC0、n−d>=1にHを使うだけである。completed ideal でも、固定次数への射影は
有限和なので同じ式が成り立つ。したがって

```
Q[n]=Jinf                    (すべての n>=4).   (5)
```

特に `w:D_n→D_(n+1)` は n>=4で同型。未証明の saturation、関係の除算、
自由性、flatnessを使わない。

## 3. w-torsion は四つの有限 R0-module で尽くされる

`T=H0_(w)(D)={d : あるkについてw^k d=0}` とする。(4)–(5)から

```
T_0=(C0 intersect Jinf)/L0,
T_n=Jinf/Q[n]                       (n=1,2,3),
T_n=0                               (n>=4).     (6)
```

degree0のC0とdegree>=1のHを取り違えない。wの反復は、degree0についてC0→H、
それ以降はHの係数を同じまま送る写像なので、そのkernelがちょうど(6)である。
全Tはこの四成分の直和であり、**w^4 T=0**。

全moduleはR0上有限である。さらに係数 `(a,b)` による `a+b tau` の表示により、
`C0,H,Jinf,Q[n]` は全て **R0^2 内の有限submodule** として扱える。
HやC0自体をR0-free rank2だと仮定する必要はない。

Tをnilradicalと同一視しない。wは被約R枝上でも消え得るため、w-torsionには
被約枝由来の元も含み得る。この還元はfull nilradicalの決定ではない。

## 4. 深さ0/1の正確な四module判定

AはR0[[w]]上有限で、Dもその局所商なので `rad((m_R0,w)D)=m_D`。
これはR0[[w]]→Dの単射を仮定しない。最大ideal torsionはw-torsion内に入り、
Tはw^4で消えるので

```
H0_(m_D)(D)=H0_(m_R0)(T)
            =direct_sum_(n=0..3) H0_(m_R0)(T_n). (7)
```

例えばm_R0^aとw^4が元を消せば `(m_R0,w)^(a+4)` も消すため、逆向きも成り立つ。
有限R0-moduleのH0が非零なら、その有限長部分は非零socleを持つ。
従って次は同値である。

```
depth D=0
  iff どれかのT_nについて Hom_R0(C,T_n)!=0.
```

どれも零なら `depth D>=1`。既知上界1によりこの場合 `depth D=1` である。
元のcone深さはそれぞれ19/20、projective深さは18/19となる。

完全に有限なmodule-colon certificateとしては、`m0=m_R0`、ambientをR0^2として

```
S0=((L0:m0) intersect C0 intersect Jinf)/L0,
Sn=((Q[n]:m0) intersect Jinf)/Q[n]    (n=1,2,3)  (8)
```

を検査すればよい。四つ全てのSが零 iff cone depth20。どれかの明示的非零classは
cone depth19のcertificateとなる。colonは九つのregular parametersについての交叉であり、
未計算の式の有限jetを零とする方法ではない。

## 5. 最上段 T3 に残る高々12個の具体的 target

`Vlow=H*L0+H*L1+H*L2` とすると

```
T3=(Vlow+H*L3)/(Vlow+L3).
```

六つのcubic graph差のnormalized coefficientをc3,1,...,c3,6と書く。
H/C0がj1*tau,j2*tauで生成されるので、T3は高々12個

```
(j_e*tau)*c3,alpha mod (Vlow+L3),
      e=1,2; alpha=1,...,6                      (9)
```

で生成される。c3,alpha=a_alpha+b_alpha*tauとすれば、R0^2内のpairは

```
(nu*j_e*b_alpha, j_e*a_alpha).
```

六本のraw graph rowsは `9,39,54,84,99,144`。この12 pairの所属は最上段を
判定する具体的な次の coefficient test である。**Vlowを含む** quotientで行う。
ここが零でもT0,T1,T2を調べずdepth20とは言えず、ここが非零でもm0-socleを
調べずdepth19とは言えない。

必要なactual inputは、同じCA trace基底のnu、I0八本、J0二本、全Q generatorの
normalized pairsである。core gauge preflightは、zero15 group equationsと八base equationsの
23次unit implicit system、その残りのmonic quadratic factorにより、これらの係数を
**algebraic power series**として指定するrecipeを与える。九つのungauged raw函数の
p=0 rank7をrank8と誤認してpを反転する方法は使わない。

今回はその係数を新たに展開していない。適当な有限jetだけのmembershipは(8)–(9)の
全次数certificateにならない。実際のalgebraic functionsを保ったmodule identity、または
正当化された有限determinacy boundを次の計算の受理条件とする。

## 6. degree3を残す限り w^4 は一般には改善できない

この小さい例はactual Dの同定ではなく、仮定だけからw^3で十分とは言えないことを示す。
`R0=k[[nu]]`、tau²=nu、I0=(nu)、J0=R0とすれば

```
C0=k[[tau²,tau³]], H=k[[tau]],
A=C0⊕w H[[w]], Q=(w³ tau).
```

degree3ではQ[3]=tau*C0だが、degree>=4ではQ[n]=tau*Hとなる。
tau²はC0∩Jinfに入り、tau² notin tau*C0。
よってD内のclass[tau²]について

```
w³[tau²]!=0,       w⁴[tau²]=0.
```

三次graphを全て保持する一般の構造定理から得る一様上界4は、この意味で鋭い。
actual六graphが下位idealに入ると別に証明できればdmax<=2に改善でき、同じ証明で
w^3 T=0になるが、その所属は今回未証明である。

## 7. 保存・検証範囲

専用checkerは既受理相対modelの全input hash、正weight row140、零base9、
oddの8+2分布、三次graph6本を読み、(4)–(9)の有限schemaと上のmonomial例を保存する。
新しいsource計算・点・field・ungauged core・大きいGB・formal coefficient solverは無い。
深さ判定を実行済みだとは主張しない。

```
.venv/bin/python code/ambient_B18_four_module_depth_criterion_20261104_v63.py \
  --output /tmp/B18-four-module-depth-replay
```

code、報告、criterion-schemaは同名stemの専用保存へ置く。古い全source・paper・ledgerは変更しない。
