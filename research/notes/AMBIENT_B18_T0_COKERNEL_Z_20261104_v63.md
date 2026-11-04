# T0計算の有限cokernel Z：実際の両側annihilatorと全R0作用

2026-11-04。T0のkernelを素朴にArtin baseへ下げない。
代わりに、**異なるmodule**

```
Z=coker(D0 -> D4)
```

が正確に長さ6以下のArtin moduleになることを証明し、その実matrixを計算する。
T0の有限台・socle・元Dのdepthはこの計算から主張しない。
ただしZの正確な同定からimage B0のminimal generator数8を得るので、
mu(D0)=9との比較による**T0非零性**は最後に証明する。

## 1. Socle問題はD0で測れる

T1=T2=T3=0とdegree4安定化は独立受理済みである。
同じR0上で、actual T0はker(D0→D4)。D4はliteral W1 tailでdepth1を持つ。
像B0⊂D4もH0_m(B0)=0なので、0→T0→D0→B0→0から

```
H0_m(T0)=H0_m(D0).
```

従ってJinf∩C0のkernelを先に求めなくても、socleの有無はD0で測れる。
T0は正次数moduleへ掛けると正次数のw-torsionとなるため、全positive-degree idealに
annihilateされる。よってこの有限socleがあることと、full transverse Dのdepth0は同値。
socleが無ければ既知の上限と合わせてtransverse depth1となる。
ただし今回はその有無を決定していない。

## 2. 有限なのはZ

CA側ではdegree1以降のcoefficient ambientはH_A=R0+m9*tau_A、degree0は
C0_A=R0+I0_A*tau_A。
T1=T2=T3=0により、D1→D4は実際の同型である。
従ってZはH_A/C0_Aのquotientであり、二つのpositive odd classesで生成される。
各zero odd coefficient f_i in I0_Aについて

```
f_i*((j_e*tau_A)*c)=j_e*((f_i*tau_A)*c)
```

の右辺は同じdegree-zero imageに入る。従ってI0_AはZをannihilateする。

同じargumentをCBのactual homogeneous presentationへ適用する。
degree-zero imageは、共通physical targetのD0そのものであり、両側で同じB0。
degree1 imageもintrinsic D1=D4である。
従って**I0_Bも同じZをannihilateする**。trace basisを二族間で同一視していない。
full idealの全99+87を保持したactual graded quotientを用いる点が重要である。

既受理の全次数等式は

```
A3=R0/(I0_A+I0_B)
  =QQ[[r,s]]/(r²,rs,s²),
```

である。CA actual F9座標の先頭z7を零にし、最後のr,sを残したもの。
従ってZは二生成A3-module、長さ6以下である。
さらにD4のminimal9とC0→D4の原点rank7から、mu_R0(Z)=2が従う。
これはcokernelのresidue-field次元なのでTorの単射問題を持ち込まない。

## 3. 正当なArtin base changeと計算

同じhonest9 frameでD4をQ1の全1548関係から表示し、そこへC0→D4の九列を足す。

```
Z=coker([Q1 | C0map]:R0^(1548+9) -> R0^9).
```

Z自身がA3-moduleなので、この**完全cokernel presentation**をA3へtensorすれば
right exactnessにより同じZとなる。kernel T0やambient imageのtensorではない。
全Q1 firstjetは独立再構築済みのfull9データから再利用した。

C0 mapのfirstjetは、保存済みCA source2と、実zero7-tail prefixの式から作った。
source function bのj番目のtail quotientは、r,sと先行zero座標を残すため、
その一次係数は保存二次係数b_jkのうちk≤jまたはk≥7だけを使う。
全35個を回復し、元の原点値とも照合した。
四つの新CA g quotientはmapのdependency closureに一つも入らない。
旧CB g商の未知first derivativeは独立な四つの形式方向として保持したまま、
全map entryから消えることを確認した。g微分を0に仮定していない。
全九変数firstjetから1,r,sだけを残せば、A3上では全次数で正確な値となる。

free A3^9のQQ次元は27。各relation columnとそのr倍、s倍を入れて疎Gaussを行った。
quotientのbasis、全relation-space basis、r/s作用行列を保存する。
relation spaceがr/s作用で保たれること、作用行列が可換で全m²積が零になることも直接照合する。
actual physical R0九座標の作用も、同じCA source inverseの一次係数から移した。

## 4. 計算結果：ZはA3のcanonical module

QQ relation rankは24、従ってZの長さは**3**、minimal generator数は2だった。
quotient basisはu=H7、v=H8、e=s*H8で、作用は

```
r*u=(8/5)e,    r*v=-(2/5)e,
s*u=-(2/5)e,   s*v=e,
r*e=s*e=0.
```

二つのgeneratorとsocleのpairing determinantは**36/25**で非零。
A3のdual basisをeps0,epsr,epssとすれば

```
u -> (8/5)epsr-(2/5)epss,
v -> -(2/5)epsr+epss,
e -> eps0
```

は正確なmodule同型Z→Hom_QQ(A3,QQ)=omega_A3である。
従ってHilbert function(2,1)だけから型を推測したものではない。

R0上でomega_A3のminimal relationは全て一次であり、数は17。
七zero座標を二generatorへ掛ける14本と、r,sの三つのlinear relationである。
canonical dual basisでは後三本は `s*epsr=0, r*epss=0, r*epsr-s*epss=0`。
これらはm²倍も全て消すため、追加の二次minimal relationは無い。
Koszul7 tensor Hilbert--Burch(m²)のdualより、minimal free ranksは

```
2,17,64,140,196,182,112,44,10,1.
```

## 5. Tor rank16と、T0の非零なleading class

actual C0 mapには、honest D4 free9のE0,...,E6を文字どおり消す七つのunit columnがある。
よってfree9→Zのmapは、最後の二つE7,E8をu,vへ送るprojectionに取れる。
Zの17minimal relationは全て一次なので、Tor1(D4,k)→Tor1(Z,k)のrankは、
Jの全firstjetをこの二成分へ射影したrankそのものである。
例外的な二次relationはm*Rel(Z)へ入り、このTor1 mapへ寄与しない。

全九変数の射影18係数matrixを独立疎Gaussで測るとrankは**16**。
従って0→B0→D4→Z→0のTor完全列から

```
mu(B0)=(17-16)+(9-2)=8.
```

actual mu(D0)=9なので、surjection D0→B0が同型であることは不可能。
**T0!=0**であり、そのD0/mD0への像は一次元である。
これはC0→D4の原点rank7だけからkernel非零を推論したものではない。

その一つのvisible leading classも同じ有限計算で識別した。
C0のbasis(1,z7*tau,g_A*tau)で、D4 originへのkernelは
`g_A*tau`と`f60*tau+(1/2)f67*tau-f72*tau`の二つのconstant vectorで張られる。
各vectorのmap firstjetをTor1(Z)/imageTor1(D4)へ送ると、前者は零、後者は非零である。
従ってker(D0/mD0→B0/mB0)の唯一のlineは**g_A*tauのbasis class**。
actual T0 elementにはm-order補正が必要かもしれず、raw g_A*tau自体がkernelに入ると
主張したものではない。有限socleやcyclicityもまだ従わない。

## 6. 次の小さいmapping-coneへの接続

保存した九つの可換作用行列N_iを使い、
R0-free module R0 tensor_QQ Z上の`x_i*I-N_i`のKoszul complexがZを正確に解像する。
各変数がmonicなので、順に除くと残り変数上freeとなるためである。
今回長さは3なので初めのfree ranksは3,27,108,252。
上のcanonical-module同定を使えば、より小さい2,17,64,140のminimal解像度も明示できる。

D4の9/44/112解像度からこの解像度へ比較写像を作り、surjection D4→Zのshifted coneを取れば
B0=image(D0→D4)の解像度を得る。
次にactual C0の9/28/56 Koszul解像度からB0へ比較し、surjection C0→B0のshifted coneを
取ればK0=ker(C0→D4)の完全presentationになる。
そこへ612個のactual Q0 liftを加えればT0=K0/Q0を正確に表示できる。
z7/g_CAの全syzygyを保つので、abstract SyzJinfを要求し直すだけのcriterionではない。

この比較写像とT0のsocle計算はまだ実行していない。
Zが有限でもT0が有限とは限らない。
閉fiberの長さやmap原点rankだけでは隠れたsocleのorderは制限できず、
固定generator数のR/(x^n)でsocle orderが任意に大きくなる反例がある。
従って、この結果からblindなjet次数増加は行わない。

Code `code/ambient_B18_T0_cokernel_Z_20261104_v63.py`。
Report/folder `AMBIENT_B18_T0_COKERNEL_Z_20261104_v63`。
新しいsource IFT/jet、点、prime、GBは無し。
