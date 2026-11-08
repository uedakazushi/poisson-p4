# cusp接線・割線のwhole source：23次元と三つの24親又は新成分という有限還元

2026-11-08、PRIMARY E、`/root/log_manuscript`。純数学のproducer note。
新CAS、点、profile、operator、minor、normal jet又は既証明書の再実行は無い。
独立全読解は本稿とは別に行う。主稿/frontier/global台帳/Gitは変更しない。
本稿はknown-parent陽性曲線又は新成分の採録ではない。

## 1. 元の全sourceと全磁気の射程

座標(h,x,y,z,w)、base(h,x,y)、fiber(z,w)で、全sourceは

```
Pxz=-hx, Pxw=3x(h-x), Pyz=3(h-x)²-3hy, Pyw=-9xy,
Pzw=3xz+3R(h,x,y), Ph*=Pxy=0,
R=r0h²+r1hx+r2hy+r3x²+r4xy+r5y².
```

actual Dは全Rで正確に `-9h*x*C, C=h²y-(h-x)³` である。
`TWO_LINES_CUSP_TANGENT_SECANT_MAGNETIC_NORMAL_FORM.md` の全域actual fiber shearは
I=r5,J=r0+r2へ戻し、Iunitでは全lambdaの `R=lambda*h²+y²` へ正逆に移す。
ここで本稿のRは同noteの全括弧にある3倍を外に出した規約である。
I0は元A6の閉補集合として保持する。
この全lambdaの稠密像を一つのknown閉親へ入れれば、元A6の閉性からI0も全て戻る。
既whole lambda C39/C40陰性をI0へ逆に外挿しない。

全sourceのPGL5飽和像の被約閉包をXとする。
次節ではXの次元を点の接空間からではなく全fixed-divisor源から決める。

## 2. 固定actual Dの全源はGL2 frame＋六磁気

次のbase二次場を用いる。

```
U=-hx∂x+[3(h-x)²-3hy]∂y,
V=3x(h-x)∂x-9xy∂y.
```

全係数で

```
U(D)=-4hD, V(D)=(3h-12x)D,
det(E_base,U,V)=-D, [U,V]=-3xU.
```

Dはsquarefreeである。canonical divergence-free liftのbase関係を調べるときは
`A=U+4h E_base/5, B=V-(3h-12x)E_base/5` を用いる。
E,A,Bはlogarithmic derivationのbasisである：各primeの一般smooth点で
Cramer分子はそのprimeで割れ、squarefree性により全Dで割れるからである。
従ってAR(D)はdegree2のA/Bで自由、AR0=AR1=0、AR2=span(A,B)。

固定Dを持つ任意のrank四canonical projective tensorのbase-base skew blockは∇Dを殺す。
∇Dはprimitiveなので、そのdegree二係数vectorはdegree四の∇Dの多項式倍となり零。
mixed部分をfiber次数で展開すると、その一次/二次fiber係数はAR1/AR0に入り零。
残るmixed部分はA/Bの定数二行frameで尽くされる。rank四ならframeはGL2である。
全GL2逆とEuler gaugeでU/Vへ戻し、strict Jacobiのbase成分を比較すると
`K_z=3x,K_w=0`、従って `K=3xz+3R`、全R六係数となる。
これで固定Dの全源は上記元sourceと全fiber GL2のactual像で尽くされる。
一つの磁気代表から未知fiberを外挿したのではない。

fixed-D affine coefficient空間はGL2 frame四係数とR六係数の10次元、
projective scalarを除くと正確9次元である。
各frameのmixed係数とRの六係数は元tensorの異なる係数から一意に読めるので、
この9を単なるparameter数の上界だけとはしていない。

## 3. base divisor軌道14と全Xの次元23

Dは本質三変数なので、方向微分の二次元kernelからquotient planeをintrinsicに復元出来る。
そのkernelの選択はGr(2,5)の六次元である。
商平面内のDはcuspidal cubicと、そのcusp接線h、cuspを通る非flex secant xの組である。

cuspidal cubicのnormalization P1ではcuspと唯一のsmooth flexがintrinsicで、
それらを保つprojective automorphismはGmに入る。
非flex smooth点を通るsecantの追加markはそのGmを有限部分群へ縮める。
hとxはcubicへのcontact型(3)と(2+1)が違うため交換されない。
従ってbase PGL3 stabilizerは零次元、商平面内のD軌道は8次元である。
全PGL5 D軌道の次元は6+8=14。

Dへの写像の全fiberは前節の9次元fixed-D源で、元の全frame/shearはactualである。
よって全source飽和像は正確 `14+9=23` 次元、Xは既約23次元である。
旧M47/F48/Euler24/scale1の正確な実階数証拠とも一致するが、
新matrixやminorをこの論証に追加していない。
全lambdaでorbit次元一定とは主張しない。

## 4. 受理済み全source木から得る有限parent還元

ここでは、既存の全top-level source定理と全点合成（被約既約nonnormal五次の全被覆、
Hquartic全源の既分類と唯一未認識Ga9linear-w像、全nonplane二H+cubic、Q+C全支持等）
を原仮定のまま入力にする。Gのprimitive beta像非優越を全像上界へ置き換えない。
これらの有限source一覧の全像次元だけを用いる。

Xを含む全ambient成分Cを取る。Dが被約である条件はrank四開集合内で開なので、
Cの一般Dも被約である。三次以上の既約factorを持たないpartitionは、
全factorのproper multiplication像という閉条件を満たすため、Xへ特殊化出来ない。
したがってCの一般Dは5、1+4、2+3、1+1+3のいずれかである。

* 5：全pointwise定理によりknown成分。dim>=24ならT36,T37,T38,T39,T40だけが残る。
  その他の採録済み既約五次成分は23次元である。
* 1+4：全sourceのうちMPS該当成分/全新23以下の像/未認識Ga9linear-w像は23以下、
  唯一の24次元可能親はGa9cuspidal plane-quartic24である。
  その全componentは下記marked障害でXを含めない。
* 2+3：全sourceはknown23成分へ入るため、24以上のCは無い。
* 1+1+3：全nonplane源はknown MPS有限閉和、plane全配置の各全像は23以下。
  24以上でこの一般factor型を持つknown成分はT21だけであり、Xはその外である。

T39/T40はaccepted whole QQ[lambda] graphのconstant残式1/4,−2/3により全Iunit源の外。
T21は全二Hのproper flagとZ/Xの閉条件で全磁気源の外。
Ga9cuspidal24については、親simple H上のquarticは4pを持ち、その二次coneがH²の倍、
marked binary discriminantが三重根を持つ必要条件を使う。
本DでH=hなら残quarticは `x[h²y-(h-x)³]`、a0=0,B2=hxなのでdiscriminantはh²x²、
三重根が無い。H=xならrestrictionはh³(y-h)で第四乗でない。
この条件はproperなH/p/factor flagで閉じ、磁気係数を一切含まない。
詳しい正逆と全H選択は `GA9_CUSPIDAL_MARKED_DISCRIMINANT_FRONTIER.md` に固定する。

よって、

```
X ⊂ T36 又は X ⊂ T37 又は X ⊂ T38、
さもなくば X自身が23次元全Poisson成分である。               (D)
```

実際、三親以外のCはdim≤23であり、既約閉X⊂CとdimX=23からX=Cとなる。
これは新成分の採録ではない。最後の分岐には三親全閉包への非所属がなお必要であり、
それは旧listed uniform plane classの陰性からは得られない。

## 5. 真T37の一つの全境界sourceを使う際の恒等式

別のpositive construction入口として、completed T37でt2=0とし、
座標(A,B,C,u,v)、p=p0!=0、全t0,t1,t3を残す。
J=AB+t0*u²/(3p)と置くと、全raw八entryのEuler-Pfaffian五項を手で足して

```
D_T37 = p*u*[6p*C*v*J - 2t1*B*v³ - 3t3*J²].             (T)
```

を得る。特にt3=0なら
`D_T37=2p*u*v*(3pABC+t0Cu²-t1Bv²)`。
全correction `P_Cv=-3pCv+t3AB+t0t3u²/(3p)` を保持し、
最後の−3pt3*u*J²を落としていない。
これはT37全体の分類でなく、正確にこの三edge境界のfull source identityである。

(T)でdivisorをtargetへ合わせるだけでは足りない。
例えばsource normalをepsilon倍してtensorをepsilon倍する一様構成は、
対応planeのcoisotropyがなければnormal-normal bracketにepsilon^-1 poleを残す。
磁気R=lambda*h²+y²を残すには、その全pole消去と全二次係数での正則性を証明する必要がある。
旧uniform-plane finite classesのnegative scopeを増やすための新計算は行わない。
今後の最弱positive義務は、(D)の一親からgeneric lambdaの全sourceへ向かう
actual PGL/真period曲線を構成することである。generic lambdaの稠密包含が得られれば、
既知親の閉性と元A6源の閉性により有限lambdaとI0も戻る。
新しいanchor又は高次normal jetを必要条件にしない。
