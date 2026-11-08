# rank一・長鎖translation・非定数g：全源像の次元25と既知129成分からの分離

2026-11-08、PRIMARY E_multiple_hyperplane。**全sourceの純幾何入口。全空間成分性はまだ主張しない。** 新ambient点、M1、minor、jet又はsolverは本noteの作成時には実行していない。

## 1. 全源と係数からの逆

`MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §5 の長鎖場

```
X=∂x+x∂y, q=(z,w), r=y−x²/2
```

を固定する。同節の非線形補助座標をPGL同値として使用しない。以下の六括弧は**元のaffine座標**の二次以下の全Poisson sourceである。

gを非定数affine式、U=u∂z+v∂wをaffine場、d=div U、Kを任意の二次以下の式とし、

```
U(g)=d*g, F=(U+d+1)K, H=gK
```

を課す。六括弧は

```
Pxy=−F−(1+d)r, Pxz=u, Pxw=v,
Pyz=xu+r*g_w+(gK)_w,
Pyw=xv−r*g_z−(gK)_z, Pzw=g.                  (1)
```

全項が元座標で次数二以下であり、定数／一次／二次項を切り捨てていない。

Jacobiを全parameterで直接確認する別の純証明もある。Jacobian一の補助変更(t,r,z,w)=(x,y−x²/2,z,w)では

```
P′=∂t∧Y+Jac(f),
Y=−[F+(1+d)r]∂r+U, f=g(r+K).
```

div Y=−1かつY(f)=−fなので、三変数Jacobian bivectorの恒等式
L_Y Jac(f)=Jac(Y(f))−(div Y)Jac(f)=0から全Jacobiが従う。Pf P′=Y(f)=−f。逆変更のJacobianも一だから、元Pfaffianも

```
Pf(P)=−g(y−x²/2+K).                          (2)
```

である。modular divergenceは(1)からXと直接読み戻せる。これらは任意二次Fを自由にした別族ではなく、上のFを全域保持した同じ源である。

係数表示はinjectiveである。Pzwからg、Pxz/PxwからUとd、PxyからFが決まる。Pyz/PywからH=gKの両一階偏微分が決まり、二つの候補の差は定数となる。しかし(U+1)H=Fgであるから、その定数差も零。g≠0なのでKも一意である。従ってparameter源の像に正次元fiberはない。

## 2. 固定modular場のsourceは既約・次元12

非定数gの空間はaffine三係数空間から線形部分零の直線を除いた既約三次元開集合である。各gについてU(g)=(div U)gはaffine Uの六係数への線形条件である。

qのaffine変更でg=zとすれば、全解は正確に

```
U=αz∂z+(β+γz)∂w,
```

であり、三次元である。任意の非定数gで同じrankなので、Uは基底上のrank三vector bundleを与える。Kは全六係数を持つ。従ってsource Sは既約で

```
dim S=3+3+6=12.
```

前節のinjectivityから実際の固定X像も次元12。上界だけのparameter countではない。Pf式(2)は全点で非零なので、全Sは一般rank四である。

## 3. Aff4軌道と全source飽和の次元22

順序(h,x,y,z,w)でXの斉次化はN=h∂x+x∂y、Jordan型(3,1,1)。NのGL5 centralizerは正確に

```
[[a,0,0,0,0],
 [b,a,0,0,0],
 [c,b,a,d,e],
 [f,0,0,g,i],
 [j,0,0,k,l]]
```

という十一係数行列で、det=a³(gl−ik)≠0。h covectorはim(N*)²という**固有の一次元線**であり、centralizerのh-rowは既にa*hだけである。h-row=hを課すのはa=1の一条件で、三条件ではない。従ってAff4 stabilizerは十次元、XのAff4軌道も20−10=10次元。

このstabilizerの元は元affine座標で

```
x′=x+b,
y′=y+b*x+d*z+e*w+c,
(z′,w′)=B(z,w)+(f,j), det B≠0.
```

特にr′=r+d*z+e*w+c−b²/2、q′はaffineである。gの非定数性とKの次数二以下を保ち、S全体を保つ。式(1)への逆も同じaffine群内にある。

Aff4×S→二次bivector空間の像をAとする。divergenceはconstant-Jacobianのaffine変換に対しequivariantで、像Pのmodular場はa_*Xである。従って各fiberのaはStab(X)のcosetに限られる。逆にそのcosetの全元でa^{-1}_*P∈Sとなるので、fiber次元は正確に十。よって

```
dim A=20+12−10=22.
```

この計算は全parameterのInfinitesimal rank又は標本からの外挿ではない。λXは(x,y,q)↦(λx,λ²y,q)でXにaffine共役だから、Aは全非零tensor scalar倍を保つ。

## 4. 唯一重複Hを動かした射影像Zは正確に25次元

g=zのchartで斉次化した実五次因子は

```
D=−h²*z*Q, Q=h*y−x²/2+K2(h,z,w).             (3)
```

一般KでQはrank五のsmooth quadricである。例えばK=z²+w²でそのrankは五になる。従ってsmooth-Q条件は空でない開集合で、h=0はDの唯一の重複既約因子、重複度は二。z=0は相異なる単純因子である。

marked Hを動かすincidenceは、PGL5/Aff4の四次元に上のaffine像22次元を加える。Poisson tensorのprojective classへ射影化してscalar一を除くと次元25。一般点でHは(3)から一意に復元されるため、marked Hを忘れるfiberは零次元であり、さらに次元を失わない。以上より実際のPGL飽和射影像の既約閉包Zは

```
dim Z=4+22−1=25.                             (4)
```

全rank四部分で、Dに重複一次因子を持つ閉条件が保たれる。generic Qのrank五・最大factor重複度二という非空開条件も保持する。gの線形部分零の補集合、定数gのH3共鳴、g零全源、shortchainはZへ含めたと主張せず、別sourceとして残す。

## 5. 既知129成分にZを含むものはない

ここで使うのは固定HEAD e2e718783b6811414a7f3c4b128a04a4503e5896 の `research/COMPONENTS.json` と現行 `KNOWN_COMPONENT_AUDIT.md` の既受理dimension/rankである。全129行のprojectionを同名prior129 inputに固定する。

* generic rank二の全六成分の全点はrank≤2なので、rank四のZを含まない。
* rank四で次元24以下の全114成分は、dim Z=25によりZを含まない。
* rank四・次元28のCA、CBは全rank四部分のDが五つの一次式の積である。これはprojective multiplicationのproper像という閉条件で、両普遍源の一般D=3H+H1+H2、H+2H1+2H2から閉包へ保たれる。rank五Qを持つ(3)はその外。
* 次元25のnil31-DDH、nil31-half、nil31-thirdも全Dが一次式の積であり、同じ閉条件で除外する。
* positive-weight NEW25は全Dに三乗一次因子を持つ。これもL³*Qのprojective multiplicationの閉像であり、generic(3)の最大重複度二と両立しない。
* 残るT01、T39、T40は次元25。もしZ⊂T_iなら、閉既約集合の同次元性からZ=T_i。しかしZの全rank四部分は重複Hの閉条件内にあり、T01の一般Dは五独立一次因子、T39/T40の一般Dは被約である。従って一致しない。MPS三成分の既受理被約開だけを用い、特殊period全体の追加認識を仮定しない。

個数はrank四が84×23次元、18×22次元、12×24次元、7×25次元、2×28次元で計123、rank二が六で計129。未知成分を次元26/27へ勝手に切る主張ではない。以上は**既知129にZを含む成分がない**という全prior比較である。

## 6. 次の有限判定に必要な最弱結論

現時点の結論は既約25次元の実Poisson閉族Zと全priorからの分離まで。Z自体が全空間成分か、未知のより大きな成分の真の部分かはまだ未証明である。

構造的に一般rank四・Q rank五である一つの固定有理点について全350×150微分がrank100なら、全Euler24とscalar一を除いた全周囲接次元は25。Zの局所次元25と一致し、その点はsmoothで全成分が一意、Zは新しい25次元全成分となる。この一回のfullambient判定は、新M1やjet列ではない。

rank100が証明できなければ、別点の列へ移らず失敗値・全matrixを保存し、parameter体又は実変形の構造経路へ戻る。さらに成分性が得られた場合でも、whole g零／定数g三次共鳴／shortchainの帰属は別に証明する。
