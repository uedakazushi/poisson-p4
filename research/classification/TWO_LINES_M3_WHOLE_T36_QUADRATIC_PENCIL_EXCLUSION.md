# 二H+cubicの残m3二源：全order二零点のjetでwhole T36を除く

2026-11-08、PRIMARY E、`/root/quartic_sources`。全sourceの純手証明。新point、局所変形jet、operator、minor、profile、GB又はCASは実行しない。ここで使う二次jetは、凍結した全tensorの式を各projective zeroで座標変更して最低次数を取る幾何不変量である。positive parent、成分性又は全E完了はまだ主張しない。

## 1. 全磁気源と受理するgeneric開

`TWO_LINES_PLANE_CUBIC_THREE_MAGNETIC_QUOTIENTS.md` とその別作者査読を使う。base(h,x,y)、fiber(z,w)、h-Casimir表示で

```
P=U wedge partial_z+V wedge partial_w+(-d*w+R)*partial_z wedge partial_w,
R=lambda*h²+h*y.
```

二つの全源は

```
node:
 C=h*y²-x*(x-h)²,
 U=(0,2x(x-h),y(3x-h)),
 V=(0,-2xy,2x(h-x)-3y²), d=x-h;
cusp:
 C=h*y²-(x+h)³,
 U=(0,2x(x+h),3xy),
 V=(0,2xy,3(y²-(h+x)²)), d=x.
```

全E_base,U,Vはgeneric独立で、determinantはそれぞれ4h*x*C、6h*x*C。U,VはCを保つquadratic logarithmic fieldsである。全六磁気をここへ送るactual shear/GLの正逆を保持し、元J零境界を削除しない。

以下の個別点非所属はnodeでlambda≠0、cuspでlambda≠±1に対して示す。special lambdaに存在する追加zero-planeを一般点の陰性で消さない。

## 2. Whole T36のproperな必要条件

四次元接空間Vのquadratic homogeneous bivector Qに対し

`alpha_Q=i_(E_V wedge Q) vol_V`

をcubic Euler-horizontal one-formとする。独立二quadric F2,G2のpencilから `[G2*dF2-F2*dG2]` を作るmapは、Gr(2,Sym²V*)上のregular projective mapである。共通因子を持つpencilも含め、char零では非零である。proper Grassmannianの像R22は閉、そのaffine coneに零を加えたC22も閉GL4-invariantである。

全tensorのorder二以上のzero qをmarkするproper incidence上では、そのquadratic jetとalphaは接空間GL4と非零scalarだけで変わる。従って `alpha_Q in C22` を加えた条件も閉で、qを忘れた像はP4のpropernessにより閉である。

真T36の全period sourceについて、vertex q=[1:0:0:0:0]の四変数jetは `GA9_LINEAR_W_WHOLE_T36_QUADRATIC_PENCIL_EXCLUSION.md` §§1–3の全式(2)であり、

```
q+=(p/5)B*u+t2*v², q-=(p/5)C*v-t3*u²,
alpha_Q=(5/p)*(q-*dq+ - q+*dq-).
```

p≠0なので二quadricは独立で、t0/t1 rowsのaffine correctionはcubic以上。よって全T36 sourceは上の必要条件を持ち、全PGL像とその閉包、すなわち**whole T36**も持つ。この全period/zero-edge/任意moving basisを覆う閉性を使い、固定frameの陰性へ狭めない。

## 3. Generic lambdaでの全zeroとvertex外のfirst jet

`TWO_LINES_PLANE_CUBIC_ZERO_LOCUS_T38_EXCLUSION.md` の全zero方程式を手で再確認した。

h=1ではsectionはh-Casimir affine bivectorそのものである。nodeの全base common zeroは(x,y)=(0,0),(1,0)。後者のzero-planeはlambda=0だけであり、今の開では無い。前者はw=-lambdaのfiber lineで、Pxz=2x(x-1)のx微分が−2なので、その全点のfirst jetは非零。

cuspのcommon zeroは(-1,0),(0,1),(0,-1)。後二者のzero-planeはlambda=∓1だけ。今の開では前者のfiber lineだけであり、Pxz=2x(x+1)のx微分が−2なので、やはり全点がsimple zeroである。fiberの自由zやRの特殊値を固定してこの結論を得たのではない。

h=0,x≠0ではraw4×4 Pfaffianはそれぞれ−4x⁴、−6x⁴で非零。raw rank四はEuler wedgeというrank二bivectorにはならないのでprojective zeroは無い。

h=x=0,y≠0ではprojective zeroはz=0のlineである。このlineをy=1 chartで調べる。projective affine係数は

```
Pzw_aff=Pzw-z*Pyw-w*Pzy.
```

line上の任意点w=w0で、そのfirst jetは

```
node: h+2w0*x+3z,
cusp: h+2w0*x-3z.
```

特にz微分が±3で、全点がsimple。raw二vectorの一成分の零だけでprojective jetを決めていない。

残るzeroはvertex line `L={h=x=y=0}=P(span(z,w))` 全体である。従って今のgeneric lambdaにおける**全order二以上zero候補はLだけ**。次節はLの全qを両fiber chartで扱う。

## 4. 全qにおけるquadratic jetの正確な式

q=[z:w]=[1:kappa]でz=1 chartを取り、r=w/z-kappa、base座標もzで割る。projective affine bracketは

```
{base_i,r}=V_i-(kappa+r)U_i-base_i*(-d*(kappa+r)+R).
```

最低二次だけを取ると

```
Qq=Wkappa wedge partial_r,
Wkappa=V-kappa*U+kappa*d*E_base.                    (1)
```

base同士の括弧は零である。Rはbase二次式なので、その寄与base_i*Rはcubic以上。従ってlambdaはquadratic jetに現れない。raw fiber-pairのlinear-base項はEuler補正で(1)のradial項になり、projective first jetは零である。rawのfirst jetが非零でもprojective first jetを非零と誤認していない。

q∞=[0:1]ではw=1、r=z/wを使うと

```
Q∞=(U-d*E_base) wedge partial_r.                   (2)
```

従ってalphaに寄与するbase fieldは、全q=[zeta:omega]に対し非零scalarまで `zeta*V-omega*U modulo linear*E_base` である。特にどのqでもalpha_Qqは非零。もしある組合せがradialなら、E_base,U,Vのgeneric独立性に反する。q∞も含め、一つのfiber点だけを選んでいない。

## 5. 既約cubic invariantはquadric pencilと両立しない

Wkappaも(2)のfieldもCを保つ。r方向もCを保つので、非零alpha_Qqが定めるP3上のcodimension一foliationに、既約cubic cone C(h,x,y)=0がinvariantである。

今回はalphaとCの次数が共に三なので、degreeだけで「alphaはC上で非零」としてはいけない。全alpha係数がCで割れると仮定すれば、alpha=C*beta、betaはconstant one-formとなる。しかしEuler-horizontalから `C*i_E beta=0`。多項式環はdomainなのでi_E betaは零、constant betaの全係数も零となる。これはalpha≠0に反する。従ってC一般点でfoliationは非零であり、その接平面はCの接平面そのものである。

alphaが独立quadric pencil G2*dF2-F2*dG2に比例すると仮定する。rational first integral J=F2/G2はC一般点で定義される。degree二の分母はdegree三の既約Cを含めない。invarianceとfoliation非零より、dJはCの全接方向で零となる。char零ではJ|Cは定数cであり、F2-cG2はC上で零。degree二<三なのでこれは零polynomial、pencil独立性に反する。

よって**全q∈Lでalpha_QqはC22の外**。他にorder二zeroが無いので、§2のwhole T36必要条件を満たさない。

## 6. 受理範囲と残課題

nodeのlambda≠0、cuspのlambda≠±1の**全点はwhole T36の外**である。各全六磁気源のPGL飽和閉包SiがT36へ全体として入ることも不可能である。

special lambdaの追加zero-plane、元J零磁気境界、又は別の同次元成分との個別交差へ陰性を外挿しない。T38については別のwhole zero-surface障害がある。従来prior126のwhole-parent還元と併用した場合、第一二sourceのそのcatalogue内の残whole候補はT37だけとなるが、これはT37所属を証明したことではない。新catalogue entriesへの比較は統合者が別に確認する。

固定1PS、normal rootのalignment、conductor normal multiplicityの保存、flat normalization又は全sourceの未知classificationは本証明の仮定にしていない。新成分数、局所scheme支持、全E完了を追加しない。

## 7. 追加zero-planeのfirst jetを直接確認して全六磁気へ強化する

前節までのgeneric開だけでwhole親を除くには十分だが、special planeを直接確認すると、今回は陰性を全源へ強化できる。これは一般点からの外挿ではない。

全域actual shearの時点で止め、Jを反転せず `R=Ih²+Jhy` とする。nodeのh=1 common base(0,0)はw=-Iのline、(1,0)はI=0の時だけzero-planeである。この追加plane上でもPxz=2x(x-1)のx微分は+2なので全点simple。I,Jがともに零でも同じである。

cuspでは(-1,0)がw=-Iのline、(0,±1)はI±J=0の時だけzero-planeとなる。それぞれPxz=2x(x+1)のx微分が+2なので、追加二planeの全点もsimpleである。同時I=J=0による両planeも捨てていない。

h零,x非零側のPfaffianはRに依らず、前節の非零rank四論証が全I,Jで通る。h=x=0,y=1のinfinity zero lineもz=0のままで、そのfirst jetは

```
node: J*h+2w0*x+3z,
cusp: J*h+2w0*x-3z.
```

±3が非零だから、ここも全I,Jでsimple。vertex lineの全quadratic jet(1),(2)は、Rが二次であることだけを使い、I,Jに一切依らない。

従って**全I,J、元全六磁気係数に対して、全order二zeroはvertex lineだけで、その全jetがC22の外である。第一二の全actual sourceは点ごとにwhole T36の外**。J零境界は別の閉性の推定でなく、この直接first-jet証明と全shear逆によって扱った。T38のspecial zero-planeに対する結論を本強化へ混ぜず、陽性T37親はなお未決である。
