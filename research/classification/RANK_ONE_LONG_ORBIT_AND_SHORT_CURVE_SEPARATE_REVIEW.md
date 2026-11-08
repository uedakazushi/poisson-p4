# rank一longのgeneric軌道と全short実曲線：独立純証明査読

2026-11-08、quartic_review、PRIMARY E_multiple_hyperplane。**PASS**。
既新25成分の独立fullambient receiptに続く、全parameterの純証明である。
新計算job／point／operator／minor／profileは作らない。

## 全long非定数gの二つのprojective stabilizer場

固定modular場X=∂x+x∂yに加え、元全六括弧へdg/gをcolumn-sharpで縮約すると

```
Yg=d∂x+[dx+Kw gz−Kz gw]∂y+gw∂z−gz∂w,
d=divU.
```

U(g)=dgとH=gKから、各componentを直接手で確認した。
特にy-componentではr gw gzが相殺し、(gK)w gz−(gK)z gwは
g(Kw gz−Kz gw)となる。
gはaffine、Kは二次なので全Ygはaffine linearで、g零にも正則。
g非零上で閉一形式dloggのHamiltonian場なのでL_Yg P=0、恒等式は全域へ延びる。

gの線形部分が非零なので(gw,−gz)≠0。
Xはq成分零、Ygは非零constant q成分を持ち、両者は全sourceで独立。
斉次化したh-rowは両方零であり、非零線形結合がEuler場になる可能性もない。
char零のalgebraic stabilizerのLie次元は少なくとも二、全sourceのPGL軌道次元は22以下。
既**同一点**の独立M47／Euler24／scalar一によるorbit22下界と非零minorの開性から、
25次元成分のgeneric orbitは22、generic Rosenlicht quotientの次元は三である。
特別点全ての軌道次元又は全軌道分類は主張しない。

generic実D=h²LQのQはsmooth rank五、Hは唯一double factor。
g=z表示のp=[0:0:1:0:0]ではQ(p)=0,dQ(p)=dhなのでHはQの接hyperplane、Lもpを通る。
L≠Hは非定数gから従い、K=z²+w²の例でQ|Lのrankは四なので、
generic LがQに非接という非空開条件も確認できる。

## 全short sourceの元六括弧による実曲線

旧全short入口はq=(s,w)、r=y−xzで、全affine U、g、全F次数二／H次数三を保ち、
U(g)=dg、(U+1)H=Fgを課す。
lambdaに対しs=x−lambda z、r=y−xz+lambda z²/2としたproducerの全六式を手で再導出した。
補助自己同型(t,r,s,w)=(z,y−xz+lambda z²/2,x−lambda z,w)の
逆はx=s+lambda t,z=t,y=r+st+lambda t²/2。
Jacobianは全lambdaで**−1**であり、旧一般入口の「一」という記述を符号まで仮定していない。

∂t∧Y+Jac(rg+H)をこの逆でpushすると、
∂t→lambda∂x+x∂y+∂z、∂r→∂y、∂s→∂x+z∂y。
これより全六括弧がproducer (1)と厳密に一致する。
divY=−1、Y(rg+H)=−(rg+H)から全Jacobiが従う。
元座標のdegreeは全lambdaで二以下、lambda係数も全てpolynomial。
lambda零で元short全六式に戻り、Pf=rg+Hは全lambdaで非零なので曲線の射影点は消えない。
非線形自己同型はここでPoisson恒等式の証明にだけ使い、PGL帰属の代用にしない。

## Punctured族の真のaffine frameと全parameter逆

lambda非零で(tL,Y,q,w)=(x/lambda,y/lambda,x−lambda z,w)を取る。
det=−1/lambda、逆x=lambda tL,y=lambda Y,z=tL−q/lambda。
この**affine linear**変更はmodular場を∂tL+tL∂Yへ送る。
r=lambda rL+q²/(2lambda)を用いて元六括弧を戻すと、

```
UL=u∂q+(v+g/lambda)∂w, dL=d+gw/lambda,
HL=H/lambda+q²g/(2lambda²),
FL=F/lambda+[qu+(1+d)q²/2+Hw]/lambda²+q²gw/(2lambda³).
```

全degree制限を保持し、UL(g)=dL g、(UL+1)HL=FL gも手展開で確認した。
たとえば後者は旧(U+1)H=Fgに、gHw/lambda²、
g[qu+(1+d)q²/2]/lambda²、q²g gw/(2lambda³)を加えた等式となる。
逆の反転はlambdaだけで、零にはこのpunctured frameを代入せず、分母なし元曲線を使う。

非定数gでは旧全入口でg=sへ正規化すると
(1+beta∂w)H(0,w)=0。有限次数多項式上でこの作用素は可逆なので、**全sourceでH=gK**。
従ってHL=g[K/lambda+q²/(2lambda²)]となり、全short非定数gは新long25閉成分へ入る。
このdivisibilityを任意short源に追加仮定してはいない。

非零定数gではd=dL=0、HLの三次部はH3/lambda。
H3零の全short源は既long定数g／H3零NEW25へ入る。
三次共鳴はULの線形部分が元U1のままなので全共鳴／零係数を保持して対応するlong残余へ戻るだけで、
未解決のlong共鳴を閉じたとはしない。
g零ではUL=U、HL=H/lambda、全Fも許された次数二のまま対応long残余へ送るだけである。

各閉PGL安定族Cへpunctured曲線が入れば、元座標の実polynomial曲線と閉性により中心もCへ入る。
全sourceを一つの点、固定jet又はmodular軌道の退化だけで代用した証明ではない。
新成分増分は零。long g零／定数g三次共鳴と全Eは明示的な別未完課題として保持する。
