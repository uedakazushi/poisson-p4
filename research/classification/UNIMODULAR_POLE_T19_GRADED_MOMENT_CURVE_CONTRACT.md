# 残h5実軌道：T19の自然weight moment曲線一つの有限恒等式契約

2026-11-08、PRIMARY E_multiple_hyperplane。準備のみ、未実行。
新point／ambient profile／M1／local jetを計算しない。固定した真T19全sourceのactual curveだけを判定する。

## 純構造から固定したframe

P*=R+B、R=Jac(C,D),C=xy+z,D=xw+y+z²/2、
B=Z1∧Z2=(0,0,0,−1,−y,xy+z)。R,BはPoisson、互換、rank二、Pf(aR+bB)=ab。
Z1=x∂z−z∂w−∂y,Z2=y∂w+∂zは可換affine symmetry。
Y=x∂x−2y∂y−z∂z−3w∂wに対し[Y,Z1]=2Z1,[Y,Z2]=Z2、R/Bのtensor weightsは2/3。

真T19の二つの対角symmetryのcharactersはsource順(a,b,c,d,H)で
S=(-1,1,2,0,0),T=(0,1,1,1,0)。
P*のnil二場が生成するlength五regular algebraではcyclic basis
(1,S,T,S²,ST)に(w,−z,y,−x,−h)が対応する。
その自然weights S=e²,T=eとsource coefficients(−1,2,1,1,1)を採ると、固定actual mapは

```
h=−2e³(b+c), x=e⁴(a−2b−4c), y=e(2b+c+d),
z=−e²(a+2b+2c), w=−a+2b+c+d+H.
```

det=4e10、全eunitでinverseは

```
a=h/e³−z/e²,
b=x/(2e⁴)−3h/(2e³)+z/(2e²),
c=h/e³−x/(2e⁴)−z/(2e²),
d=y/e+2h/e³−x/(2e⁴)−z/(2e²),
H=w−y/e−z/e²+h/e³.
```

periodをr=R e³,v=r+V e²、edges t0=(r−v)/2,t1=r−vに固定する。
sourceは必要coupled correction 2t0t1 bH/(r−v)を含むtrue completed T19 JSONのみ。
source period openは−r(r−v)=RV e5、generic eではRVunitが必要十分。
q3=cd−b²,q4=ab+dHはleading −h²/(4e6),h²/(2e6)、Hはh/e³。
実D共変則より全target actualDの中心は−RV h5/2。
このPfaffianだけをtensor convergence又はP* membershipの証明には使わない。

## ONE finite literal test

全R,Vを未定元として、上記一つのactual frame／period sourceの全150canonical係数を戻す。
fullsource→H-Casimir helperの全Euler差、全inverse、元350Jacobi、全126actualPfと共変則を確認する。
negative e係数は全てR,Vに線形なので一回の全係数線形系で尽くす。
中心の全150係数をaR+bBという固定二次元空間と照合し、同じ全線形系で解く。
RVab unitの解がある場合だけ、source generic eのT19帰属と中心のactual normalizerでP*をT19へ送れる。
中心aR+bBのwholeunit正逆は、まずtensor1/a、λ=b/a、次に

```
(x',y',z',w')=(λx,λ^-2 y,λ^-1 z,λ^-3 w),tensor λ²
```

でP*へ移す。逆はλを反転、tensor a/λ²。
sourceはR=1へ全体tensor scaleを除いてもよいが、reportは全R,V線形系を保持する。

pole又は中心の全係数照合がunit domainで失敗したら、この固定curve classを非閉鎖として停止する。
追加e係数、cyclic vector、frame parameter、別point／profileを追加しない。
T19全体の非所属又は他の親の非所属へ拡張しない。
全60秒／1GiBでfailureもexclusive新reportへ保存する。
