# Actual NEW fiberの明示的Luroth生成元とJの二次写像

2026-11-04。同じbeta底K_B=k(B_NEW15)上のactual reduced Y16について、
未決だった明示的生成元を構成する。k=QQ,C。既受理のfull source次数8、tau下降後4は変えない。
以下のr_f,s_fはflow座標であり、旧sectionのshape変数r=u²/1600とは異なる。

## 1. full tensorを保つ四つのsource変換

A=256、r_f=A+v(384X+3072W−1)、s_f=A+v(384X−3072W+1)とする。
F5=k(u,K,Iv,Ic,Ib)は既受理のderivation定数体である。
generic full stabilizerの実行列Eを評価点へ適用すると、X'=-X−4/(3v)、
W'=1/1536−Wとなり、r_f,s_fは両方反転する。Eはscalar1、affine determinant1なので
c,Bは変わらない。modular Gaでp_y=0へ戻す操作はX,Wを変えない。

独立したr_fだけの反転は、既受理のfull sigma equivalenceのsign+1を使う。
target v'=-32u−vで、行列のX,W成分へ代入するとr_f'=-r_f、s_f'=s_fとなる。
full tensorのscalar規約からc'=-cv²/v'²、PfaffianとdeterminantからB'=Bv³/v'³である。
ここでもunique Gaでp_y=0へ戻せる。全五つの定数はこれらの変換で恒等的に保存される。
従って二変換は同じfull targetを保つsource変換である。

raw-pair matching後のgroup変数hは、固定したbeta sectionとactual betaから一意に決まる。
このため、単にSp4同値なbetaを得るだけでなく、full tensorを保つ上記変換の下でhも固定される。
full E/sigmaの全係数恒等式は既存の独立証明書を再利用する。

## 2. 四次deck群の固定体

既存の有理parameter tに対し、D=Kt²−K+1、N=Kt−K+1とすると

```
r_f=A(D−2N)/D,  s_f=A(D−2N)/(D−2tN).
```

両符号反転とr_f反転はそれぞれ

```
t -> (Kt−K+1)/(K(t−1)),   t -> (K−1)/(Kt)
```

となる。二つは異なる可換involutionで、genericにV4を作る。
これはnormal form parameter同定のV4(tau,sigma)とは異なるsource deck群である。
r_f²は固定され、tについての次数は4だから

```
F5(t)^V4 = F5(r_f²).
```

## 3. tau下降と明示生成元

tauはuを反転し、Kを1/Kへ、r_fとs_fを交換する。
conic関係から

```
q=(r_f²−A²)/(s_f²−A²)=1+(K−1)r_f²/A²
```

であり、tau(q)=1/qである。従って

```
w_fiber = u*(q−1)/(q+1)
        = u*(K−1)*r_f²/(2A²+(K−1)r_f²)
```

はtau固定である。逆式は

```
r_f²=2A²*w_fiber/((K−1)*(u−w_fiber)).
```

よって(F5(r_f²))^tau=E5(w_fiber)である。これはu²∈E5、F5=E5(u)と
上の逆式からも直接分かる。full sourceを保つ群の次数は4×2=8で、既受理のactual
full source次数8と一致するため、包含は等号となる。結論は同じ埋込みで

```
k(Y16)=k(B_NEW15)(w_fiber).
```

actual Y16の既存の有理性を、初めて具体的なfiber生成元で実現した。
sourceのz=u(t−2)/tからの写像は依然として正確に次数4である。

## 4. full PGL不変量Jはactual fiber上で次数2

J=(v/u+16)²はfull tensorの有理函数である。
R=r_f²、S=s_f²=A²KR/((K−1)R+A²)とすれば、第一積分Ivから

```
J=[Iv*S*(R−A²)+8u*(R+S)]²/(u²*R*S).
```

Rの上の逆式を代入するとw_fiberについて互いに素な次数2の有理式を得る。
専用証明書には約分した式を保存した。独立にtについての次数8とsource次数4のtowerからも
actual fiber上の次数は2となる。これはNEW定義方程式の次数やspectral次数10/20/40とは別である。

全てgeneric函数体の主張であり、特殊fiber、全域finite morphismやbranch locusは未決。
新点・高次不変式・source26・大きいmatrix計算は行わない。
