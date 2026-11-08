# H＋二既約quadric・joint SS rank一：全実衝突の独立査読

2026-11-08、quartic_review、PRIMARY E。**PASS**。
元の全括弧表、全source unitと零補集合、真のT16/T17全period、全actual frame正逆、
全六括弧のpole相殺と中心、全Pfaffianを別作者として手で再導出した。
新しい計算job、operator、profile、point又はjetはない。

## 真のnamed parentを直接戻す

T16 raw catalogueの座標(x0,…,x4)=(A,B,C,D,H)で、
L=(0,(rho1−rho2)B/2,(rho1+rho2)C/4,−(rho1+rho2)D/4,0)を引く。
Praw−E∧Lの全六括弧は
F=pAB+uH²、G=qCD+vH²、cross=k(A∂A−B∂B)∧(C∂C−D∂D)。
全periodへの逆は

```
(rho0,rho1,rho2,t0,t1)=(q+2k,−p−2k,p−2k,u,v).
```

ここでu,vはRAW constant係数である。
旧chord helperのF=p(AB+uH²)等に対するt0=pu,t1=qvという別normalizationと混同しない。

T17は未完成の単純edge和を仮定せず、完成JSONのt0t1H²/rho1項まで読む。
L=(0,−rho1 B,−(rho1−rho2)C/2,−(rho1+rho2)D/2,0)を引くと、
F=pAB+uH²、G=qCD+vF、cross=(A∂A−B∂B)∧(aC∂C+bD∂D)、a+b=pとなる。
全period逆は

```
(rho0,rho1,rho2,t0,t1)=(q−p+2a,p,p−2a,u,pv).
```

唯一のsource分母pは今回のpunctured族で非零である。
GのvFを独立な三係数へ緩めていない。
全六Jacobiはa+b=pによって零、Pfは両parentともFG。
named source認識と閉PGL不変性は既受理の全MPS source定理を固定入力とする。

## B1／B2：全自由係数を保った一次衝突

元reportのB1を直接読み、F=n0h²−2n5xy、alpha=−n8/(2n5)、beta=n9、
zeta=z+(n3/n5)hと確認した。Pzetawは**alpha F+beta w²**であってn8Fではない。
元unit n0,n5,n8を保持し、beta零ならactualDがhF²で相異二既約quadricではない。

T17 frame C=w−epsilon zeta、D=w+epsilon zeta、
p=−2n5,u=n0,q=−2epsilon beta,v=−2epsilon alpha、
a=−n5+epsilon n4,b=−n5−epsilon n4を代入する。
det=−2epsilon、逆zeta=(D−C)/(2epsilon),w=(C+D)/2。
全六inverse括弧はproducer note §3のpolynomial表と一致する。
特にPxy=F、Pxzeta=−x(n5zeta+n4w)、
Pxw=−x(n5w+epsilon²n4zeta)、Pzetaw=alpha F+beta(w²−epsilon²zeta²)。
y側の符号と全零境界n3=n4=0も同じ恒等式に含む。
中心は元の全六parameter源、actualPfはhF times最後の二次式である。

B2のn0零branchは元の全Pfにxy因子があり対象外。
n0非零では元n6,n9 unitを保持し、
w'=w+(n5/n6)h、z'=z−(n4/n6)w、d=n3−n4n5/n6を使う。
逆はw=w'−(n5/n6)h、z=z'+(n4/n6)wであり、全n3,n4,n5を保持する。
F=n0h²−n6xy、alpha=−n9/n6、beta=n8とすると元Pz'w'=alpha F+beta hw'。
beta零は同じdouble-F対象外である。

C=h+epsilon z',D=w'、p=−n6,u=n0,q=epsilon beta,v=epsilon alpha、
a=−epsilon d,b=p+epsilon dをT17へ入れる。
det=epsilon、逆z'=(C−H)/epsilon。
全六括弧は、Pxz'=−dx(h+epsilon z')、Pxw'=(p+epsilon d)xw'、
Pz'w'=alpha F+beta(h+epsilon z')w'と対応するy側、Pxy=F。
全てpolynomialで、元full表と中心で一致する。

## B3：三次poleと全branchを保持する

Q=n4h²+n6(w²−2hz)、Qepsilon=Q−epsilon²n6z²。
frame H=h,A=h+epsilon w+epsilon²z,B=h−epsilon w+epsilon²z,C=x,D=yの
detは−2epsilon³、逆はw=(A−B)/(2epsilon)、z=(A+B−2H)/(2epsilon²)。

n1非零branchの元unit n1,n6,n7を保持し、T17へ
p=2epsilon n6,u=−2epsilon n6−2epsilon³n4、
q=−n7,v=−n1/(2epsilon³n6)、
a=epsilon(n3+2n6),b=−epsilon n3を入れる。
厳密にF=−2epsilon³Qepsilon、G=−n7xy+(n1/n6)Qepsilonとなる。
parent crossの係数aはC=xに、bはD=yに作用することを保って計算すると、

```
Pxz=−(n3+2n6)wx,  Pxw=−(n3+2n6)(h+epsilon²z)x,
Pyz=n3wy,          Pyw=n3(h+epsilon²z)y,
Pxy=(n1/n6)Qepsilon−n7xy, Pzw=Qepsilon.
```

従って三次period poleを含め、全六括弧で正確に相殺する。
actualPfはhQepsilon[(n1/n6)Qepsilon−n7xy]。
中心で元の全表へ戻る。

n1零branchは別の全unit n0,n6,n7を保持し、T16を使う。
A=x,B=y,C=h+epsilon w+epsilon²z,D=h−epsilon w+epsilon²z、
p=−n7,u=n0,q=2epsilon n6,v=−2epsilon n6−2epsilon³n4,k=−epsilon n3。
G=−2epsilon³Qepsilon、Cartan C∂C−D∂Dのinverseは
(w/epsilon)∂z+(h/epsilon+epsilon z)∂w。
全六括弧は元表のPzwをQepsilon、Pxw/Pywのhをh+epsilon²zへ替えるだけ。
Pxy=n0h²−n7xy、actualPf=h(n0h²−n7xy)Qepsilonとなる。
ここでも全epsilonでpolynomial、n3又はn4をunitと仮定しない。

## Whole Eへの矢印と範囲

旧五supersourceの尽くしとJacobi枝は既独立認証済み。
A1及びA2非零ss枝は旧factor排除、A2零ss枝は既全plane被覆へ戻す。
B1/B2/B3で反転したsource係数は元unit n5/n6だけ、追加beta unitは
相異二既約quadricという対象条件から得て、その零補集合を明記した。
B3のn1零／非零の両方を別parentで覆っている。

各punctured族は真の全period named sourceのactual可逆PGL像で、
全target括弧がpolynomial、中心tensorは元rank四点として非零である。
従ってclosed T16/T17から中心の全source点も含まれる。
全元normalizer逆を戻し、指定lambda零・joint SS rank一の全H＋相異二既約quadric源は既知MPS有限閉和へ入る。
新成分増分は零。joint SS rank零や他factor型の完了をここだけからは主張しない。

本receiptは全六恒等式と正逆の手による別作者再導出であり、
新しい数値／CAS算術runを実行したとは表示しない。
