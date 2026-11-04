# NEW25 の全PGL moduli：二つの不変量と二次の体の挟み込み

2026-11-04。実際のfull Poisson族P(u,1,v)のPGL5 orbit-invariant fieldについて、
既知の同時符号変換を除いた二parameter体と、その明示的な二次部分体との間に挟む。
これにより**どちらの端点かを未決のままでも、generic moduli体はQQおよびC上で
純有理・超越次数2**と分かる。未決の追加同値はsigma:(u,v)→(u,−32u−v)一つへ絞られる。
sigmaが実際のfull PGL同値か否か、全特殊parameter、有限stabilizerは分類しない。

## 1. firstjet座標から得るglobal不変量A=u²

既受理のactual B15 invariant fieldの指数1と有限逆recipeにより、primitive βから
shape座標rを有理関数として回収できる。これはsource parameterをそのまま不変量と
呼ぶ議論ではなく、target上の実際のrational functionである。

元のfull source六変数をt=2 sectionへ戻すと、r=u²/1600である。uは元source derivation
δの第一積分なので、この式はt=2だけでなく全六sourceで成立する。さらにn1正規化と
modular Gaによる評価点y=0への移動は全tensorの有限同値で、u=n0は評価点やPoisson
scalar cに依存しない。従って一般のrank4評価点pについて

\[
1600r\bigl(\beta(P(u,1,v),p)\bigr)=u^2.
\]

任意の固定pに対する値の有理symplectic正常化とperspective除去の間の違いはSp4作用で
あり、primitive rはそれを消す。scalar cを変えても右辺は同じである。よってこの関数は
projective Poisson tensorに定義される。

full成分上の一般tensorΠに対し固定点e0でこの正常化を行いA(Π)=1600r(β(Π,e0))とする。
Π=g*P(u,1,v)と書けるdense sourceでは、これは元Pをg⁻¹e0で評価した値に等しい。
前式は評価点に依存しないためA(Π)=u²となる。さらに任意のPGL変換を合成しても同じである。
dense full source上のこの有理恒等式により、Aは**全成分のrational PGL5 invariant**である。
全評価点・境界でregularな関数へ延長したという主張はしない。

## 2. intrinsicなH、conicと投影されたfoliation

既受理のfull anticanonical divisorは3H+Qで、一般にHは唯一のmultiplicity3 hyperplane、
Qは滑らかなquadricである。Q|Hはrank3のconeで、その唯一のvertexは
[x:y:z:w]=[0:0:1:0]。従ってH、Q、このvertexはfull projective tensorからintrinsicに決まる。

Hのvertexからの線型投影はP2(x,y,w)を与える。元の全二次tensorのh=0部分は

\[
P|_H=\partial_z\wedge V,
\]
\[
V=64w(x+vy)\partial_x+
 \left(-\frac{uv}{2}y^2-uxy-w^2+\frac{x^2}{64}-64wy\right)\partial_y
 +(x^2+vxy)\partial_w.
\]

これは保存済みfull Pの全degree2係数から直接照合する。Vはzに依存せず、Euler方向を
除いたline fieldとしてP2へ降りる。従ってこのfoliationもintrinsicであり、full PGL同値は
そのprojective同値を誘導する。tensor scalarはline fieldを変えない。

Q|Hから得るconicをscalar倍して

\[
q_2=-x^2+64w^2+64uxy+32uvy^2
\]

とする。正確にq2=72Q|Hである。行列のdeterminantは−2048u(32u+v)なので一般に滑らか。

## 3. 二つのinvariant lineのintrinsicな一意性

l+=x+8w、l−=x−8w、F=x+vyと置く。全係数で

\[
V(l_+)=8F l_+,\qquad V(l_-)=-8F l_-,\qquad V(q_2)=-uF q_2
\]

が成り立つ。従ってfoliationはclosed rational one-form

\[
\omega=d\log\frac{q_2}{l_+l_-}
          +\frac{u}{16}d\log\frac{l_+}{l_-}
\]

で定義できる。三divisor上のresidueは1,u/16−1,−u/16−1であり、degreeを掛けた和は
2+(u/16−1)+(−u/16−1)=0だからprojectiveに定義される。

**genericな幾何的parameter**で、これ以外のinvariant line Lがあると仮定する。
Lがpolar divisorの成分でないので、その正規化P1へのωの引戻しは零である。各点でのresidueは

\[
\operatorname{ord}_p(q_2/(l_+l_-))
 +\frac{u}{16}\operatorname{ord}_p(l_+/l_-)=0.
\]

二つのordは整数であり、generic uは有理定数ではないので両方零である。従ってL上で
l+/l−もq2/(l+l−)も定数になる。最初の条件からLはO=(l+=l−=0)=[0:1:0]を通る。
しかしq2(O)=32uv≠0であり、l+l−はOで零だから二番目の比はOにpoleを持つ。矛盾。
従って幾何的generic parameterでinvariant lineは**ちょうどl±の二本**である。
特殊な有理uに対する共鳴・追加代数葉の不存在を主張していない。

この一意性により、full PGL同値がunordered pair {l+,l−}を保つことが分かる。
特殊値を除く条件を無根拠に一つのZariski openと呼ぶ必要はない。幾何的generic fibreで
一意に決まるunordered pairなので、対応するrational constructionは函数体上で下降する。

## 4. dual conicによる第二global invariant

q2の対称行列Mとline covectorell±を用いて

\[
\boxed{J=256\,
 \frac{(\ell_+^{\mathsf T}M^{-1}\ell_-)^2}
 {(\ell_+^{\mathsf T}M^{-1}\ell_+)
  (\ell_-^{\mathsf T}M^{-1}\ell_-)}}
\]

とする。Mのscalar倍、各lineのscalar倍、二本の交換の全てに不変で、projective座標を
変えても同じである。unordered pairの上で対称なので、lineを個別に有理選択できなくても
Galois下降し、full component上のQQ-rational PGL invariantとなる。

実際の族では

\[
\ell_\pm^{\mathsf T}M^{-1}\ell_\pm=\frac{32u}{v+32u},\qquad
\ell_+^{\mathsf T}M^{-1}\ell_-=-\frac{2(v+16u)}{v+32u},
\]

なので

\[
\boxed{J=(v/u+16)^2}.
\]

これらは3×3 inverseによる全有理式の恒等式であり、点の評価ではない。

## 5. 二次の挟み込みとmoduliの有理性

k=QQまたはCとする。full PGL invariant fieldをEfullと書く。PGL×二parameter族から
actual成分へのdominanceにより、そのpullbackはk(u,v)へ単射である。既知の**full tensor**
有限同値(u,v)→(−u,−v)を保つので

\[
E_{\rm full}\subset k(u,v)^{\pm}=k(A,b),\quad A=u^2,\quad b=uv.
\]

前節までの二つのglobal invariantを合わせると

\[
\boxed{k(A,J)\subset E_{\rm full}\subset k(A,b)},\qquad
J=(b/A+16)^2.
\]

A,Jは独立であり、bの多項式は

\[
b^2+32Ab+(256-J)A^2=0.
\]

discriminantは4A²Jで、k(A,J)のJ-adic valuationが奇数なので平方ではない。
従って右端/左端の体次数は正確に2。この拡大の中間体は二端点しかなく、

\[
\boxed{E_{\rm full}=k(A,b)\quad\text{または}\quad k(A,J)}.
\]

**どちらも純有理2変数体**である。従ってendpoint未決でもactual generic full PGL moduliは
QQおよびC上でrational2である。この論証はX25自体の有理性を用いていない。

残るfield involutionはb→−32A−bで、uの符号を既知tauで調整すれば
sigma:(u,v)→(u,−32u−v)である。generic parameter同定は既知tauだけか、それにsigmaも
加わるかの二択へ絞られる。sigmaのfull tensor有限同値行列を得た、あるいは不存在を
証明したわけではない。次に必要な有限判定はこの一つのsigmaをfull tensorで比較すること。
同一parameter内の追加有限stabilizerの分類とは別問題である。

## 保存・再現と範囲

`code/q_new_family_global_two_invariant_preflight_20261104_v63.py`と同名report directoryに
全H-leading tensor、projected V、conic matrix、cofactor、dual pairings、二次式を保存する。
point、新しい高次数不変式、GB、sigmaの有限行列探索は使わない。

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_global_two_invariant_preflight_20261104_v63.py \
  --output work/global-two-invariant-replay
```

上限60秒/2GB。概念的global invariance・generic line uniqueness・体の挟み込みは
専用の独立監査へ渡す。全boundary、全特殊u、全finite PGL stabilizerの分類へ拡張しない。
