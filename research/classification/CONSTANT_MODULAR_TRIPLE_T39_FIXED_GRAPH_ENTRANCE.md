# constant TRIPLE の一 actual orbit：固定 T39 全graph認識の入口

2026-11-08。実行前の有限契約。新 arithmetic はまだ無い。

## 1. 全sourceから同じcanonicalへ

`CONSTANT_MODULAR_TRIPLE_LINE_CUBIC_WHOLE_FINITE_REDUCTION.md` §§4–5の全四parameter源を使う。
物理座標(t,x,y,r)、f0=x³+3xy−3y/2+r として

```
f=f0+p*x+d,
Yx=x/3+y−l,
Yy=−xy+l*x+y/6−r/3+v,
Yr=−3y²+xr+L*x+(3l−p−5/4)*y+r/2+w,
v=2p/9+l/2−L/3, w=d+l*p+p/3+3*l/4−L/2.
```

全p,d,l,L、元の全excluded-zero補集合を残したnecessary/sufficient源である。
既accepted actualnormalizerはまず

```
r'=r+p*x+d, t'=t+(p/3)*y
```

で p,dを零にし、l'=l+p/3、L'=L−d−p/6 を得る。
次のf0-preservingmapを使う（新座標を大文字で示す）：

```
X=x+A, Y=y−A*x,
R=r−3*A*x/2−3*A*y−A³,
T=t+β*Y+γ*R+α*X,
β=−A²−A/2, γ=−A/3.
```

全四sourceの全六括弧を保ち、l→l−2A/3、L→L−2A³−2A²−A−3α。
A=3l/2、α=(L−2A³−2A²−A)/3 で全四を零にする。
逆はX,Y,R,Tの順で元x,y,r,tを解くtriangular polynomial、全det一、tensor一。
parameter unitで割らず、全四sourceそのものがONEactualorbitとなる。
この全normalizerの両逆、full150source復元は元ambient reportでも検査済みである。
本契約は別sourcepointを探していない。

## 2. 同じcanonical tensorの固定translationとrational Darboux frame

r_old=r_new+1、他三座標を不変とするactual affine変更を一度だけ施す。
逆 r_new=r_old−1、det一、tensor一、Euler零。
新方程式は f=f0+1。
canonicalの全六括弧から定数項は、順序(t,x,y,r)で

```
K=[[0,0,−1/3,1/2],
   [0,0,1,3/2],
   [1/3,−1,0,0],
   [−1/2,−3/2,0,0]].
```

Pf(K)=−(−1/3)*(3/2)+(1/2)*1=1、detK=1。
従ってこの**同じ**actualtensorにfaithful S61 adapterを適用出来る。

旧reusable `darboux` の決定的な最初のpair選択(0,2)に従えば、
Gram–Schmidt rowsは次のSとなることを手で確認した：

```
S=[[1,0,0,0],
   [0,0,−3,0],
   [3,1,0,0],
   [0,0,1/2,1/3]].
```

正確にSKS^t=J（J12=J34=1）、detS=1。
Darboux coordinates(a,b,c,d)=(t,−3y,3t+x,y/2+r/3) の逆は

```
t=a, y=−b/3, x=c−3a, r=3d+b/2.
```

元物理r_oldに戻す時は最後に+1する。
全source parameterで変わるframeでも、generic pointを次々変える操作でもない。

## 3. 一つの既全graphを読む最弱含意

変更しない旧caller `code/classification_modular_source_graph.py` が
全150Euler quotient復元、全134originalS61equations、全61座標、
exact Darboux/perspectiveと全45Q51 graph residualを記録する。
Q51は既証明のwhole closed T39/S61 graphそのもの。
全45が零なら、このcanonicalはactual PGL逆を通じてT39に入る。
§1の全四sourceは同じactual orbitなので全てT39に入り、
元fullfinite reductionの他known branchesと合わせ全highest-tripleがknown closed unionへ入る。
scheme smoothness、local uniqueness、ambient tangentの別rankは一切不要である。

旧callerはQ52も固定で読むが、それは旧入口を変更しないconsistency controlであり、
本研究の目標は既wholeT40排除後に残るT39だけ。
旧reportの歴史的field `point_result_requires_smoothness_to_transport_to_whole_source`
は今回の全actualorbit proofには不要であり、このfieldを根拠にしない。

一回だけ、wall60秒/RSSbudget1536MB。
陰性又はresource失敗ならその全exact residualを保存してSTOP。
origin変更、別point、period/profile/minor/jet、graphの再試行へ自動継続しない。
全陽性には別作者のfaithfulframe/全45恒等式と全sourcebridge査読が別に必要。
