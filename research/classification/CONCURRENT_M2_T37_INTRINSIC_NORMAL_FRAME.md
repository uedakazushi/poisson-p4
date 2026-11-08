# concurrent m二smooth24親：零secantとnormal三座標のT37 frame

2026-11-08、PRIMARY E。保存した完全graphの整数e=1だけを用いる、計算前の純frame導出。literal fit成功は後続の全150/Pf receiptへ区別して保存する。新しいpoint/profile/operator/jetを加えない。

## 零直線と全二特殊点

整数e一tensorはactual zero line

```
ell: a=0, u=−96b, v=−256c
```

を持つ。その全raw括弧はE∧Rellと一致し、

```
Rell=(−24576c,24768c,−54b,−10340352c,20736b)
```

である。例えばPab=24576bc、Pac=24576c²、Pau=−2359296bc、Pav=−6291456c²、Pbc=−54b²−24768c²で、残五つも同じEuler wedgeで一致する。

c一chart、normals(A=a,U=u+96b,V=v+256)におけるfirstjetの三行πab,πbu,πbvは

```
[[25920−162b²,                 −128,                0],
 [3773952+151632b², 25920−162b²,                       0],
 [202176b,                    −432b,162b²+33984]].
```

detは4251528(b²+1888/9)³、二特殊点ではrank一となる。θ²=−118を用いればb=±4θ/3である。無限遠b一,c零では同じ三normalの行をπac,πcu,πcvに変えるとmatrixは[[162,0,0],[-151632,162,0],[0,0,-162]]であり、追加特殊点は無い。

## Canonical strict liftと三次元normal block

全projective quadratic Poissonの一意divergencefree liftをPbar=P−E∧δP/5とする。ここでδP_j=Σi ∂i P_ij。これは既受理canonical liftであり、新しいglobal stabilizerの延長を仮定していない。

normal coordinatesを

```
A=a, U=u+96b, V=v+256c,
L=U+594A, eta=U−468A,
xi=L−2theta V, zeta=L+2theta V
```

と取る。Pbarに対するnormal三座標の全括弧は、κ=432thetaとして

```
{xi,eta}=κ(−6xi eta/5),
{xi,zeta}=κ(−2xi zeta/5−2eta²),
{eta,zeta}=κ(−6eta zeta/5)
```

となる。導出にはℓ上のEuler fieldを全体の一次場へ延ばして一度引き、閉normal三変数blockを取る。そのcanonical blockのmodular場は

```
(768V/5, 359424V/5, −1026432A/5−1728U/5)
```

である。etaの固有値零、xi/zetaの固有値は±3456theta/5。全normal bracketを比較すると上式のκと二次eta²係数も一致する。固有値だけから全tensor同定を推論しない。

## 残るframeは一回のlinear十係数問題

既canonical index二T37 modelの座標を(xi,eta,zeta,d,e)とする。その三normal括弧は上記κ一のもの、残七つは

```
{xi,d}=2xi d/5+4eta zeta,
{xi,e}=6xi e/5−2eta d−2zeta²,
{eta,d}=3eta d/5+3zeta²,
{eta,e}=−3eta e/5−3xi²,
{zeta,d}=−6zeta d/5+2xi²+2eta e,
{zeta,e}=−2zeta e/5−4xi eta,
{d,e}=−de/5+2xi zeta+3eta².
```

d,eを元五変数の全一次式（十未知係数）とする。最初の六mixed equationsはd/e両方を保持しても全て未知係数にlinearである。全六×十五係数を一回解き、唯一invertible frameが得られる場合に限り、最後のde式、左右逆、全150raw係数と全signed Euler-Pfaffianを検査する。失敗時はexact残差を保存し、別frame pattern又は点へ反復しない。

さらに完成済みT37全correctionのall-one源からindex二modelへの実frame

```
(xi,eta,zeta,d,e)=(2x3,x0,−2x4,−12x2,12x1)
```

と合成し、最終的には元のcompleted T37と保存整数e一tensorの直接全150/Euler/Pf恒等式を証明する。modelの名前やdivisorの同一性だけには依存しない。
