# HF-A全源のT36実同定

2026-11-08、PRIMARY E。計算前に保存した `HF_A_T36_EIGENFRAME_AND_VERTEX_REDUCTION.md` の全frameを、同じ整数e=1点で完全検証した。別作者レビューは別receiptに委ねる。

Q(√6)上で、b非依存Ore sourceをraw十括弧から再構成し、modular場と四固有座標を確認した。vertexの全20×8線形系は唯一解ell=−3c+2v、f=144a+60uを持つ。

```
x0=b−3c+2v
x1=−72a+5u−15√6(9c+v)
x2=−72a+5u+15√6(9c+v)
x3=u+√6(v−3c)
x4=u−√6(v−3c).
```

実PGL行列式は−311040であり、左右逆行列を保存した。完成T36の全correctionを含む値は

```
(p0,t0,t1,t2,t3)=(300√6,−1800,−1800,−300√6,300√6).
```

tensor scalar1で、全raw恒等式は

```
M_*Π_HF-A(e=1)=Π_T36+E_X∧M_*Rtotal,
Rtotal=(1500c+150v,144a+60u,40a,360c−240v,144a−40u).
```

元全150二次係数の残差、全10のE∧差、およびsigned Euler-Pfaffian差は全て零である。五次式の符号も

```
D_T36(X)=−311040 D_HF-A(e=1)(M^{-1}X)
```

として検査した。reportに全M/逆、全core/W、固有場、vertex線形系と唯一解、元/先の十括弧、明示Euler場、両五次式、およびinline実行sourceを保存した。gate PASS後の一回、1.077秒である。shared6022/helperは不変であり、新しい点・ambient rank・minor・profileを作っていない。

既保存の完全有限graphはanchorでsmooth unique24次元成分を与える。上の点は閉T36上にあるので、その成分はT36である。graphのe≠0はトーラスで一軌道、e=0も閉T36に入り、anchor23軌道閉包はT36に含まれる。HF-Aの全PGL×A10像閉包がこの23軌道閉包に等しいという前段の全源証明から、HF-Aの全点がT36に入る。generic higher parameterだけの帰属ではない。
