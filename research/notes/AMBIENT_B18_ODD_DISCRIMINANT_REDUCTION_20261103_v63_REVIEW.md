# B18の実odd idealの判別式上での像：独立監査

2026-11-03。CA/CB共通底とは次元の異なるNの正則底を別に確認する。独立監査済みN固定旗基底と独自の150係数Lie作用から全source Jacobianを作り直した。27実target行のminor−1152/5と、核に沿う全有限source曲線を係数表現の冪零指数関数で照合した。全150二次補正は保存値と一致し、実target row27は−3となる。従って28個の実target関数が長さ2のparameter fiberを作り、正則sourceのCM性によりN正規化環は自身の正則底上free rank2である。

三者それぞれについて、traceを引いた正規化環は `B=R[t]/(t²−ν)`。Bの正則性からνはRの正則parameterである。実像Aはactual target関数の係数を使って `R+I*t` と書ける。Rをtarget座標として含むことが重要で、形式的なnormalizerのtensorではない。

Rの座標をνと残りに取り、source座標をtと残りに取ると、target関数のt微分はtを法としてそのodd係数となる。従って全最大minor idealは正確に `(t,I)B`。その商は `R/(ν,I)` である。tは非零線形部分を持つためcritical idealの最小生成元の一本であり、除くとminimal generatorsはCA/CBで11−1=10、Nで9−1=8となる。

正則ambient `R/(ν)` における商のminimal Cohen presentationの一意性から、Iの判別式上への**像**は、三者別々の形式座標変換により六つの線形座標とCA/CBでは `(wu,wv,uv,u³)`、Nでは `(wv,v²)` で生成される。各小monomial集合は互いを割らず、六線形と合わせ最小生成数10/10/8を与える。実Iの最小生成数にも同じ下界がある。

ここで像 `(I+(ν))/(ν)` を `I/νI` と同一視しない。CAとCBの座標変換が共通Rの中で整合するとも主張しない。全odd idealのν方向へのlift、実CA/CB交叉の全完成環・nilradical・Ass・depthは未解決である。
