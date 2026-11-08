# Ga9 linear-w固定D全源の独立査読

2026-11-08、integrator `/root`。新算術operator/点/minor/jetを実行せず、producerの全源定理を手で検証した。親への包含は未証明。

D=-4h(4h³w-8h²y²+y⁴)、U=(4h²-y²)∂y+(16hy-4wy)∂w、V=-hy∂y+(-4hw+4y²)∂w。base(h,y,w)、fiber(x,z)。
[U,V]_y=-4h³-hy²、[U,V]_w=-16h²y-4hwy+8y³を独立微分し、[U,V]=-hU+2yVと照合した。det(Ebase,U,V)=h*((4h²-y²)(-4hw+4y²)+hy(16hy-4wy))=D。U(D)=-4yD,V(D)=-4hD。FはC[h,y]上primitiveなw一次で既約、gcd(F,h)=1、D被約。

A=U+4y Ebase/5,B=V+4h Ebase/5はDの二次syzygy。det(Ebase,A,B)=Dで、被約Dの一般smooth点でlog derivationのCramer分子は各primeを割る。従ってSaito basisとなり、Euler分解からAR(D)=C[h,y,w]A⊕C[h,y,w]B、AR0=AR1=0。全Dを仮定する定理で、既約quinticに限定された旧分類結果は使わない。

canonical lift Π#dD=0におけるbase-base skew blockのベクトルは∇Dの有理倍。∇Dはprimitive（squarefreeとEuler恒等式）ゆえ有理倍はpolynomial、degree2対degree4で零。mixed各列をfiber次数で展開すると、次数1/2部はAR1/AR0ゆえ零、次数0部は定数行列CによるA,Bだけ。detC=0ならbivectorは一つのfiber方向でくくれてsquare零だからrank四ではdetCunit。actualfiberGL2で全Cを戻せる。

Euler gauge E∧((4y/5)∂x+(4h/5)∂z)でmixedはU,V、h行零。このliftでprojectiveJacobiのhを含む四添字がh倍strictJacobiなのでstrictも零。三添字式は[U,V]=Kx U+Kz V、U,Vはfraction field上独立だからKx=-h,Kz=2y。従ってK=-hx+2yz+R2(h,y,w)、全R6を保つ。

actual x'=x+l,z'=z+m,l=ah+by+cw,m=dh+ey+fw、det1、逆shift符号反転。全有限pushは
ΔR=h*l-2y*m+U(m)-V(l)=(a+4e)h²+(2b-2d+16f)hy+5chw-(4c+3e)y²-6fyw。
各monomialを直接展開して確認した。所望五係数にc=rhw/5,f=-ryw/6,e=-(ryy+4c)/3,a=rhh-4e,b=0,d=(16f-rhy)/2が全域逆。分母は非零定数だけで、零磁気も保持。kernel l=by,m=bhはcanonical Ga field。
fiber同時κscaleとtensor1/κはmixed/線形Kを保ちR→κR。従ってRw²の係数unitは4へ、零はR0へ、全固定D源は正確二actual軌道。

D_w=-16h⁴,D_h|h0=-4y⁴なのでreduced singular supportは唯一S=(h=y=0)。上の代表をSに制限するとR_w²*w²∂x∧∂z。wunitの射影chartでrank二/零を確認出来る。ゆえにactualD=PGL(Dtarget)かつS一般rank二がGa9の必要十分認識条件となる。

この条件に親からのactual正則非零tensor曲線が一つ加わればwholeGa9軌道/全Runit閉包はそのknown親に入る。D-only極限、genericnormal根の極限、Srank零、旧localnormal支持だけではその条件を補えない。T36/T37等の親と全Eは依然未解決。
