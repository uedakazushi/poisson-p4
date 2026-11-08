# generic triple・R2L・d(c/a)非零の全排除

2026-11-08、producer純証明。
constant-root obstruction及びfullweight-one identityは `PLANE_M3_CONSTANT_ROOT_OBSTRUCTION.md`。
以下の全係数は同補題author plane_collision_reviewが別に手計算し一致を確認した。
新CAS計算又はchosen higher sliceを使わない。

f=a²b、h=a^e b^(2−e)、e=0,1,2、G3=h z³。
全nonlog diagonal character pにd(c/a)が非零ならactual c shear/scaleにより

`p=(dlogf)/3+αdlog(a/b)+d(c/a)`

へ直せ、全offdiagonal rは未固定のまま残す。
q=2p−dlog(f/h)は

`q=(e−2/3+2α)dloga+(5/3−e−2α)dlogb+2d(c/a)`。

constant bad column vector Cについてq(C)零、C非零が必要。
全homogeneous polynomial coefficientsを比較するとCa=Cc零であり、
Cb非零からα=5/6−e/2が必要。
従ってCを非零scaleで∂bとしてよい。

raw P0^{i,z}=zX_i、P0^{i,w}=zZ_i+wY_iと書くと、全X,Yは正確に

`X=(-a,2b,−c+3αa), Y=(-2a,4b,−2c+a)`。

これはω=3fβ−dfのtracefree inverseで得られる。
Zは全rを含み、以下のuniversal coefficientには現れない。

normal K=z(k20z+κw)とすると、原projective Jacobiのw³項は

`Ebase∧[(DY)C+κC]−(2Y−X)∧C+P_{C(f)}=0`。

DY C=4C、P_{C(f)}=P_a²=2a∂b∧∂c。
ab coefficientは(κ+7)a、bc coefficientは−(κ+7)c+(4−3α)a。
従ってκ=−7、α=4/3が必要となる。
しかしe=0,1,2のα=5/6−e/2と一致しない。

したがってd(c/a)非零の全diagonal sourceは、全r/K/higher sourceを保ったまま空。
logarithmic性を仮定せず、reduced-R residue比較を流用していない。
この結論をd(b/a)非零又はlog-only sourceへ拡張しない。
