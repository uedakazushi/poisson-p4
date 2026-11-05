# 六因子base上の非空な二tail十分条件gap

2026-11-05。独立受理済みの事実は、Hが選択二分子の厳密共通GCDでsquarefree、
gcd(H,adFBCJ)=1、そして同じ六因子環でH∤D12である。
最後の事実はROOTの別方式Newton値peer `041d9234...` により確定した。
本sourceは新しい値を計算しない。

U=adFBCJ、D12=N/uとする。ただしuは反転因子の積であり、D12はR=QQ[a,e,d,1/U]の元。
複素幾何点でV(H)∩D(U*N)が空と仮定する。char0ではH-squarefreeがgeometricにも保たれるので、
NullstellensatzによりU*N∈sqrt((H))=(H)。gcd(H,U)=1をUFDで消すとH|N、従ってH|D12となる。
これは独立非除算証明に反する。よって

\[
V(H)\cap\mathrm{FULLU}\cap D(D_{12})\ne\varnothing.
\]

Hの少なくとも一つの成分でU*Nは非零なので、その成分の非空なopenは次元2を持つ。
Hの全成分がD12 openに会う、あるいはHが既約という主張ではない。
次元2は(a,e,d)という三次元parameter base上のdivisor部分の次元だけを指す。

選択二tailは(Q0,Q2)=H(A0,A2)だから、この場所では両tailが零。
したがってこれは現在の選択二tailから得たfull57十分条件が判定しない、実在するgapである。
実cubic matrixがrankを落とす場所やbad Fitting schemeだとは言わない。
別のminor/source行がfull57を回復する可能性を保持する。

次の固定候補は既存source original32 column3のT3である。genuine current15 minorに関する
`A_anchor*D15_3=-T3*D12` とGLOBAL橋渡しは、完全な六因子環で成立する。
同じH_axis四次商環でD12とT3を別々に非零と調べるだけでは、zero-divisorのため不十分。
結合積D12*T3の余りが非零であり独立peerを得た場合にだけ、同じH成分の次元2 openで
このgenuine minorとactual full57が非零になると結論する。
第三tailと結合積の値はまだ実行していない。全H/FULLUの被覆・quartic239・尽くしは未解決。
