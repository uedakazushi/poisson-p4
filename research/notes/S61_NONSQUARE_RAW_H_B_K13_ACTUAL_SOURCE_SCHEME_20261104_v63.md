# actual B-rankbad schemeの全source・開条件ゲート

2026-11-04。対象はR=QQ[a,e,d,1/(adFBCJ)]上の**actual 7×3 block B**の最大minor ideal I=I_3(B)である。full Poisson cubic matrixやcoker Bの長さとは区別する。

元の22個の正規化minor生成元は合計2963項、(a,e,d)次数は高々(7,11,7)、異なるmonomialは283個である。e=t(a)、d=s(a)は固定済みのKa13=QQ[a]/g13の全係数を使う。新しいaの値・点・素数を選ばない。各monomialを一度だけ計算し、すべての乗法についてden_out*num_left*num_right−den_left*den_right*num_out=g13*quotientという整数多項式恒等式を保存する。その後の22個のsource和は、次数12以下のQQ[a]係数の線形相殺として検証する。巨大な未還元t^11*s^7は展開しない。

source openのFULLU=adFBCJは同じsource49項から評価する。受理済みNAMED_SATURATIONの七つの整数Bézout証明を再検証し、H4(d)=(d-r)Q_r(d)+H4(r)とH4(s)=0から(s-r)^−1=−Q_r(s)/H4(r)をr=0,rB,rJに適用する。先頭単元の逆元を掛けてFULLU(s)の逆元を得る。Gの既存整数Bézoutも再利用する。新しい逆元探索はしない。すべてKa13の全体で行い、既約性・体・étaleを仮定しない。

scheme同型の証明は両方向の写像である。既受理の全Q/L閉成分の除外はGをR/Iで単元にする。normのsource membershipによりH,h1,k1,m1はIに入り、source bridgeがe=tとg13=0、d-Euclidがd=sを与える。従ってKa13からR/Iへの全射alphaがある。逆に元の全22生成元が(a,t,s)で零で、FULLUに逆元があれば、R/IからKa13への単位的写像betaがある。beta*alphaは生成元aを固定するので恒等写像であり、alphaは同型である。これは長さだけからの推測ではない。

このゲートが成功すればactual Fitting bad algebraは正確にKa13、QQ長さ13、非空である。一方、reducedness・étale・既約性・幾何点数は別問題として残す。heavy blockのrankbadからfull Poisson cubic rankの低下、quartic239、generic exhaustionを結論しない。

一回60秒・2GB・32MiB・明示150000項の上限で実行し、各node後にcheckpointを保存する。再現: `.venv/bin/python code/s61_nonsquare_raw_H_B_K13_actual_source_scheme_20261104_v63.py --output work/K13-actual-source-fresh`。新しいsource gateの独立レビューは別途必要である。

実行は1.264秒で完了した。全22生成元は零で、全49項sourceから得たFULLUの値と逆元の積は厳密に1だった。monomial nodeは312、保存した明示項は6932、monomial係数は最大2091bit、商環整数演算は417回。FULLU値は13項・最大分子1610bit/分母1557bit、逆元は13項・1950bit/1918bitである。全artifactは約6.55MBで、巨大未還元式や新しい逆元探索を必要としなかった。

従って受理済みsource gateの下で**actual I_3(B)の商環はKa13と同型、QQ長さ13で非空**である。これでこのBブロックがFULLU全点でrank3という目標は否定されたが、他のpure/genuine方向を持つfull cubic matrixの階数について何も否定していない。次段はこの全13次閉代数上で、元のfull cubic sourceを評価する問題である。
