# B63 の socle と delta63 の lower-cache 安全性：独立 SOURCE 検証

独立検証は PASS した。報告 SHA256 は
`195299be1b82211f4c26bdafa48cba579327531a79e36904da1072ab75d88b2f`。
元105入力の既存 cache、受理済み全105勾配零のsource、今回のB63 schemaを固定した。
deltaやその他scalar DAGの出力、head29/30応答は評価していない。

B63は同じR0座標 `(z0,…,z6,r,ell)` のmonomial quotientで、基底は
plane次数≤4の15個、zi*plane次数≤2の42個、z0/z4/z5の二次積6個からなる。
独立 reviewerは九変数指数から基底を作り、全3,969積を再構成した。
下方閉性と全250,047三重積の結合則を確認した。

socleは正確に32次元であり、plane四次5個、zi*plane二次21個、active zz6個からなる。
最大idealとの全288積とsocle二乗の全1,024積は零である。
残る31スロットは plane≤3の10個とzi*plane≤1の21個で、既存入力で全て決まる。
未提供の真のimplicit係数を零だとは主張していない。

plane10およびB27への二つのquotient mapを、各3,969積で独立検証した。
どちらも未知socle32を全て消す。
既存joint105のordinary `(r,s)` 係数を `s=r+ell/5` で変換し、
既存Q20 B27との630個のplane≤2 overlapを確認した。
全105個の63スロットlower extensionも、出力を合成せず入力配列だけで一致した。

同じ原点値・既知31スロットを持つPhiとlower extension Lの差は、
`Soc(B63)^105` に属する。原点で単元の分母を持つ各元rational DAG fについて、
\[
f(Phi)-f(L)=D f(Phi0)\,(Phi-L)
\]
が成立する。積の差で現れる高次項はm*Soc=Soc²=0により消え、
逆元の差も `−a0^(-2) Delta_a` となるため、DAGの有限帰納で証明できる。
元delta63=`20144−40922` の105勾配が全て零である独立証明を適用すると、
**delta63のB63像全体**が既知lower cacheだけで決まる。

この安全性はzero-gradientであるdelta63についての結論である。
bA/bB/nu等の一般の出力は、両quotient mapが未知入力を消すことから、
**plane10またはB27へ射影した出力だけ**を信頼できる。
それらのfull B63値や、lower extensionがfull105式を満たすことは主張しない。
後続workerは出力の種類と許されるスロットを区別する必要がある。

独立SOURCE所有上界45,036、入力変換のbefore-operation実測45 bit、
実行0.0602秒、55秒・2 GB内で終了した。
新しいimplicit解・physical jet・scalar出力を求めていない。
head30のfull coherent式、将来の正のheight/所有量、worker/custody、ROOT GOは別のgateである。
元H0・Kbar型・nil三乗の追加結論はない。
再現コードは `code/ambient_B18_C_connecting_B63_safe_output_source_review_20261105_v63.py`。
