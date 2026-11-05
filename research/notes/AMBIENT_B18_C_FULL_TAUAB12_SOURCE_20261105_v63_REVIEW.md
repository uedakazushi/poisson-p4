# 全 tauA tauB 四次の形式 source・実装監査

独立報告 `183f613e...` は主 SOURCE `2debebb7...`、全実装 gate `fabe4a9b...`、
固定 worker `2e189705...` を実値前に確認した。
39個の既知 scalar 入力は、出典照合だけ行い、応答の計算では独立形式変数として扱った。
954-node 回路を自前の QQ 係数多項式へ展開し、別に physical 式から作った全 WA/WB、
18個の full NF、15個の residual 多項式と係数ごとに一致させた。

A の原始四次は -nuB2*bB2、B は -nuA_B2*u2-delta2_mixed*E3 であり、
nuA_B2=nuA2-nuA_z beta2 による beta² carry と E3 の負の mixed carry を保った。
H 三次の全四係数を、(delta30*r+delta21*ell)P0、lambda12*ell P1、
mu03*ell P2 で消す。最後の項と mu03 の二乗も省略していない。
未知 H 四次を二次 pivot で消す補正の A/B carry は五次からなので、
delta4・mixed source3・beta3 は不要である。

残る九つの nonpivot source は、実 bB2 と u2 の零、全 plane delta3 の既証零、
active mixed delta2 の文字通りの零によって四次 class が零となる。
inactive mixed 葉は零と仮定せず、対応する E_i3 が零であることを使った。
second-order の純平面 z rewrite は diag(0,-beta_i2,0) で可換であり、
z を含む補正は四次の曖昧さに五次以上でしか寄与しない。

展開前の整数 support 上界は6544、全所有量113222で、独立展開の実ピークは4611。
形式係数の演算前最大72bitであった。最初の preflight では Path 変数を loop の整数で
上書きして hash 保存前に停止したため、失敗を ATTEMPT1 に保存して修正した。
この失敗で応答や形式展開は実行しておらず、主証拠は変更していない。

値 worker は五ファイルを一つずつ読んで39葉を検査する。
全954正上界を独立再構成し、H86・回路171bit、固定 minor・両側逆行列と
全四列 Cramer multiplyback の演算前1553bitを確認した。行列式の大きさの下界は使わない。
全体所有量は静的45000以下、hard60000、全体150000以下である。

親の条件付き一回実行条件は本 PASS で満たした。本監査では応答を評価していない。
residual 非零なら全 HF4=4、b=2、e は三又は四。
ell⁴ 商にも残る行1、3、4の residual が非零なら、下二行を任意のまま保つ
四列 minor v0*Wk-vk*W0 が非零となり、C4 長二十三・e=3 となる。
全零なら既知 HB と十二 source の陰性だけで、全四次の完全性を主張しない。
