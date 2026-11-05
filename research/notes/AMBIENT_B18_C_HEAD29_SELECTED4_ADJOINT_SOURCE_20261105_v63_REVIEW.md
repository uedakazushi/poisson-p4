# delta63 の原点勾配零と四次 implicit 解からの独立性

独立SOURCE検証は PASS した。報告 SHA256 は
`d3af22d22e56cf914ce271b5f08b9c6b7332e7996a231ca1ae2dbaecd45ddd4c`。
primary `0a5fffe3…` と全所有量補完 `d7f00043…` を固定した。
delta63_r²ell²そのもの、RHS4、lower-input四次値、bB63_ell³、head29応答は求めていない。

実 physical63 の deltaは同じ元105/F9 frameの scalar root20144から40922を引いたものである。
gA505938やf637790との取り違えはない。元C2のroot locatorと符号も別のsource読みで照合した。
未知105個の元node順序、式105個の元node順序、実原点値、CA target
`(0,…,0,r,s)` は受理済み joint cacheと一致する。

authorの reverse ADとは異なる方法として、元scalar DAGの必要4,387 nodeを、
全105方向の sparse first-jet algebraで前向きに評価した。
各nodeには原点値と独立な105変数の一次係数を持たせ、last-useで解放した。
元の有理数定数はn/dとして読み、逆元nodeの原点が1であることを確認した。
最終的に、20144−40922の**全105個の原点勾配が厳密に零**となった。
同じ保存済みJ105と逆行列の両積を独立に確認したため、rho=grad(g)*J105逆も零である。
新しい逆行列やimplicit解は計算していない。

Phiの二つの持上げが三次まで一致すると、その差は四次から始まる。
原点で正則な同じ有理関数gのTaylor式から、四次の差は
`D_Phi g(Phi0)*Phi4` である。全勾配が零なのでこの差は零となる。
従って delta63_r²ell²は、既知Phi≤3だけを入れたlower-input値に等しい。
これは**真のPhi4が零だという主張ではない**。
そのlower-input値も、このレビューでは評価していない。

own実装は、未使用の18,584-record tableを必要closureだけへ落としてからADを行い、
ADのpoint/gradient/last-use bufferを解放してから固定Jの積を確認した。
前向きphaseの事前所有上界は75,524、実liveは128 node、
before-operation有理数上界は実測16 bit、実行0.269秒だった。
producerの旧123,655は主算術配列の値であり、wholeと解釈しない。
別のcounter補完が借用recordのmetadata、全controller、temporaryまで課した
149,030<150,000を与えることを確認した。元primaryのbyteは変更していない。

今後lower-input値を求めるには、独立したpositive source/profile、全後処理上界、
immutable worker/custody、ROOTの限定GOが別途必要である。
元H0・Kbar型・nil三乗についての追加結論はない。
再現コードは `code/ambient_B18_C_head29_selected4_adjoint_source_review_20261105_v63.py`。
