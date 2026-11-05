# native 二十三 pivot の二次九十二係数：独立 SOURCE/code 検証

独立検証は PASS した。報告 SHA256 は
`31668953eb6faf273d29f9698f0e66bf4a840f30535fc68b6306012f27dd1844`。
SOURCE `ab2e8d90…`、未実行 worker `910b6d8f…`、custody `836b2126…` を固定した。
このレビューでは新しい符号付き native pivot 係数を一つも評価していない。

元の二十七 source 列のうち OLD 四列604,613,622,631を分離し、
native 二十三列を `640,644,648,652,656,660,664,752,…,767` と照合した。
32 native source の最初の28列は、八 G 生成元に対する全二元 Koszul 列
`G_j e_i−G_i e_j` である。各符号と順序を全組合せから独立に作り直した。
残り四列は free source8,…,11の標準列である。
二十三列は前七 Koszul 列の mask0と四 free 列の全四 maskに一致する。

実 native 行列は `Res4(N^t) Res4(S12x32)` である。
N の添字は `N[source12][physical12]` であり、転置するのは十二行列だけで、
四 trace 成分の乗法 block を転置しない。
元 Q48 列の実 source 結合、`Q Res4(N^t) Res4(S)=0`、同じ current27/P27/Pi3を固定した。
Gamma は元56行の covectorを保持し、実際に row55 成分を含む。
native X の下八行が零であることにより当該積でその成分は寄与しない。
Gamma 自身の支持が48行以内だとは仮定しない。

N の point/r/rr は元 B11、ell/re は元 C16の受理済み cacheから照合した。
両 cacheの point/r と元 physical C4 adapterの pointが一致する。
G8 は独立確認済み canonical C68 の同じ branch・free4・canonical30 frameを用い、
352 B11および128 C16係数の全 overlap証明を固定した。
plane ordinary 基底は `(1,r,ell,r²,r ell)` である。
monic の初めの三係数は零で、二次は
`nuA=(-25/192)r²`, `nuB=(725/36)r²+(41/9)r ell`。

ordinary 係数での rr は `L0*Srr+Lr*Sr+Lrr*S0`、
r ell は `L0*Sre+Lr*Se+Le*Sr+Lre*S0` の和である。
個別列で S0 を零に置かず、free16列では S が定数なので最後の項を残す。
L と S の双方で monic が変化するため、L二次の N0*nu二次 collisionも保持する。

`tau^a tau^c tau^b` の全64 formal routingを両結合法で独立照合した。
各 monic 因子の指数は高々1で、nuA*nuB は位数4以上なので二次には寄与しない。
worker の routing 関数だけを ASTから抽出し、無関係な synthetic monic を用いた
全80 base/mask fixtureも確認した。実 N/G/nuの積や応答値は計算していない。

各入力 class の共通分母 D と整数絶対値 U を独立に再構成した。
一応答の経路数は粗くても `2*12*4*4*5*5*2=19200` で、SOURCEの2^20以内である。
Gamma/N/G/高々一 monic の全経路と部分和は同じ共通分母の約数を持ち、
正の U 積と2^20で抑えられる。これにより committed280 bitを得た。
二項有理数演算の未約分一時整数は561 bit以下、保守的三因子・最終 envelope比較も
841 bit以下で、hard384/1153を変えない。
worker は各計算後、mapや出力へ commitする前に280 bitを検査する。

正の支持だけから GammaL の所有量312、current intrinsic S の最大36を再構成した。
source payloadは metadata/boolを含め3,883、二重保持7,766。
parsed N/G/nu/Gammaは2,880/160/10/32。
92応答と322内訳の入力・出力複製、cursor・scalar・table予備を加えると、
通常14,024、guardの prefix保存14,372で、21,034以内となる。
worker は標準 libraryだけを使い、全 literal hash、source peer、self hashを検査する。
範囲外時は GammaL/current S/完了応答を保存し、自動 retryや cap変更を行わない。
RAW上界2,360,252 bytes、60秒・2 GBの既存条件を維持する。

初回の reviewer assertion は、独自 counterが boolを除く3,856を返すのに対して、
producer counterの3,883を直接期待したため停止した。27 bool labelの違いを診断した。
停止時 codeと diagnosticを `FIRST_COUNTER_SCOPE` に保存し、
boolを含む数え方を明記してから全 SOURCE/code検査を完了した。
元 primary、worker、custodyは変更していない。応答値はどちらの監査実行でも求めていない。

独立 SOURCE所有上界は42,638、実行0.187秒だった。
本 PASS は同じ一回の限定 jobの準備条件だけを確認する。
全 matrix response、honest quotient、OLD79、outer OLD四 pivot、N三次、617 upper、
HF・長さ・元 H0について新しい結論を与えない。
再現コードは `code/ambient_B18_native_pivot23_second_source_review_20261105_v63.py`。
