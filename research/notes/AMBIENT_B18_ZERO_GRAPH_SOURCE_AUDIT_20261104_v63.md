# 四つの零 graph の source 監査と、未解決の正確な membership

2026-11-04。source/theorem 調査のみ。高次 source jet・core 三次・GB は行わない。

## 1. 普遍的な物理恒等式だけでは四 graph は消えない

元の零 weight canonical target は、40個の ey-wedge quadratic field 係数に四つの
divergence 線形条件を課した36次元空間。15個の実共通切断を固定すると21座標が残る。
その literal な座標系は R0 九個、零 odd 出力八個、graph rows64,68,73,74 の四個である。
既存の4x4 divergence unit adapter から、選択17座標の変分が全て零で、四 graph の
一つだけを1動かす四つの独立 canonical 方向を構成した。
全方向は divergence zero、ey-wedge、y に独立な二次 field であり、Jacobi と射影 rank2 条件も満たす。
従って Euler・divergence・一般の ey-wedge Jacobi/rank2 恒等式だけを使って
四 graph を選択17座標の関数にすることはできない。
これは CA/CB 特有の source 方程式を加えた問題に対する反例ではない。

## 2. 個別 graph 表示と共通 graph 表示の違い

各側では既に

    q_ij=a_ij(R0)+b_ij(R0)*tau_i,
    b_graph in I0_i=(z7,g_i)

が証明され、四 graph は八個の odd 係数を使って正則表示される。
この表示の係数は七 coordinate Hadamard と、非零三次 g_i による正則商を含む。
g_i は中心で単位ではない。二族の四商の origin は零で、個々の一次微分は
零ではない場合もあるが、実全九座標 adapter 後の差は既存の四次監査で相殺された。
これを全次数の相殺に拡張した証明はまだない。

一方の graph 表示を他方に代入したときには、even 部分の差と graph 表示係数の差が残る。
個別の I0-membership は、両族に同じ ambient graph extension が存在するという主張ではない。
凍結済み二族の graph 商を中心で単位反転する操作も正当化されない。

## 3. 消える determinant と、捨ててはいけない even defect

二正規化の自由四基底を (1,tauA,tauB,tauA*tauB) とし、出力差を

    d_i=a_i+b_i*tauA+c_i*tauB

と書く。ここで c_i は CB の odd 係数に負符号を付けたもの。
rows18,19 を i,j とし、任意の graph k に対し

    Delta=b_i*c_j-b_j*c_i,
    alpha=b_k*c_j-b_j*c_k,
    beta=b_i*c_k-b_k*c_i,
    R_k=Delta*a_k-alpha*a_i-beta*a_j

とおくと、普遍的に

    Delta*d_k=alpha*d_i+beta*d_j+R_k

が成り立つ。今回、この恒等式と四つの実 DAG 入力を保存した。
Delta の中心値は零であり、これを単位として割らない。Delta 自身の非零性もここでは判定しない。
R_k を零と置いたり odd ideal に入ると仮定したりもしない。
generic coefficient field の消去では、この二つの問題を隠してしまう。

## 4. これ以上小さくしても Tor を失わない有限完了目標

二つのモニック関係は既に自由四係数代数へ組み込まれている。
八つの odd 出力差とそれぞれの四基底倍を列に持つ、実4x32行列 M を
レビュー済みの全4x48表示から literal に選んだ。
四 graph 出力差の自由四列を G_k とする。全 graph 冗長性の正確な完了条件は

    M*c_k=G_k  (k=64,68,73,74),   c_k in R0^32.

元の純CA48+CB48係数系・同じ R0・全行列ノードを添付に固定した。
これは新しい source 求解や既存 graph 商の大きな展開を必要としない。
ただし完成局所環 R0 内の係数が必要で、分数体での解では十分ではない。
今回この四 membership は求解していない。

Npair 因子 L を用いる正則九環 S では、八 odd 差が L*(五単位座標と三 core 関数)、
四 graph 差が L*G_graph である。S は整域なので L を消去できる。
従って上の実 membership は、四 core graph 関数が三実 core 関数
(K67,K72,K63) の ideal に入ることと同値である。
これが判明して初めて、全 core が三生成になるかを議論できる。
さらに height3 が必要であり、三生成だけから CI1/CM1 を宣言しない。

節点付き初期曲線の generic chart への単射は、関数の零性を調べる助けにはなるが、
この membership・maximal-ideal saturation・socle を自動的には決めない。
node の二 preimage の gluing と上の局所係数条件を保持する必要がある。

結論: 普遍的物理恒等式だけによる shortcut は否定したが、CA/CB source 特有の
四 graph 冗長性は未解決。新しい係数を推測せず、実4x32・四列の完了目標までを固定した。
