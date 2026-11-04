# H15 全体の full239 検査に対する独立レビューの準備

2026-11-04。**準備段階**。下記の実係数・selected98 matrix・rank gate は
未受領であり、full239 の独立 PASS receipt ではない。

## 固定する対象と source frame

対象は monic H15 の全有限 algebra であり、今回の K7W による rank150
cover の二次 frontier だけではない。D0 の SOURCE preflight は、元の
319行＝quadratic products174、cubic products112、raw true quartics33を
全て保持する。未構成の70+15+4 kernel の完全性や、21個の pure quartic
root が既に計算されたとの仮定を置かない。

物理 Tn frame の全239係数を joint 座標へ移した**後**に受理済み NF117
を適用する。その出力を H²,H,pure の順に \(f_2,f_1,f_0\) とする。
\(\chi=u\widehat Q_{77}+v\widehat Q_{131}\), \(U_2=(50aBC)^2\),
\(\delta=-64000C^2\), \(\Pi_{\mathrm{act}}=D_{\delta7}\Pi_{\mathrm{orig}}\) とし、

\[
\mathrm{Section}=L_{\mathrm{act}}^t\Pi_{\mathrm{act}},\qquad
\Gamma(f)=f_0-(f_2/U_2)\chi_0
 -(f_1-(f_2/U_2)\chi_1)\mathrm{Section}.
\]

L は実 delta7 frame の112×23右逆である。元112 pure tails から Section を
作るときには先に delta7 源変換を適用する。\(\chi_0\)、全 cubic pure
tails、f0 を捨てない。Qhat は `(50/B)^2 F(T34 diag(Tden) v)` で、cubic の
row-unit は FrameAdapter の指示どおり **除算**する。H² source と23個の
H1 section source はいずれも実319 source module に属する必要がある。

従って、全319個の Γ 行が98列を生成すれば、117個の既知 quadratic unit
pivots、H² unit 一個、H1 section23個を戻して元239列が生成される。
この source/target の復元は全行を検査対象とする。特定の selected subset
だけが完全であると先に仮定しない。

## 有限 integral lift の論理

\(V=\mathbb Z_{(32003)}\), \(S=V[z]/(H)\) は V 上自由15。
既存 factor inventory は Hbar の七因子を全て保持し、重複度一、次数
`4,1,1,1,2,3,3`、全積一致を記録する。全 coefficient denominator は
QQ の32003-unit だけでなく、必要なパラメータ分母も S の単元でなければ
ならない。各因子上で逆元の積が1であるとの証明を確認する。曲線の五つの
既知 localizer と全因子の互いに素性も保持する。

各因子の**同じ実 Γ**について98独立列が確認されれば、CRT により
\(\bar\Gamma\) は \(\mathbb F_{32003}[z]/\bar H\) 上全射。
実 Γ の finite S-module cokernel は有限 V-module であり、mod32003 が
零なので Nakayama により零となる。これを QQ へ移せば、標数零の H15
の全有限候補上で元 quartic239が得られる。体性・既約性・被約性を標数零
で仮定せず、特殊ファイバーから因子や点を除外しない。

正の結論の範囲は **H15 の有限候補上**。これだけで P100 曲線全体239とは
しない。一因子の欠損・分母失敗・上限超過から QQ の rank failure は
推論しない。

## 受領後の固定 witness 検査

準備済み helper
`code/s61_nonsquare_P100_H15_fixed_rank_review_helper_20261104_v63.cpp`
は既存 prime、七個の既存因子、著者が保存する各因子の元 selected98×98
matrix と固定 pivot 列だけを入力とする。独自の係数配列 field arithmetic
で、指定順の前進消去・逆元積1・全98pivot非零を検証する。pivot search、
source query、係数生成、factor choice、prime choice は行わない。
factor の既約性・全積・source matrix binding は別に検査し、helper の
非零判定だけを源係数の検証の代わりにはしない。

現段階で helper は実 matrix に実行していない。最終独立 receipt は、
著者の固定係数・分母・source proof と rank witness の受領後に別途凍結する。
