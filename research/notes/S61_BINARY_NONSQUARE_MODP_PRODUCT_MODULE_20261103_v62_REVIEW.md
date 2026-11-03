# 非平方product残余の有限体right inverse：独立監査

2026-11-03。凍結 `S61_BINARY_NONSQUARE_MODP_PRODUCT_MODULE` を、**残余行列の標数32003の証明**として承認する。
元のproduct41、QQでの全射、full cubic57やquartic65、ambient除去はこのpacketから結論しない。

## 元QQ残余との係数ごとの同一視

94,478,611-byteの元QQ26×60残余をstreaming SHA256で固定した。巨大ファイルはread_bytesで一括読込みしない。
先頭10行の各係数を別のFraction parserで読み、全分母が32003の単元であることを確認した。
保存された53列との全530entryが直接の係数還元と一致する。省いた7列は**この10行prefixで**元から零である。
後半16行や、それより前の31個のQQ pivotのprime整合性をここで検査したとはしない。

native Mは10×53・205156項、Lは53×10・183182項で、全entry・全exponent・末尾EOFを検査した。
Lの最大total degreeは18である。53行全体が存在し、途中の一行だけを全行列とみなしていない。

## 別native処理による全100多項式恒等式

独立C++ checkerは、作者のlex順序と100個別entry和に対し、degrevlex順序と53個のcolumn outer productを使う。
FLINTのpolynomial演算を用いるが、SingularのGB/lift/discovery処理は呼ばない。
右辺はpretty-string parserを使わず、aP³、C=100a+18e²-9e-2、J=5P-6e-4を個別に構成して掛ける。
全100entryで厳密に

    M L = a P³ C J² I10

を確認した。有限個の点への評価やrank標本ではなく、Fp[a,e,P]の全係数の等号である。
従ってこの残余MはD(aPCJ)上で全射であり、特に作者の指定open aPFB0CJ!=0上でも全射である。
変数Pはserializeでdと呼ばれているが、物理的にはP=adである。

## 壊れたexportと論理的限定

INVALID_first_row_onlyとmaxdeg衝突logは変更せず保存されている。
独立parserはINVALIDを完全native行列として拒否し、別の2×3 per-entry fixtureを正確に読み直した。
修復後の完全53×10データだけで全100恒等式を検証している。

最後の残余の分母がprime単元でも、前の31 QQ pivotsが全てFpへ可逆に特殊化することは自動ではない。
従って元のFp product41へのbridgeは未監査であり、結論しない。
またFpのlocalized全射からQQでの全射を外挿しない。
QQ cofactorや他の厳密証明は別タスクであり、この監査ではGB・CRTを新たに実行していない。

review code: `code/s61_binary_nonsquare_modp_product_module_review_20261103_v62.py` と
`code/polynomial_modp_right_inverse_review_20261103_v62.cpp`。
receipt/log: `reports/S61_BINARY_NONSQUARE_MODP_PRODUCT_MODULE_20261103_v62_REVIEW/`。
