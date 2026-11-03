# 一つの既存 NEW25 点による mixed D 非零：独立監査

2026-11-03。対象 `Q_FIRSTJET_SPIN4_NEW25_POINT_20261103_v63` を独立受理する。
検査する点は、受理済み36点 packet の `actual_QQ_source_points[0]` **一つだけ**である。
source parameter は n=(7,9,9)、scale=1、base=(-1,1,1,-2)、保存された恒等 Sp4 chart。
他の35点への評価、新しい点の探索、次数の追加 probe は行っていない。

## 点の再構成を保存 beta16 から独立にした

既存 dense NEW25 map の6079-node DAG を Fraction 演算で新たに評価し、全61座標を得た。
value pivot=2051/12、value Pfaffian=117727/2592 はともに非零であり、元の rational map の定義域にある。
次に既存の正確な primitive-jet/S61係数写像から24係数を直接作った。
これは author point program の「保存 beta16 と primitive24×16 basis の積」と異なる入口である。
得た全24係数は保存 beta16 からの値および今回の author certificate と全て一致した。

この点が NEW25 の first-jet image に属することには、受理済み dense map の元全Poisson族への
同定と、その defined source 点であることを使う。今回だけで新しい成分同定を仮定したり、
一般式の恒等性を一点検査から推測したりはしない。

## 二次体算術と物理 frame

cofactor normal と Gram determinant は、別の Fraction Gaussian determinant で再計算した。

`q=−10025632819168382287872/1019784129324739375=16p0<0`。

よって p0 は非零で、QQ[s]/(s²−q) は真の二次体である。
author の K class、代数的数値近似、point producer は使わず、有理数対 (a,b) の加算・乗算
`(a,b)(c,d)=(ac+qbd,ad+bc)` と norm による逆数を独立に実装した。
Gamma_n²=qI、二つの projector、選択列の symplectic pairing と normalization、frame、
g=frame^-T の全係数を再構成した。g^TJg=J、det g=1、固有値(s,s,−s,−s)を厳密に確認した。

元24係数の pullback は、二形式の行列 loop を再利用せず、

`beta'_{ij,k}=sum g_{lk} beta_{ab,l}(g_{ai}g_{bj}−g_{bi}g_{aj})`

を直接足した。全24係数が一致し、閉形式・primitive の8条件、最初のPsi行零、
cross determinant=s を確認した。12個の二形式係数は、先に独立受理した24×12埋込みの
別の左逆から抽出し、author の係数抽出式と一致することを確認した。

## 判別式の完全な QQ 非零反例

F5 の六係数は二次体の convolution で作った。V の quadratic resultant は4×4 Sylvester
determinant、quintic discriminant は9×9 Sylvester determinantを leading coefficient で割って
独立に計算した。この一点では leading coefficient も resultant も非零である。
したがって保存170-node回路の評価だけを二重に繰り返した監査ではない。

R、Disc(F5)、D=Disc(F5)/R² は全て author certificate と厳密に一致した。
D の s係数は正確に0、QQ成分は正確に非零で、分子420 bit・分母319 bitである。
完全な numerator と denominator は独立 certificate に保持した。
全170節点も別の有理数対演算で再評価し、保存値と共役作用の全一致を確認した。
有限体での非零や数値近似を、QQ反例の根拠にしていない。

既受理の normalizer 偶性により、frame の選択は降下した D の値を変えない。
p0!=0 のこの NEW25 点で D!=0 なので、どの k>=0 に対しても、もし p0^k D が多項式へ延長すれば
その値は非零である。従って **この混合判別式候補と、その極を p0 の冪だけで払った numerator は
NEW25 hypersurface の方程式ではない。**

これは候補の厳密な否定であり、混合判別式・Spin4切片の構成を否定しない。
pole order k、別の判別式、NEW25方程式、定義次数への新しい下界は決定しない。
一点での消滅が仮に起きた場合には方程式の証明にならなかったが、今回の非消滅は候補否定に十分である。

再現入口 `code/q_firstjet_spin4_new25_point_review_20261103_v63.py`。
60秒・2GB上限、author point codeを import/実行せず、元の凍結証拠も変更しない。
