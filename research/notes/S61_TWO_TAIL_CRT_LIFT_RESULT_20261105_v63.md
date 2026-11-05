# 固定二tailの94項共通因子：厳密な二恒等式の一次結果

2026-11-05。ROOT実行前監査後に一回だけ固定40prime CRT pipelineを実行し、完了した。
一次報告 `S61_TWO_TAIL_CRT_LIFT_ENTRY_20261105_v63` の SHA は
`f2ff75ad00931757c3f0eade5d3250d6d265199f12d86e69e10fd1327e4c32b6`。

元の4500項N0と2317項N2に対し、primitive整数多項式Hと整数A0,A2を得た。
\[
N_0=H A_0,\qquad N_2=H A_2.
\]
二本とも整数除算の全余り零に加え、元三変数入力へのfresh全係数掛け戻しが零である。
単なるunivariate像の除算や正gcdからの推論ではない。

Hは94項、最大40 bit、変数次数(6,9,4)、全次数11、lex先頭単項式a⁶ed⁴で、先頭係数1160156250。
primitive content1と正の先頭係数を保存し、Hのmonic32003還元は既存seed94項と一致する。
A0は1795項、最大252 bit、次数(10,23,14)、全次数38。
A2は724項、最大229 bit、次数(6,16,10)、全次数27。
共通の元D197は変更していない。

値は `factor-sparse.txt`（a,e,d,integer）と `factor.bin`、二つのcofactor binaryへ保存した。
binaryはmagic8、input index uint32LE、count uint32LE、各key uint64LE＋sign byte＋length uint32LE＋little-endian magnitude。
keyのlegacy部はkey mod2²⁶、W部分はkey div2²⁶=2a+e。
H text SHA `241488fafa9b45e5353f0b707f75902bc5bc32e85ee093b9f2e0a5f07c7d571c`。
H binary SHA `82b792eb63e301fe912ca7550e679094643bfe6ee30da06a9e314bc1469c52ec`。
A0 binary SHA `d6f597c7cc58d518695735e379fdbf51a5f2500a99ae4ad9fc34d021a11aa762`、
A2 binary SHA `289c3cc9b386b9d47a9b2ddcb1745b591afb196d1216a58f1eaf8799a2647178`。

全40primeを受理し、uint modular2638783719 products、bigint632130、peak52055/60000、
高さgate664/555 bit、RAW1137878 bytes、約10.21秒で完了。追加prime・retry・code変更はない。
共通因子の一次結果までを記録し、独立恒等式peerは別に要求する。
既存modular quotient-coprime peerからの最大性lemmaは、独立identity受理後に別途扱う。
Hの既約性、実bad Fitting、全FULLU57/239、補集合消滅や新しい一般成分は主張しない。
