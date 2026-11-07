# 確定H94の直接三変数QQ factor：実行前sourceと条件付き掛け戻し

2026-11-05。準備のみ。因子値はまだ計算しない。対象は確定した同じH94一つ。
三変数fmpz_mpoly_factorを一回だけ使い、別prime/d値/点・modular factor farmを追加しない。
以前の固定modular像の因子をchar0へ推測してliftする方法ではない。

## 入力とcanonical出力の高さ

Hはprimitive整数多項式、94項、入力40 bit、変数次数(6,9,4)、全次数11。
bound専用のKronecker写像a=t50,e=t5,d=tは、このboxの各項と各因子boxで単射。
H像の次数は309、L1 normは2^47未満。
primitive整数因子もこの一変数像のprimitive整数因子になるので、Mignotte boundにより
canonical primitive factorの各係数は357 bit以下という保守的上界を持つ。
ここでunivariate factorを計算・decodeすることはない。値backendは直接三変数である。

installed FLINT3.1のfmpz_mpoly_factorはintのSTATUSを返す。0を「零多項式の因子結果」と扱わない。
get_baseとprimitive_partを用いて整数primitiveな正のlex先頭係数の候補を出す。
canonical係数357 bit、各factor350項、11因子/総3850項までを出力前にcheckする。
factor APIのscalar unitを掛け戻したと仮定せず、positive primitive因子の全積を元Hと厳密比較する。

## opaque部分と384 bit部分を分ける

factor APIとlibrary primitive normalizationの内部高さ・内部supportは計測できない。
ここはnative1.9GB address limitとshared45秒soft/60秒hardで制限するopaque library phaseと明示する。
内部まで384 bitで認証したとは言わない。停止時には元入力・code/pins・phase journal・logを残し、
library内部の再開可能checkpointや未出力因子の存在を主張しない。
canonical外部出力は357 bit。Pythonで行う全verified identityはhard384 bitを保持する。

掛け戻しの各stage直前に、実際のcurrent部分積の最大bit数bpと次factorのbqを測る。
一係数の寄与数より保守的な全350²積を使い、`bp+bq+18<=384` の場合だけそのstageを行う。
single-factorの最大bitだけを二倍する誤ったgateは使わない。
各product/accumulatorをsource絶対上界で囲み、384 bitをcheckしてからcommitする。
gateが通らなければ、全canonical候補とcurrent/next/pending/cursorを保存してfrontierで止める。
capraise・並べ替えfallback・別factor callはしない。

三変数boxは高々350係数。各partialの実degreeと次factorのdegreeの和を確認し、
packed basis indexによるcarry aliasを禁止する。このindexは三変数monomialのaddressingだけで、
univariate因子decodeではない。全非零94係数と残るboxの零係数まで元Hと比較する。

## 全job資源と結論

native accessible cellsとPython全候補/JSON/current/next/metadataの上限50000、
verification pair上限1.5百万、RAW4MiB、shared45秒soft/58秒driver/60秒hard/2GB。
native opaque phaseの1.9GBと、source上50MB未満のdriverを合わせても2GB以内とする。
失敗・timeoutには再実行せず、whole sourceと存在する候補/verification prefixを保持する。
privateAPI内部は復元不能である点をfrontier reportにも明記する。

全掛け戻し成功で得るのは直接QQ factorizationのPRIMARY結果である。
各baseのQQ既約性の完全性にはFLINT factor API契約を使うので、独立peerが必要。
一因子ならQQ既約までで、絶対既約とは言わない。複数の非定数factorの厳密全積ならHのQQ可約性を示す。
実bad Fitting、full57/239、補集合消滅は別問題である。
fresh code/custody/ROOT監査と一回GOの前に、backendの--runを呼ばない。
