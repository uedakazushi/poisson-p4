# 固定六根・62完全積workerの実行前custody

2026-11-05。準備のみ。根・二次factor・Hensel/root coefficient・candidate値は実行しない。
ROOTのsource/content/絶対既約性定理監査
`S61_TWO_TAIL_H_ABSOLUTE_HENSEL_ROOT_SOURCE_REVIEW_20261105_v63`、
SHA `523ac8adfaec65208092de73dcb724693b2b084f4c07df8a4062a29f51865113` を固定する。

値workerは外部factor APIを呼ばない。既存f0と既存degree4 gcdを再びexact monic divisionし、
固定四次に対するFp全要素のHorner根探索を一回行う。この範囲は有限体内部のfactor algorithmで、
新しいe/d/prime、physical/beta sourceを呼ばない。四根の全productと二次cofactorのproductを確認する。
二次discriminantのLegendre値−1、p=32003の2..178 trial division、六つの相異根を全て保存する。

Fp²の積は五つのbase products。Karatsuba crossの両和は必ずmod pにしてから掛ける。
従ってfield一時値30 bitを維持する。二次fieldの各inverseはconjugate/norm/base inverseと
積1を記録する。根seriesでは六つの固定derivative inverse以外の未知係数を反転しない。
alphaと55の単元性、六root identity、L*product6=Fbar mod e10も全係数で確認する。

全64 subset productを保存するが、62 candidateのG/Q/fullproductは一組ずつreleaseする。
producer自身は全133 coefficientを比較し、各失敗について非零coefficient witnessを一つ保存する。
一致する一組があれば固定Fp²像のfactor pairまでで停止し、char0可約性とは呼ばない。
全62組が不一致ならmodp絶対既約性とHのgeometric既約性のPRIMARY証拠であり、独立peerが必要。

実装はpure Pythonの明示finite-field演算で、単一whole-job counterを使う。
hard2百万base products/24kowned/4MiB/45秒soft60秒hard/2GB。
ownedはstate内のcoefficient・tuple/dictionary/controlを数え、6000のlocal/parser/serialization reserveを加える。
全stateはstreamで一つのprefixに保存し、pending scalar product、root/order、poly convolution indices、
current candidateと全過去subsetを保持する。失敗時のreportと全prefixを保存し、自動retryを行わない。
原Hの40bit入力とsource formal remainderの整数はfield30bit域と別で、以前のglobal384上限を変更しない。

実行前にworkerの--helpだけを確認する。合成field fixtureや実seed値は計算しない。
fresh ROOTの全コード/custody監査と一回GOの前には、このworkerのmainを呼ばない。
irreducibilityの結果によってactual bad Fitting/full57/239/尽くしを結論しない。
