# R7絶対既約性の14候補：既存四根jetだけを再利用するsource

2026-11-05。準備のみ。新しいseed・Hensel・root coefficient・candidate値は計算しない。
独立確定したH=R7*S4のうち29項R7を対象にする。同じd=1、prime32003、
既存Fp²=Fp[w]/(w²+25984w+7499)だけを使う。

## content・次数・旧Qとの全係数binding

R7のa先頭係数は140625*d*(80−25d+66ed)。a定数係数は非零のeだけの多項式で、定数−6408。
従ってdはcontentを割らず、d次数1の既約なlinear因子もeだけの定数係数を割れない。
R7は代数閉包上でもaについてprimitiveである。
d1制限はa次数4、e次数高々6、a先頭係数140625*(55+66e)を保持する。
modpのcontentは別に、a定数係数の6e+5によるformal pseudo-remainderが非零であることを確認する。

Rbar=R7(a,e,1)/140625と置く。sourceの全35係数を、既存六根packetで一致したcofactor Qと
scalarを払って比較する。Qのe7..9係数が全て零であり、残る全係数がRbarと一致することも保存する。
Rbarのa先頭係数は正確にL=55+66e。
旧matched subsetはmask3、cofactorはmask60なので、再利用する根は旧index2,3,4,5の四本。
その各々のe0..6係数だけをimmutable inputへcopyする。根を再構成しない。

## 一意性と14完全積

値workerはcached二次fieldのLegendre証拠を再検証し、primeと初期四根の相異性を確認する。
新しいroot/seed探索は行わない。四つのRbar(root,e)=0 mod e7と、
L*product4=Rbar mod e7を全係数で確認する。このe0 identityと初期根の相異性から
root derivativeはunitであり、formal root/factorの一意性が従う。

任意の幾何学的因子はa-leadingを1またはLに正規化でき、a-degreeは正。
根はこの四本の部分集合に限られ、formal monic係数はFp²に属する。
各因子のe次数は高々6なのでe7までで正確なpolynomial候補を復元できる。
一根を含むproper subset七個とleading allocation二通り、計14組だけを比較する。
producerはa次数4/e次数12のFULL65 coefficient productを全て確認する。
14組不一致ならRbarはmodpで絶対既約であり、global geometric a-primitivityと次数保持により
R7もchar0で絶対既約となる。positive PRIMARYには別の独立identity/theorem peerが必要。
一致した場合は固定modular sliceのfactorまでで停止し、char0 factorと呼ばない。

## 計算と範囲

rootjet評価4900、全16subset積29400、14fullproducts30870、leading allocationその他を含め
source base products75000未満/hard100000、owned12000、RAW2MiB、45soft60hard/2GBを提案する。
Karatsuba crossの和はmodpにしてから掛け、field一時値30 bitを保つ。
全旧六根再構成や32003-element seed探索は行わない。旧jetはinput、root residualは有限identity check。
ROOT fresh code/custody監査と一回GOの前に値実行しない。

R7絶対既約性が独立受理された場合に初めて、Hの二つのgeometric因子がR7とS4だといえる。
各成分と追加D12 openの関係や、S4が旧Theta/pに相当するかは別のsource inventoryとする。
full57/239やactual bad Fittingの主張はここに含めない。
