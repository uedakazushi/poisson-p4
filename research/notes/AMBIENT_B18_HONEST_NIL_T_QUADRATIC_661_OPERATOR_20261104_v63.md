# 二次層の正確な残余operator：六行・661 source

2026-11-04。既知二次初期span174と、d3-only153本の六商像零を結合する。
真のHF2は未決だが、未計算のactual scalar functionalを正確な6×661に縮める。
新しいsource/gauge jet、X3、m2、scalar値は一切求めない。

実M1は849列・rank35で、その35選択basisは全て最初の681 X列にある。
X-only kernelは646次元。d3一次行列のrank15のbasisについて、X側35basisを引いて
打ち消すmixed sourceが15本ある。さらにd3-only kernel153本がある。
この646+15+153=814本は、free-X・pivot-d3・free-d3という三つのidentity blockで独立となり、
ker M1の全体を尽くす。単にd3-only153が零だから168本全てを捨てることはしない。

d3-only153像は既に実sourceから零と確認した。従って、保存した六quotient dualを
全次数actual minimal matrixの二次係数と残る661 sourceへ適用する有限operatorにより、

  HF2(T_nil)=6−rank(actual evaluated 6×661 operator)

と正確に書ける。形式係数はactual matrix entryへ作用させてから取り、finiteNFを微分しない。
全180×849 scalar係数を展開したわけではない。661がrank最小source数という主張でもない。

六dualはwedge部分をr二乗だけで読むので、最高sourceの候補はrow20 X3の
z0r二乗、z4r二乗、z5r二乗へ限定される。ただしunit pivot補正とm2row7の輸送が残り、
既存plane/univariate cacheがこれらの全方向を与えると見なしてはいけない。
新しい係数計算前に、正確なcallback相殺または既存cacheへの完全bindingと全予算を要する。

完全known module GBを得ていないので、有限annihilator・長さ上界やnil三乗の結論は追加しない。
