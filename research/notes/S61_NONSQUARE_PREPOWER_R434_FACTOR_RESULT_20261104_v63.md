# R434 因数分解の保存結果

固定された434項のRを一回だけnative QQ factorに渡した。FLINTの完備factor routineは
20次・434項・重複度1の因子一つを返した。constantを含む全積を元のRの全434係数へ
掛け戻して一致した。従ってこのprimary algorithmの結果はRのQQ既約性を与える。
独立な既約性証明や別アルゴリズムのfactorはまだ行っていない。

既存pの全次数は6で、返却された唯一の既約因子は20次なのでgcd(R,p)=1となる。
これは全pによる一度の非exact除法だけでは分からなかったproper factor共有も除く。
ただし(R,p)が多変数環で単位idealとは結論しない。Delta5、h22、Sauxとのgcdは
計算せず、Rの零集合を既知例外へ同定もしない。

実行はcompile込み0.681秒、native0.0032秒、accessible support最大18060項。
library内部支持はinstrumentしておらず、全processの2GB/60秒制限を保持した。
最初のlaunchはdynamic loaderでlibflint.so.19が見えずmain前に止まったためfactor呼出しは0。
失敗logと旧driverを保存し、既存Cloud activateを実行launchにも適用して、未使用の
一回だけのfactor予算を実行した。数学計算の再試行ではない。

primary report SHA256:
`d7b711aec63e1c3335b476078d8e44ae33ce38b80f7dc12ced97143d686de01d`。
native SHA256:
`9a8b88511b83759a1deb17e378be4e714f27364a612ff3edaa09f95dba7e4c61`。
source/全coefficient/重複度/積は同名 `S61_NONSQUARE_PREPOWER_R434_FACTOR_20261104_v63`
report folderにあり、pre-main failureの全hashは追加のlaunch-provenance.jsonに固定した。

先行 direct P5＋H5は別のsource恒等式によりD(Delta5*Saux)上でRを反転しない。
因数分解結果から5×12のrank5やfull FULLU cubic57を結論しない。
