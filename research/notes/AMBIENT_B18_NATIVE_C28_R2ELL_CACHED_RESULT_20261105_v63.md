# canonical native の rA²ell 係数三十四個の独立確認

元の canonical unit-IFT 枝、actual kappa9、CA F9 の物理 R0 と ell=5(sA−rA)、trace basis (1,tauA,tauB,tauAtauB) を保持する。P7 の基底は1,rA,ell,rA²,rA ell,rA³,rA²ellで、raw33からの全1089積を確認した。二 monic のplane定数・一次が零なのでbase/trace全重み一の最大基底重みは五、radical六乗が零である。従って元 rational root の P=5+d² の逆は三項で正確である。

下表は ordinary coefficient であり factorial を掛けない。元 C132 cache だけを使い、新たなphysical source/gauge解は求めていない。要求された八 G8 の四maskと二scalar monicだけを一回の固定workerで求めた。

| G8 | 1 | tauA | tauB | tauAtauB |
|---|---:|---:|---:|---:|
| v4 | 0 | -8113/1296 | -383023/20736 | 0 |
| v3 | 0 | 330529/25920 | -109197329/1658880 | 0 |
| v5 | 13751/1152 | 0 | 0 | -634559/5760 |
| v7 | 0 | -2395/72 | -593/144 | 0 |
| v8 | 0 | 0 | 0 | 1066493/20736 |
| Q0 | 0 | -2189/648 | -729395/10368 | 0 |
| Q1 | 0 | 1429781/1296 | 163762769/82944 | 0 |
| s*r_core² | -11275/72 | 0 | 0 | -33200869/69120 |

二monicの rA²ell 係数はともに零である。SOURCE の正上界37 bit committed・112 temporary、全保持14403を維持し、実最大31 bit・保持6466・12556積・0.054秒で完了した。

ROOT の独立checkerは展開済みproducer polynomialを使わず、元の五rational root、geometric inverseの冪、core preparationと別のtrace積表を用いて全34を再構築し一致した。値の再実行はこのreceiptでは行っていない。

個別 native pivot の L二次S0、N三次、Xc三次、honest g/h quotient、617 upper、HF、元H0・nil三乗・T2全長への結論は含まない。次の coefficient-only native pivot SOURCE は別packetであり、まだ値を求めていない。

証拠: `reports/AMBIENT_B18_NATIVE_CANONICAL28_PLANE3_VALUES_20261105_v63/report.json`、独立確認: `reports/AMBIENT_B18_NATIVE_C28_DIRECT_REVIEW_20261105_v63/report.json`。
