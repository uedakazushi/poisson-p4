# 完全な実三次判定operator：3×917のsource-only構成

2026-11-05。値を計算する前の完全なsource operatorである。Mはaccepted all-order minimal5×849、Mdはその普通の斉次Taylor係数。
追加source/gauge/jetやM2/M3値は計算しない。

## 定数と一次sourceの全経路

M1のrankは35、定数kernel Cの次元は814。既存661 sourceと153 d3-only sourceを併せ、15 mixedも保持する。
Π6M2Cのrank3を消した定数kernel Dは811次元（658+153）。
二次NFを消すことと二次ベクトル全体を消すことは異なる。
後者はL=(M1R35)·S1の174次元像に戻し、固定H174で280個の一次source係数へ持ち上げる。
さらにLの106次元kernel K106を全て保持する。

三次known closureの商は3次元で、実sourceから得た三dualをΔ3とする。
その代表はtarget0/3/4のrA三乗（600座標119,479,599）で、Δ3はこの三代表上でI3になる。
完全追加operatorは

`[ Δ3(M3 C D − (M2 R35)·H174(M2 C D)) | Δ3(M2 R35)·K106 ]`。

左811列と右106列、計917列である。積の後半はdegree2 relationとdegree1 coefficientの普通のmonomial convolution。

任意source a=a0+a1+a2+…に対して、低次零条件はまずa0=CDd、次に
`a1=−R35 H174(M2CD)d + R35 K106 e + C v1`。
C v1の三次像は既知二次relation×変数に入り、M1a2も既知一次relation×二次monomialに入る。
従って上の917列以外に三次initial relationを与える経路はない。
これはM2/M3の値を既知としたという意味ではない。

## 完了条件

`HF3=3−rank(operator)`。既証length16とrA-reductionから、rank0/1/2なら全HFは順に5533／55321／553111となる。
rank3は残り長さ3と矛盾する。この有限判定を解けば、さらなる高次HFのためのjetを自動的に要求しない。

## 最高source/cache gateはまだ別

δM3の最高線形項とδM2の二つのbilinear項、および106carry項をexact typed operatorとして保存する。
primitive sourceへは同じbranch/gaugeを通してcoherentにchain ruleを適用しなければならない。
raw topにはX3、Hが一つziで割るwedgeにはnaiveにrow20 X4が現れる。
それらやm3を零に置かず、全917列のprojected感度で消すか、保存済みの選択functionalへ明示的にbindする必要がある。

accepted JOINT135は5個の第三monomialに対する152個の選択adjoint（83 gaugeを含む）であり、
旧619用のimplicit −rhoRHS項だけである。full135 source vectorでもnative M3全体でもない。
同じ100次元という名称でも、joint B25はbase20+第三triples、gauge120 B25はbase20+z6 first+四second pairsなので、tupleで照合する。
このpacketはそのprimitive bindingを未決gateとして明示し、新しい係数jobを許可しない。
