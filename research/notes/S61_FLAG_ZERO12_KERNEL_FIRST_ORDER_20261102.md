# zero12実kernelの一次底厚化は自由333次元

2026年11月2日。指定固定β点の実kernel M、R2=QQ[[v0,v1]]、b=(v0,v1)を用いる。
[111実生成元](S61_FLAG_ZERO12_KERNEL_GENERATORS_20261102.md)と
[QQ七次初期一致](S61_FLAG_ZERO12_SEVENTH_INITIAL_20261101.md)から

    M/b²M ≅ (R2/b²)^111,     dim_QQ M/b²M=333

が従う。全MのR2自由性とは区別する。

111個の生成元の重みは1から6であり、F=⊕R2 m_i→MのkernelをKとする。
先頭項の完備除去で、この射は重みfiltrationに関してstrictである。
QQ七次証明は元135式の七次までの全閉包を検証し、J6からの新しい初期関係がない。
従ってE=gr F→gr Mは次数7以下で同型である。

Kがb²Fに含まれなければ、ある成分の係数に底次数0または1の非零項がある。
その成分の重みは高々6なので、関係の最初の非零重みは高々7である。
これは上のgraded同型のkernelに非零元を与えて矛盾する。
よってK⊂b²Fであり、M/b²M=F/(K+b²F)=F/b²Fを得る。

さらに既証明のfiltered交叉により、M/b²Mは
S/(b²S+m_S^8)内の111候補とそのv0,v1倍の張る部分に正確に一致する。
Dの同じ切断の長さは98である。
実際、gr DをP2上で見たHilbert級数は(1+4z)/(1−z)だから、
b²を殺すと(1+4z)(1+2z)/(1−z)となり、その七次までの係数和は98である。
従って

    length S/(b²S+m_S^8)=333+98=431.

これは有限体の標本長431からQQへ移した主張ではなく、既証明のQQ七次と
実平方零商を使った結論である。

[残余作用](S61_FLAG_ZERO12_RESIDUE_ACTIONS_20261102.md)によりtM⊂bMである。
したがってこの自由333商でのt作用は一意な二つのQQ行列T0,T1により

    t = v0 T0 + v1 T1  modulo b²

と表せる。現在、同じ原式からのmod32003計算では両行列が零となったが、
本証明書はT0=T1=0をQQでは主張しない。
その係数確認は次の有限計算であり、全Mの冪零性への自動的な橋でもない。

```
.venv/bin/python code/s61_flag_zero12_kernel_first_order_20261102.py \
  --output /tmp/zero12-kernel-first-order.json
```
