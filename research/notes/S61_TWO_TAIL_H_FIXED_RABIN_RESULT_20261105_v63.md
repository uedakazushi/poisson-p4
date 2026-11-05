# 確定Hの固定六次Rabin検査：完了したnegative候補

2026-11-05。ROOTの一回GOに従い、固定e=0,d=1,prime32003だけを実行した。
一次報告 `S61_TWO_TAIL_H_FIXED_RABIN_ENTRY_20261105_v63/report.json` のSHAは
`8628312b3b266d5dd12107db3f996492eb98f0607b781411c21225bc838fee58`。
`computations_completed=true`、`all_passed=false` であり、guard失敗ではない。

特殊化のa先頭係数18248は非零で、次数6は保たれた。monic像は係数を低次から並べて

```
[5457,11750,555,10369,14038,10960,1]
```

である。六Frobeniusと二つの完全Euclid chainを保存した。
`x^(p^6)=x mod f` は成立したが、固定Rabin条件の二GCDはそれぞれ

```
gcd(f,x^(p^3)-x)=[2668,18535,17561,16979,1]  (degree4),
gcd(f,x^(p^2)-x)=f                           (degree6)
```

となり、1ではなかった。従ってこの固定像は可約であり、この候補からHの既約性は証明できない。
これはQQ[a,e,d]上のHが可約であるという結論ではない。
既受理の厳密共通GCD・二整数恒等式・六因子FULLU上での非単元性は変更しない。

field積8740、実field一時値30 bit、owned2550、約0.018秒、入力込みRAW16129 bytes。
既存H係数40 bit、mod前の入力変換整数scratch最大41 bitというscopeは別の
`S61_TWO_TAIL_H_FIXED_RABIN_CONVERSION_SCOPE_20261105_v63`、
SHA `e88bb64ab174109c60cddf461dd47291aa720442e310343823c752050d860f4a` に固定した。
30 bitを全入力変換の上限とは記さない。global384等の別jobの上限は変更していない。

別のprime・特殊化・factor familyへ切り替えず、この結果で停止した。
Hの既約分解、実bad Fitting、global full57/239、補集合消滅について新たな主張はない。
