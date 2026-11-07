# 固定二変数像の因子対：六根workerの完了したnegative結果

2026-11-05。ROOTの全コード/custody監査
`S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_ROOT_CODE_REVIEW_20261105_v63`、
SHA `c68e2870124f99016edcda9d7a4b21b2e2b8946095a53a63480d9a929add117e` 後、一回だけ実行した。
一次報告 `S61_TWO_TAIL_H_ABSOLUTE_SIX_ROOT_ENTRY_20261105_v63/report.json`、
SHA `3028597f392c2c16d3a1aa51065b59cb5702401df3ed2f93f045742dc2ad4fa6`。

固定d=1、prime32003のFbar=H(a,e,1)/17578125に対し、prime trial証明、
seed四線形＋一二次の全product/二次Legendre、六derivative inverse、全54 root updates、
全60 root jet係数、L*product6 identity、全64 subset productまで完了した。
新しいe/d/prime/physical sourceを追加していない。

最初の二candidate pairは不一致で、mask3・LC allocation0の三番目がFULL133係数で一致した。
a次数2と4の因子対で、両方の係数はFp²の虚部が全て零、すなわちFp上にある。
次数2の側は

```
G=a²+(12770+12794e)a+18179+19971e+11005e²+9216e³  (mod32003).
```

補因子Q、全GQ product、seed/root/subset/candidate状態を
`complete-six-root-certificate.json` に保存した。
SHA `0a9dbf9d3e438ab58c6514cd618fdbf02a1d0b827c555222c94f60765fa6335f`。
この固定像の因子対が見つかった時点で止めたので、62組不一致の証明ではない。

224922 base products、owned13782、field一時値30 bit、約0.197秒、入力込みRAW78103 bytes。
guard停止、retry、別のprime/d値、factor familyの自動追加はない。
これは固定modular二変数像の可約性だけを示し、char0のH因子や可約性とは呼ばない。
Hの絶対既約性は未解決であり、既存のsquarefree/厳密GCD/H∤D12は変更しない。
bad Fitting、actual full57/239、補集合消滅、一般成分の尽くしを結論しない。

次の候補は、既に確定したH94そのものに対する一回のexact QQ factorizationを別予算で検討すること。
このnoteはその値実行を許可せず、新たなfactor/GCD/CRTを走らせていない。
