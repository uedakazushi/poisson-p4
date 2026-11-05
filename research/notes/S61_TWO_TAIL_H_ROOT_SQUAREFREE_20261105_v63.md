# 確定共通因子Hの無重複性と選択二tailの余次元一部分

2026-11-05。ROOTの独立微分Euclid検証
`S61_TWO_TAIL_H_ROOT_SQUAREFREE_20261105_v63/report.json`、
SHA `49990187b3b0b8edaf8939e4d13ac297f4f2ecba7f8b9f5326639089fa043ce2` を受理した。
同じ固定e=0,d=1,prime32003だけを使い、producerのRabin/Frobeniusとは異なる方法で
fとf'の完全Euclid chainを保存し、gcd=1を確認した。新しい点・primeは使っていない。

Hのa最高次係数は17578125*d^3*(80-25d+66ed)。保存されたa定数係数のe=0制限は
`-1025280-12816d` であり、dにも既約なL=80-25d+66edにも割れない。
従ってHはQQ[e,d][a]でprimitiveで、a次数0の非定数因子を持たない。
a次数6とその先頭係数が固定還元で残るため、Hに重複因子があればその正のa次数も残り、
fに重複因子が生じる。しかし独立微分gcdは1である。
よってHはQQ[a,e,d]でsquarefreeである。これはHの既約性の証明ではない。
先のRabin negative結果は「可約だが無重複な固定像」と両立する。

受理済みの厳密最大公約因子はHであり、選択二tailのidealは、完全な六因子環

```
R=QQ[a,e,d,1/(adFBCJ)],   (Q0,Q2)=H*(A0,A2),   gcd(A0,A2)=1
```

である。Hは六つの反転因子のどれとも共通因子を持たない非単元なので、
V(H)はこのbase上の非空なreduced divisorである。
height1で(A0,A2)は単位idealになるため、選択二tail idealの余次元一部分は
正確にこの無重複Hで与えられる。残るcofactor idealはproperならheightが少なくとも2。
scheme全体のreducednessや、cofactor零点の存在はこの主張に含めない。

この非空性は六因子baseに限る。追加open D(D12)との交わりは、まだ結論を出していない。
選択tail条件が失敗することを、実bad Fittingやfull57のrank低下と同一視しない。
全FULLU cubic57/quartic239、補集合消滅、一般成分の尽くしは未解決である。
本publicationは検証を再実行せず、独立ROOT報告と既存sourceのhashを結ぶだけである。
