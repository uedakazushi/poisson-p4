# Top12のE0二次成分は零：新しいjetを使わないexact pointer証明

2026-11-04。pendingなJinf最小generator論に必要なorder gateだけを証明する。
このpacket単独でT3=0や元B18のdepthを主張しない。

正確なHONEST10 top12列を読むと、唯一のunit pivotが消した旧free generator5の成分は、
十二列全てで**literal DAG node0**である。
単位Schurの新E0は

```
new_E0 = raw_E0 - (pivot_E0/unit)*old_component5
```

なので、変換後のtop E0は全てraw E0と厳密に同じnodeである。
全十二pointerの等号を確認した。pivotの微分を求める必要すら無い。

元のtrace pair積の定義から

```
raw_E0(top_(alpha,e)) = nu_CA * j_e * b_cubic_alpha.
```

各積についてDAGの可換積signatureで三factorを確認した。
nu_CA、二つのpositive odd coefficient j_e、六つのrelative cubic odd coefficientは
全て原点値0である。従って **全top E0はm_R0³に入る**。
とくにE0の二次jetは零であり、新たなsource三次やhonest二次計算は不要である。
より強いm4 orderはこのgateでは必要とせず、主張もしない。

既存FULL9 FIRSTJET証明書の全90係数と四g感度を読み、top12の全成分がm²に入ること、
全2069 relationのE0一次成分も零であることを同時に確認した。

従って得られたtyped orderは

```
top12 ⊂ m²*R0^9,
E0(top12) ⊂ m³,
E0(Kinf) ⊂ m².
```

最後の式はKinfが2069関係とtop12で生成されることによる。
これらの事実を、独立に検証するJinf kernelの最小generator43本の一次部分と
残る一つの二次部分へ結び付ける手順は別packetで扱う。
一般にはベクトルのorderだけからmodule membershipを推測できないので、
その最小性・残り二次classの独立性の証明を省かない。

Source `code/ambient_B18_T3_E0_order_gate_20261104_v63.py`、
`reports/AMBIENT_B18_T3_E0_ORDER_GATE_20261104_v63/`に全pointer/factor/originと依存hashを保存。
既存の有限DAG・証明書を読むだけで、GB、jet、点、source再生成は無い。
