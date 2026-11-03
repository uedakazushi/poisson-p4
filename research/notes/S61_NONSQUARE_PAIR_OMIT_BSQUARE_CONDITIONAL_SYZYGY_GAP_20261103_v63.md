# Strong pair から五式linear-syzygyへの未解決の段階

2026-11-03。凍結済み `S61_NONSQUARE_PAIR_OMIT_BSQUARE_20261103_v63` の§2の推論と§3への適用は、独立監査では受理しない。全aPJのactual2×5 rank2定理は別に検証済みで変更しない。

五式のlinear syzygyで、G2とpair平面への制限により Q0,Q3 のcofactor K0,K3がspan(A,X,W)に入るところまでは正しい。S²,SU,U²の係数には scroll三式Q1,Q2prime,Dの寄与がない。しかしB*S,B*U係数には、scroll三式のcofactorのS/U部分が寄与する。例えば

`S*Q1=A*B*S/6075+W*C1*S`

である。従ってこれらを落として「retained2×5 matrix × (K0,K3)=0」とは言えない。実際のscroll moduleでの消去または別の恒等式が必要であり、まだ提示されていない。これはactual five-syzygy命題の反例を主張するものではなく、証明が不足しているという判定である。

§3の、四式quadratic-cofactor syzygyをG1/G2とKoszul減算で五式linear syzygyへ写す構成を用いても、上の未証明段階は残る。従ってG1/G2だけからraw84×126 rank78を昇格しない。もし五式のlinear syzygyが二つのscroll syzygyのみと別途証明されれば、残る復元の検討に進める。

著者は2026-11-03の独立監査中に寄与の欠落を認め、§2の結論と§3への適用を撤回した。元の公開bytesは保持する。既存のproduct41の証明は五式の**定数関係**を排除する別の補題を用いるため影響を受けない。追加のternary必要locusの有限性も、仮に受理されてもこのscroll-module段階を自動的に補わない。
