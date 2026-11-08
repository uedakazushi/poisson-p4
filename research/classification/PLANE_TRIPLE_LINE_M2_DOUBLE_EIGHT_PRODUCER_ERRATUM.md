# 八源producerの使用範囲訂正：非零G2とP2のaffine particular

2026-11-08。root、別作者quartic_reviewの指摘を受けた訂正。
旧producer/report/proof note/manifestのbytesは変更しない。
対象report SHA256は `5dbe47a0b0ec2266c7809bf36f817e057fdd6b9b20a52685ecebf0fb08cdf161`。

1. 各recordの`actual_P0_zero=true`はm三modeからコピーした誤labelである。
   実assertは`act(P0)==profile['G2']`、保存G2は非零で、全行列の算術はこの非零P0から作った。
   当該flagは数学の根拠として使わない。正しい状態はactualD(P0)=指定nonzeroG2である。
2. NIL0の`actual_D_has_proper_a_factor=true`及びそれを使う`G24_closed=true`は、
   fullP1 kernelとhomogeneous K2だけのa-rowを比較した判定であり、全sourceの証明ではない。
   固定P1上のP2解はJ2 P2=−Jac(P1,P1)のaffine translateで、particularがK2内とは限らない。
   **NIL0全空というproducer partial noteの結論は撤回し、NIL0を未閉へ戻す。**
   次は全P2九を戻したOre又はactualD幾何で処理する。

sem六源のu因子は、全P1 bad u-row零と**全P2九にnormal columnsが無い**ことを使うので有効。
この訂正後、元reportが証明する全source空性は六sem源だけで、NIL0/NIL1の二全源が残る。
全350×30/9、完全核、非零minor及びgroup行列という線形算術自体は別作者が再構成中である。
粗28を実成分とせず、旧false/誤flagを成功へ書き換えない。
このproducer packetはまだaccepted ledger1002又は主稿のG完了に加えていない。
