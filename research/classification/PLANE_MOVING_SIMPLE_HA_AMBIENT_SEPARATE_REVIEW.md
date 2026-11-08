# HA_1,t=2/15：同一canonical点の全ambient独立確認

2026-11-08。producer `/root/quartic_review` とは別作者 `/root/generic_exhaustion`。
**全算術と二十三次元全成分性を独立PASSとして受理する。**
全W六のactual一軌道定理は別のpure依存、全prior128との相違及び新加算は統合担当／producerの別監査であり、本receiptでは確認していない。

## 独立構成と結果

producerのFraction/q engine、実行program又はsource helperをimportしない。
既存の別算術 `classification_three_hyperplanes_quadric_independent.ambient_matrices` と
`classification_rank3_mixed_component_review` のinteger monomial Jacobi／SymPy QQ matrixを用いた。
座標(a,b,c,z,w)で同じactual点の70倍を、reportの150vectorからコピーせず次の八括弧として直接入力した。

```
Pab=70w², Pac=−70a²+270bw, Paw=−45z²,
Pbc=140ab, Pbz=70z², Pbw=140zw,
Pcz=−140az−42w², Pcw=28aw+378bz,
Paz=Pzw=0.
```

全150 coefficient order／全350 equation orderとactual点を照合し、元全350 Jacobiが零であることを再構成した。
全150方向の微分D、Euler25、実GL5作用25、tensor scale一のM51を定義から再構成し、
**全matrixの全entryがproducer V2 reportと一致**した。
そのうえでproducerの同じ102/24/48 minorのexact nonzero determinantをSymPyで再計算し、
全三行列の完全right kernelを元全rowへ代入して零、free-coordinate identityから独立性を確認した。
新point、normal、minor/profile選別、jet、generic数値rankは用いない。

```
rank_QQ D350x150=102,
rank_QQ Euler150x25=24,
rank_QQ M150x51=48,
D*M=0 (全350×51).
```

fullEuler-Pfaffianも独立再構成し、70²を除いて

`−3w(w²+3az)²/5−27bw²z²/14−27abz³/5−9cz⁴/14`

に一致し、producerの全六monomial sparse coefficientにも一致した。

## 成分性の厳密な矢印

元raw150 Poisson schemeの接次元は48、Euler translations／GL5／nonzero scaleの実作用軌道も48次元。
従って局所Krull次元と接次元が48で一致し、この点はregular/reduced/smooth、唯一のfull ambient component上にある。
Euler quotientは常に24次元線形fiber、非零projective scaleは一方向なので、射影局所次元は23。
actual PGL5 orbitの射影次元もM48−24−1=23である。
同次元のactual orbitはその成分でdense open、従ってorbitclosureがfull23componentになる。
これは選んだquadratic slice又は一部のJacobiだけからの成分性ではない。

全W六の点がこの同じactual orbitへ入るというpure正逆定理は
`PLANE_MOVING_SIMPLE_TWO_WHOLE_ORBITS.md` に別に証明されており、統合担当が受理した依存である。
本receiptの新独立算術は同一canonical点だけで十分であり、二つ目のgeneric点を追加しない。

## 有限実行identity

新jobのgate PASS後、180秒/3GiB上限で一回実行。
1.228525秒、最大RSS51,172KiB、全PASS。
exact inline reviewer bytesをinput archiveとreportに保存した。
shared6022 source/current jobsを編集・rebindしない。

保存reportは `CLASSIFICATION_PLANE_MOVING_SIMPLE_HA_AMBIENT_SEPARATE_REVIEW_20261108.json`。
同一実行で、producer V2の全三fullmatrix、same minors、complete kernelsとactualDを確認した。
V1 producer metadata failureを成功扱いせず、ae5e81…のV2だけを照合した。
新成分ID／prior128相違／全moving-simple点被覆の受理は、この成分性に加えた別のscopeである。
