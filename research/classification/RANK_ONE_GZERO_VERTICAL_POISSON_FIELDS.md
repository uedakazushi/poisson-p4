# rank一long・g零の全源：affine vertical Poisson場

2026-11-08、純全source補題。新point／rank／minor／operatorは実行しない。ここでは全軌道の**上界**を証明し、数値点で得ていない下界又はgeneric値を仮定しない。

g零源をU=u∂z+v∂w、d=div U、(U+1)H=0、全quadratic Fで書く。

```
Pxy=−F−(1+d)(y−x²/2), Pxz=u, Pxw=v,
Pyz=xu+H_w, Pyw=xv−H_z, Pzw=0.
```

任意affine式ℓ(z,w)についてYℓ=ℓ∂yと置く。Lie微分の唯一残る成分はxyで

```
(L_Yℓ P)xy=−(1+d)ℓ−U(ℓ).
```

他の成分は係数のy独立性、Yのy以外のcomponent零及びPzw零により零。従ってYℓはPoisson場であることと

```
(U+1+d)ℓ=0
```

は同値である。全FとHには追加条件を付けていない。X=∂x+x∂yも全源のmodular Poisson場で、非零Yℓたちとは常に独立。これらは全てprojective linear fieldsへ延長する。

元finite routingのcenter／eigen frameで、次の独立fieldsが全sourceに存在する。

| whole source | X以外のfields | 全PGL軌道上界 |
| --- | --- | ---: |
| M1、B1、B2、B5 | 今回の一様下界には不要 | 23 |
| M11 | ∂y | 22 |
| B3 | z∂y | 22 |
| B4 | ∂y | 22 |
| B6 | w∂y | 22 |
| V1 | ∂y、w∂y | 21 |
| SC | z∂y、w∂y | 21 |
| M21 | z∂y | 22 |

centerを戻す場合はz,wを対応するaffine中心座標へ、全q-affine及びPGL frameを戻す場合はそのactual pushへ移す。M11/M21の全period、V1の全p(w)、SCの全binary T、全Fを保持する。例外係数でさらに多くのfieldが現れる可能性を除外していない。

今後同じ既固定点でM rankが上界+25に達した場合だけ、非零minor開と既約全sourceによりそのgeneric軌道値が上界と一致する。下界が届かなければ、このnoteだけでgeneric値を宣言しない。これは追加pointを許可する規則ではない。
