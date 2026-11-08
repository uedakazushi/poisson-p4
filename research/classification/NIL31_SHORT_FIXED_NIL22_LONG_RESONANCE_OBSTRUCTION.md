# shortnil の simple-H 模型：fixed nil22long 源の既約quartic閉包にも入らない

2026-11-08、PRIMARY E、純証明。新計算・source/operator・profileは無い。
全r2/3 E所属は未解決のまま保持する。

`NIL31_SHORT_SIMPLE_FACTOR_SPECTRAL_ENTRANCE.md` の全tensor同一性を、
公開 nil22long 全20係数源へ手で戻すと

    n1=n5=−3d/7,n8=−9t/49,n11=1/7,n14=n19=−2/7、他全零。

これは六括弧の全係数で一致する。divはy∂x+∂y+w∂z。
既potential座標 t0=y,u=x−y²/2,v=z−wy の値は

    Φ=(9t/49)w²v,
    Vu=−2u/7+27t v/49+3dw/7,Vv=−3v/7,Vw=−2w/7。

したがって Hessian rank0、q=n14=−2/7、actualF=Φ。
既証明 `HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_SOURCES.md` は、
**同じ全20係数 Poisson source** の actual irreducible quartic を全て

    q∈{−1/2,1/2,−1,−1/3}

へ必要還元した。rank2 Hessianはfactor、rank1全resonancesとrank0 q−1の全補集合、
q−1/4の単位理想まで独立V2で確認済みである。
従ってその全既約quartic源の閉包も閉条件

    (q+1/2)(q−1/2)(q+1)(q+1/3)=0

を満たす。target q−2/7はこの像外。
regular parametercurveを同じfull20sourceの中で動かしgeneric actualquarticを既約にすることは不可能。
固定 source の q がgeneric有限集合へ入ればwholecurveで同じ閉条件を満たすためである。
generic parametersだけを標本化した否定ではない。

これは **fixed cofactorとfixed full20 coordinates** の全閉包に対する排除。
cofactor、infinity、spectralflagsが動く全PGL飽和親の非所属へは拡張しない。
reduciblequarticの別全family認識もこの補題だけでは与えない。
