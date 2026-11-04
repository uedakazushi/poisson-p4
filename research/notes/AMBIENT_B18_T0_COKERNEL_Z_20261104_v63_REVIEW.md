# Actual有限余核Zの独立監査

対象は `Z=coker(D0→D4)` である。`T0=ker(D0→D4)` と異なるmoduleであり、Zの有限性や長さからT0の有限性・socle・元交叉環の深さを推測しない。

## 両側annihilatorと正当なbase change

CAのdegree-zero coefficient環は `C0_A=R0+I0_A τ_A`、正次数のcoefficient ambientは `H_A=R0+m9 τ_A` である。受理済みのT1消滅により、degree-oneからstable coefficient module D4への射は同型である。従ってZは `H_A/C0_A` のquotientとなる。`f∈I0_A` と `h=e+jτ_A∈H_A` に対し、`fh=fe+fjτ_A∈C0_A` だから、全次数の恒等式として `I0_A Z=0` である。

同じ説明をCBのactual homogeneous presentationに適用する。degree-zero imageは共通のphysical target D0からD4への同じ像であり、degree-oneは同じintrinsic D1である。全CB kernel99とgraph差87を保持しているので、CBでも `I0_B Z=0` となる。この証明はτ_Aとτ_Bを同一視しない。

受理済みのactual等式 `R0/(I0_A+I0_B)=QQ[[r,s]]/(r,s)^2=:A3` により、ZはA3-moduleである。完全なD4 honest9表示にC0の九像列を追加すれば、全余核表示は `[Q1|C0map]` となる。余核のright exactなtensorにより、この全表示をA3へ移すことは正当であり、結果は元のZ自身である。kernelやambient imageのtensorでTorを落とす手順ではない。

C0→D4の原点rank7とD4のminimal9からμ(Z)=2が従い、A3の長さ3により `length(Z)≤6` となる。ここまでのannihilatorとbase changeはfirstjetの観測ではなく、全次数のmodule恒等式である。

## 独立算術

専用checkerは著者producerをimport・実行しない。既に独立受理された二つのpure57 source二次証明書を固定し、疎dual numberで元sourceの一次値を再評価する。全CA792・CB784のprefix quotient一次値を、保存二次単項式のprefix制限差から独自に復元する。新CA35 tail quotientも同じ方法で復元する。旧CBの四つの未指定g微分は独立な形式方向のまま残し、全map成分から厳密に消えることを確認する。新CAの四つのg商はmapのdependency closureに入らない。

Q1の全1548列とC0の九列を、完全な九変数firstjetからA3へ制限する。A3上では高次項は消えるので、これは全次数で正確な係数である。27次元QQ空間に全列とそのr倍・s倍を入れ、著者と逆順のpivotによる独立消去で関係部分空間を求める。著者の保存基底は部分空間の包含と次元から照合し、保存r/s作用は各basis vectorの像を元27次元空間へ戻した合同式で検証する。可換性、最大idealの平方零、および同じphysical R0九座標の作用を確認する。

実際の長さ、関係rank、実測時間・メモリ、現在の依存hashはreview reportに記録する。新しいsource IFT solve・jet・GB・標本は使わない。

## socleへの接続の範囲

D4のdepth1から `H0_m(D4)=0` であり、その部分moduleであるD0の像にも有限台torsionは無い。従って `0→T0→D0→image→0` の左完全性により `H0_m(T0)=H0_m(D0)` は正しい。しかし、この等式も今回計算したZも、H0の値そのものを決めない。残るT0/socleとfull Dの深さは別の計算課題として保持する。

## Canonical moduleと、別のTor論証によるT0非零

計算されたZは長さ3・二生成で、基底u,v,eに対するr/s作用は、二つのgeneratorからsocleへ向かう行列

```
[[8/5,-2/5],[-2/5,1]]
```

で与えられる。行列式36/25が非零なので、ZはA3のQQ双対、すなわちcanonical module `omega_A3` に正確に同型である。この同型の全作用を確認した。R0上のminimal resolutionは、二変数でのcanonical Hilbert--Burch resolutionと七zero座標のKoszul resolutionのtensorであり、Betti多項式は `(2+3t+t²)(1+t)^7`、特にβ1=17である。17本のminimal first relationは全て一次であり、これはgraded Hilbert関数だけからの推測ではない。

C0の七列がD4 free9のE0,...,E6そのものだから、D4→Zのfree liftは残るE7,E8への正確なprojectionである。全Jinf relationのfull9 firstjetをこの二成分へ射影したrankは独立に16だった。Q1だけでも同じrankで、全Jinf firstjetとのspanの一致を確認した。Zのminimal first relationが全て一次であることから、このrankは `Tor1(D4,k)→Tor1(Z,k)` のrankそのものである。高次relationやm倍は、このTor空間では寄与しない。

`B0=image(D0→D4)` とすると、Tor完全列から

```
mu(B0)=17-16+9-2=8.
```

既受理のmu(D0)=9との相違により、**T0は非零**と分かる。これはZが有限であることやC0 mapの原点rank7だけからkernelを推論したものではない。D0のresidue fiberへ見えるT0の像は一つのlineであり、basisのg_A*tau classに一致する。もう一つの原点kernel vector `f60*tau+(1/2)f67*tau-f72*tau` のboundary classは、上のrank16の外側にあることも照合した。raw `g_A*tau` そのものがactual kernel elementであるとは主張しない。

T0の非零性は、有限台socleの非零性と異なる。T0の有限性、cyclicity、socle、元交叉環の深さはなお未決定である。
