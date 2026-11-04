# κ・追加tail2の固定全QQ係数計算

2026-11-04。凍結preflight73d433...の一列だけを計算する。全四signed cofactorと全四source
columnを保持した `B*kappa*Mcol2` であり、対象はprojected prePOWER15×22である。
元raw H7のeta/HZ/HW/Cdelta liftは別依存として未展開のまま保存する。

全240typed/238live nodeをQQ/a-Laurentで実行し、四cofactor、四B-cleared column、最終
polynomialの全係数をstream保存する。各演算のlast-use、全preloaded atoms10456、scratch、
15000work reserveを数え、予測119119、上限150000を維持する。

核心DAGの値を解放した後、保存された答えから四つのcofactor×column productを逆順で
引き、正確な零を別の加算順で確認する。入力cofactor/columnは積の直後に解放する。
全部分差は同じ21740格子内であり、再乗算確認の上界112416も予測119119以下となる。
これらは新しいcolumn queryではない。

κA=0は実cofactorを定義する普遍的4×3 determinant identityを、別の小さい整数多項式
エンジンで全三列確認して保存する。既知六pure列への零性は、独立受理済み4×6全90二次
minor零、非零anchor polynomial、domain cancellationで全点恒等式になる。
anchorを特殊点で反転する条件は入れない。実四cofactorと列の全係数も独立peerへ渡す。

結果が非零でも新しいsource fragmentの非零性だけである。κ係数のchart単位性、全pure
rank、full cubic57、原rawH7への同cokernel、finite bad locusを主張しない。零やguard停止も
全source付きで保存し、別列・prime・再試行へ自動的に進まない。
