# 残九源ambient batchのB6整数化訂正とpartial receiptの扱い

2026-11-08。元nine caller／job／points／false reportのbytesは変更しない。

最初のONE九batchではM1,M11,B1,B2,B3,B4,B5の七点について全計算が完了した。B6の整数係数guardで停止し、B6のJacobi／ambient及びV1はまだ実行していなかった。B6ではd=−3/4、Pxyのr係数は−1/4、r=y−x²/2なので、tensor scale十二ではx²係数が3/2となる。整数のみを受ける既helperへ渡す**前**にguardがこの点を検出した。

数学的source条件の違反、別のprojective point又はpure25次元入口の反例ではない。rootの明示許可により、新しいcontinuation identityで**同じB6の全tensorをさらに二倍（scale二十四）**し、元と同じV1 scale十二だけを実行した。元七点を再計算せず、別点も選んでいない。

continuationの全Jacobi／全150／350／126／D,V,M／minors／kernelsは完了。B6はD100、V1はD99だった。元partial false receiptを後からPASSへ変更したのではなく、元七件の`arithmetic_complete=true`と新二件を別々に参照するcomplete-point indexを保存する。

元partial batchは11.516864秒、新二continuationは2.983330秒。第一reportのoverall false、B6未完record、例外全文、全元sourceをそのまま残す。完全に終わった七件のmatrix証明書を、整数化実装の後続失敗だけを理由に捨てたり、未実行の全九PASSと呼んだりしない。
