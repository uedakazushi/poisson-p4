# K92：global cleanupは既存残余四角の一行倍である

既存の三つのSchur pivotはT7の行列順P=[0,6,1]、残余はJ=[2,3,4,5]である。generic cleanupはT7の行2をdelta*(行2−59375 a^{-68}行1)にする。行2は残余側、行1は既存pivot側なので、P3自体は変わらない。pivot行を引いた部分のSchurは零であり、新残余はdiag(delta,1,1,1)倍の同じ残余四角である。従って既存七pivotを保ち、別の五pivotからやり直さない。

Q2の列1はe0、列6はe6である。N0とN6はいずれも既に消去したpivot列の成分であり、新しい残余方向の階数増加を意味しない。正しいP順でのpivot row callbackは[N2,N6,a^{-68}N0]。残余側は[N1,N3+18N2/7,N4+147N2/25,N5+126N2/25]である。これらの全pivot補正を保持する。

N6の同じformal代入値は零だったが、generic rootを零へ置き換えない。次の意味のある計算対象は、四角残余の実一行またはそのborderであり、安いcallbackを順に求めることではない。新しいcoefficient jobは実行せず、全chart補集合・全域57・行列式の未決を保つ。
