# 全 G3 一次写像の固定 rank・kernel 証明の独立監査

元の source gate とその独立監査で、418 source のうち後半180成分には固定射影が対角 ±1 として作用し、前半238成分には零として作用する。したがって全定数 kernel は前半238の active 写像（252行）の kernel と正確に等しい。

本監査は新たな消去法・rank 探索を行わない。提供された62本の元 source 表示を別の Fraction 算術で行列へ掛け、正規化された像を全て照合した。像の相異なる先頭行と対角1から独立性が従う。元 source 表示は選択列順で三角、対角は各固定 pivot の逆数なので、元62小行列式は非零 pivot 積に等しい。

残る176本についても行列積が厳密に零であること、各 free source 座標が1で他の free 座標が零であることを検証した。選択62列と free176列の分割が全238列を覆うため、これらが完全な kernel 基底となる。active rank62、全 rank242、全 kernel 次元176を受理する。

元 producer・source・note・evidence のハッシュを固定した。producer は import せず、保存済み有限係数の後検証のみを行った。G3二次係数、T2の実HF2、全長、H0指数については結論しない。内部 unit transport の二次式は別 gate である。
