# κ17 独立監査の packed Laurent source gate

旧checkpointを開かず、凍結V2 literal20 atomと49定義、既受理四κを用いる。producer/nativeはimportしない。各Laurent単項式をa+64e+4096dで符号付き整数keyにし、nodeごとに固定共通分母と整数係数辞書で表す。全source/product boxのa幅が64未満、e範囲も64未満なので、この符号化は単射であり負aにも使える。

全49 sourceのbox・次数・support・共通分母・絶対L1 boundを別の正の算術で再生する。四積は別のproduct bufferを作らず、一つの共通分母を持つaccumulatorへ直接加える。sourceはlast-useで解放し、κは一つずつ読む。raw literal records、packed keys、係数、分母、headers、root照合の一時JSON、metadata/整数temporaryを全て計上する。

符号付きsourceやdot値を求める前に、旧75410 slot・548bit capを上げずに収まることをreportで確認する。新しい物理source、point、prime、rank探索はない。実workerは全49 sourceと四source rootのQQ一致を確認した後、固定順序のκ内積一つだけを検査する。κ17零から言えるのは列17のx4零のみである。
