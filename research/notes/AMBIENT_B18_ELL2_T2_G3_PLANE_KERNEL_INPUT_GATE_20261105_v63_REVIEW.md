# G3 plane-kernel 拡大の source gate 独立監査

既知 L48/ell² cover の固定 normal form を独立再帰し、六 target・七 z 座標・九変数の全378積 `e_t z_i x_j` が零となることを確認した。これは次数二初期加群の事実であり、元 T2 の全次数で z_i mT2 が零という主張ではない。

238 active source の実42 target一次写像から、r/ell の全84係数と、定数商 PiY による z の21係数を独立構成した。省略された projected r/ell 六行は最初84行の明示定数結合である。従ってこの105行条件の kernel は、実 minimal6 で一次零となる十分な family を与える。全238列を照合したが rank は求めていない。既受理 full active rank62 の線形射影なので rank は62以下である。

元 G3 列の bottom30 は全次数で零である。最初27消去はこの列を変更しない。残る99消去を含む target 写像を F(R) とすると、定数部分は `(PiY0;0)`、未知一次部分 F1 には商射影と bottom 単元の微分が含まれる。条件内の c について G3_1 c は z の一次項だけを持ち、その定数射影は零である。従って二次項の未知寄与 `F1 G3_1 c` は全て `e_t z_i x_j` であり、上記378恒等式から既知二次 cover 内で消える。未知微分を零としたのではない。受理済み constant coherent graph frame もこの説明と両立する。

これは過去の unrestricted projected393 kernel の再採用ではない。r/ell の **full**84条件が追加されている。新しい rank/kernel worker は75684 numeric slots、raw18167424 bytes の source bound 内で、384-bit fail-closed guard を保つ必要がある。この高さは予測保証ではない。

後続二次 source gate では q=D31 I3 c の z 射影を保持する。一般 q 自体が固定182列の span に入るとは仮定せず、`q−Efree42 Pi42 q` を入れる。これは K0 の最後42列が Efree42 なので、選択126座標の抽出を変更しない。現 gate は二次値、実HF2、全長を決めない。
