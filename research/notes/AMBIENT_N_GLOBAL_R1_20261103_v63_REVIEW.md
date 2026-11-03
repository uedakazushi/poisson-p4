# 全射影 N の R1：全旗ファイバーを使う独立監査

2026年11月3日。frozen QQ producer を別directoryで再実行し、秒数以外一致した。
独立の sparse150物理計算で全N19係数と全GL5の微分行列、非零28次minor6912/390625、
intrinsic factor行列rank4・核L=∂yを作り直した。
canonical field Ctilde=(h²,−4hx,−15w²+hz,15z²+hw)と旗条件の全係数恒等式も再計算した。

全幾何的旗ファイバーが一つである証明を受理する。
点はcovector Casimirを持たずNfへ入らないので、全proper模型の持上げはNgである。
unique intrinsic factorがLを戻す。自己cofactor条件でlambda0は四成分の独立性に反し、
lambda非零ならx²,zw,z²,w²係数がB=C=D=0を強制する。
従ってH=Ch。第二covector条件はk∈span(h,x)を強制し、独立性からW=span(h,x)。
固定旗でのcanonical係数写像は単射なので、係数点も一意である。

源は滑らかなprojective27次元で、全affine像rank28からprojective微分は単射となる。
これはtargetの接空間を源の像と同一視したものではなく、**fiber tangentが零**を意味する。
全fiberの台が一点なので、有限型Artin fiberのm/m²=0とNakayamaからscheme fiberは長さ1。
propernessから近傍の射は有限で、被約閉像Nの局所環A→有限Bは単射である。
B/mA B=Cから再びNakayamaでB=Aとなり、その点でtarget自身が滑らかと分かる。

Ngは既約26次元でNの唯一のrank≤2境界因子、Nfは25次元でcodim2。
従ってNgのgeneric pointは滑らかであり、受理済みrank4開部分のR1と合わせ、全射影NがR1となる。
この結論は従来のS61範囲だけのR1をここで正当に拡張する。
正規性やCM性は推論せず、指定rank4点Qの非正規・非CM性とも両立する。
再現コードと結果は ambient_N_global_R1_review_20261103_v63.py、同名 REVIEW report に保存する。
