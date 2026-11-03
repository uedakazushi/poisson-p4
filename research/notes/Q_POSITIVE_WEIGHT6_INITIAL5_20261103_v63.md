# 正 weight6 軌道の完全な五次 initial と残る三平面

2026-11-03。実八変数形式模型の**全 ordinary initial ideal の次数5部分**を、
すべての低次数 cancellation を含めて厳密 QQ で求めた。
二次5本・三次8本・四次8本に対し、新しい五次15本で全次数5を生成する。
実横断局所環の Hilbert function は

    (1,8,31,72,126,171)

である。全五次までの斉次 ideal は次元3、次数5であり、radical は

    (z1,z2,z4,z6,z7) ∩ (z1,z3,z4,z5,z6,z7)

である。前者は自由座標(z0,z3,z5)の三平面、後者は自由座標(z0,z2)の二平面である。
実形式横断環の次元上界3、全 projective Poisson 局所次元上界25は改善されなかった。
三平面の実形式的持ち上げ・新しい25次元成分・次元の等号はまだ証明していない。
次元24以下も未証明である。

## 全350式の五次形式解

前の凍結三次・四次模型の同じ95単元 block と八方向を使い、全150 tensor の五次解を求めた。
Jacobi は物理 tensor に関して二次なので、五次の未補正項は

    Bilinear(R1,R4)+Bilinear(R2,R3)+Bilinear(R3,R2)+Bilinear(R4,R1)

の全係数である。95選択式の係数へ −B^(-1) を掛け、元の全350式へ戻すことで五次残余を得た。
選択95行は全係数で零になる。全350五次残余と95形式補正を中間ファイルへ先に保存した。
一つの178秒・2.2GiB制限の計算は115.57秒で完了し、保存を終えた後に制限へ達したものではない。

## 全 cancellation を小さな48関係に還元

二次独立な五つの actual 形式生成元を A、三次の新しい八つを C とする。
degree3 source (S1*A2,C3) の48行はすべて独立である。
したがって元の345本 G の全三次係数をこの48行で一意に表し、同じ polynomial 倍の actual A・C を引くと、
残る337本 H は次数4以上から始まる。これは ideal を保存する三角的な形式 row 操作である。

H の四次部分から八つの新しい D を選ぶ。
最小 degree4 source (S2*A2,S1*C3,D4) は252行・rank204であり、全 QQ kernel は48次元である。
全337本 H4を204本の独立な行で表し、同じ actual 形式級数の polynomial 倍を引くと、
Dを除いた329本 E は次数5以上から始まる。この row 操作も全形式 ideal を保つ。
全48 kernel と、これらの表示係数の全係数恒等式を保存して検査した。

これにより全 degree5 initial 空間は正確に

    S3*A2 + S2*C3 + S1*D4 + span(E5)
      + carry5(kernel(S2*A2,S1*C3,D4))

となる。最初の三項が952行、E5が329行、kernel carry が48行で、計1329行・792列である。
degree2で消えるには A の定数係数が零、degree3で消えるには上の48行の独立性により
A の一次係数と C の定数係数が零となる。
次の cancellation は252行の kernel48に尽くされるため、ほかの degree5 source はない。
これは元の350形式生成元の任意の形式級数係数による組合せを次数ごとに扱った論証であり、
各行の最低次数だけを集めた計算ではない。

## QQ での完全性と radical

有限体32003で rank621の証人を取り、新しい五次候補15本を選んだ。
この有限体 rank を QQ の上界として使っていない。
その15本と前の21本で作る小さな QQ 八変数 ideal の GB を求め、
**全329+48=377本の厳密 QQ 五次式がその GB で零へ還元されること**を追加検証した。
逆包含は15本が実際の1329行の選択であることから従う。
よって生成する次数5の空間は正確に同じであり、QQ Hilbert function の値171から rank621も厳密に従う。

既証の二次 radical は二つの五平面に尽くされる。
一方(z1=z5=z7=0)への制限は五変数(u0,u2,u3,u4,u6)を使い、
その minimal primes は (u6,u4,u3) と (u6,u4,u2) である。
他方(z1=z7=0,z4=−2z6)への制限は五変数(v0,v2,v3,v5,v6)を使い、
minimal primes は (v6,v2) と (v6,v5,v3) である。
全て QQ 上の線形素 ideal なので、複素数への拡大でも素性は保たれる。
重複と包含された平面を除くと、冒頭の三平面と二平面だけが残る。
全377式の QQ 還元とこの radical 分解の追加計算は約3.8秒で終わった。

## 保存と再開

`code/q_positive_weight6_initial5_20261103_v63.py` が全五次と cancellation を構成する。
native modular rank helper は同名 `_rank_` C++、厳密 QQ completion・radical split は同名 `_cone_` Python、
保存 manifest は同名 `_publish_` Python である。
`reports/Q_POSITIVE_WEIGHT6_INITIAL5_20261103_v63/` に実 construction receipt、全350五次、全1329行、
48 kernel・329 E の cofactor 表示、選択15本、QQ還元・radicalの M2入力と成功ログを保存する。

次は残る三平面の actual 形式持ち上げ、またはその一般方向を除く構造的 obstruction を調べる。
ordinary initial の三平面だけから実三次元族を推論してはいけない。
六次 jet や別の大規模 GB は本 packet では実行していない。
