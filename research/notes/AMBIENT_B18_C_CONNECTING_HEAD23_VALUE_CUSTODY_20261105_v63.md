# head23 cache 値 worker の未実行 custody

固定 worker は root38648 の六係数を last-use で合成し、既存五係数と逆元両積を
照合した phase1 を先に保存する。その後、保存済み全 coherent polynomial を
全項評価し、二高座標・canonical座標と五 head の低二行を保存する。
root 段階83bit、post段階207bitという正の事前上界に対し、parse/add/mul/逆元
diagnostic 全てで Fraction 演算前に検査する。hard2048bitは変更しない。

code は stdlib のみで importでき、--help と source値を触らない guard/B6逆元
fixture を確認した。実 value folder は存在しない。独立した全 coherent SOURCE
および single-root/code/whole-budget peer と、親の ONE GO が揃うまで実行しない。
postが停止した場合に既計算 root の phase1 を保存し、再計算を自動では行わない。

単純 target の h0/h1 minor が非零の場合だけ ker(M7→G4)=J2⊕k と判定する。
他の head の高い座標、元 H0、K torsion の指数、nil 三乗は計算対象に含めない。
